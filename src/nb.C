/**************************************************************************************
 * Workhorse of the program. This file contains functions to perform Naive Bayes
 * Haplotype Calling of individual reads and to assign genotype scores using the
 * haplotype calls.
 * 
 * Haplotypes are called only for reads that contain at least one overlapping known snp
 * with previous read.
 * 
 * Genotypes are assigned separately for candidates that lie on reads that were assigned
 * a haplotype or not:
 *  1. For reads that are assigned haplotypes: genotype score includes haplotype
 *     probability in calculations. Genotype scores are assigned only if >90% of the
 *     reads that contain the site were assigned haplotypes, such that these reads
 *     are also contiguous having no gaps of unassigned reads.
 *  2. For reads that are not assigned haplotypes: genotype scores do not include
 *     haplotype probability. Genotype scores are assigned to all sites that are
 *     covered by a read.
 * ***********************************************************************************/
 

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <assert.h>
#include <math.h>
//#include <boost/math/distributions/skew_normal.hpp>
//#include <boost/math/distributions/binomial.hpp>
//#include <boost/math/distributions/beta.hpp>

using namespace std;

#include "read.h"
#include "snp.h"
#include "nb.h"

extern vector<SNP*> snp_list;
extern vector<READ*> reads_list;
extern double errate;

int PROXIMAL_READ_CT = 10;

//#define DEBUG

// Works on known snp list
// Given t-th read in list, returns base qual of snp with position pos.
int GetPosQual(int t, long pos)
{
	READ *rd = reads_list[t];
	int snp_count = rd->GetKnownCount();
	SNP ** snp_list = rd->GetKnownList();
	for(int i=0; i<snp_count; i++) {
		if(snp_list[i]->GetPos()==pos)
			return rd->GetKnownQualScore(i);
	}
	return 'N';
}

// Works on known snp list
// Given t-th read in list, returns allele of snp with position pos.
char GetPosAllele(int t, long pos)
{
	READ *rd = reads_list[t];
	int snp_count = rd->GetKnownCount();
	SNP ** snp_list = rd->GetKnownList();
	for(int i=0; i<snp_count; i++) {
		if(snp_list[i]->GetPos()==pos)
			return rd->GetKnownAllele(i);
	}
	return 'N';
}

void NaiveBayes(int snp_start, int snp_end, int T)
{
	int i, j, t;
        SNP **known_snp_list = new SNP*[100];
        int *known_index = new int[200];
	int reads_list_size = reads_list.size();

	for (t = 1+reads_list_size-T; t <= reads_list_size-1; t++) {
		int known_snp_count = 0, flag = 0, direction = 1, hap = 0, mismatchsum = 0;
		double norm = 0.0;
                double emission_list[3], haprob[3], emission_set[3][100];
		READ *pd = reads_list[t-1];
		READ *nd = reads_list[t];
		int plen = pd->GetLen();
		int nlen = nd->GetLen();
		int ppos = pd->GetPos();
		int npos = nd->GetPos();
                int rd_start = ppos < npos ? npos : ppos;
                int rd_end = ppos+plen > npos+nlen ? npos+nlen-1 : ppos+plen-1;
		
		//mismatchsum = GetMismatchSum(nd);
		//if(mismatchsum>=50)
		//	continue;
		
                GetKnownSnpList(known_snp_list, &known_snp_count, known_index, t);
                // This adds the count of known overlapping snps between adjacent reads
		(*nd).AddKnownCount(known_snp_count);
#ifdef DEBUG
if(known_snp_count>0)
cout << "Read " << nd->GetPos() << " ; "<< t << " and " << t+1 << ", " << rd_start << ", " << rd_end << " with " << known_snp_count << " known snps." << endl << endl;
cout << "SnpPos  R A LL..\tAA_gen AA_obs AA_genob\tAB_gen AB_obs AB_genob\tBB_gen BB_obs BB_genob\tprob\n";
cout << "Throwing emissions from read: " << nd->GetPos() << endl;
#endif

		for(j=1;j<=2;j++) {
#ifdef DEBUG
cout << "Haplotype " << j << endl;
#endif
			emission_list[j] = 0.0;
			for(int count=0; count<known_snp_count; count++) {
				double emission;
                       	        SNP *sp = known_snp_list[count];
#ifdef DEBUG
cout << sp->GetPos() << " " << sp->GetRef() << " " << sp->GetAlt();
#endif
				// Obtains naive bayes score for the count-th overlapping known snp 
				// between the adjacent reads for haplotype, 'j'.
				// Working with relative haplotype here.
				emission = compute_new_emission(known_snp_list, count, t, known_index, j);
					if(emission==0.0)
						continue;
					emission_list[j] += ((emission));
					emission_set[j][count] = emission;
			}
#ifdef DEBUG
cout << "Total Emission = " << emission_list[j] << endl;
#endif
	  	}
		if(emission_set[1][0]<emission_set[2][0])
			direction = 2;
		for(int count=1; count<known_snp_count; count++) {
			if(emission_set[1][count]>emission_set[2][count]&&direction==2)
				flag = 1;
			if(emission_set[1][count]<emission_set[2][count]&&direction==1)
				flag = 1;
		}
#ifdef DEBUG
cout << "Discordance = " << flag << endl;
#endif
		norm = emission_list[1] + emission_list[2];
		if(norm!=0.0) {
			haprob[1] = emission_list[1] / norm;
			haprob[2] = emission_list[2] / norm;
		} else {
			haprob[1] = 0.5;
			haprob[2] = 0.5;
		}
#ifdef DEBUG
cout << "Happrob[1] = " << haprob[1] << endl;
cout << "Happrob[2] = " << haprob[2] << endl;
#endif
		// If haplotype probabilities are equal, randomly assign haplotype
		hap = haprob[1] > haprob[2] ? 1 : haprob[1] == haprob[2] ? t%2+1 : 2;
		double happrob = haprob[hap];
		// convert relative haplotype to absolute haplotype
		hap = (*pd).GetHap() == hap ? 1 : 2;
		nd->assignHaplotype(hap, happrob, flag);
		flag = 0; direction = 1;
	}

	delete [] known_index;
	delete [] known_snp_list;
}

// Assign genotype/somatic scores
void FindSomaticMutations(long snp_start, long snp_end)
{
	int known_hap_count = 0, gap;
	for(vector<SNP*>::iterator snp_it = (snp_list).begin(); snp_it != (snp_list).end(); snp_it++) {
		if((*snp_it)->GetKnown())
			continue;
		long pos = (*snp_it)->GetPos();
		if(pos<=snp_start)
			continue;
		else if(pos>snp_end)
			break;

		int g=0;
		double post[3], prior[3], happrob[3];
		double norm = 0.0, snp_rate = 0.000003;

#ifdef DEBUG
cout << "Prior het prob for som " << (*snp_it)->GetPos() << " = " << prior[1] << endl;
#endif
		// Get posterior probability from haplotype probabilities
		gap = somaticHaplotypeProbability(snp_it, happrob, &known_hap_count);

		prior[0] = 1 - snp_rate - snp_rate*snp_rate;
		prior[1] = snp_rate;
		prior[2] = snp_rate*snp_rate;

		for(g=0; g<3; g++) {
			norm += prior[g]*happrob[g];
		}

		for(g=0; g<3; g++) {
			post[g] = (prior[g]*happrob[g])/norm;
		}

/*
		Should not run
		known_hap_count=1;
		if(known_hap_count==0||gap==4) {
			post[0] = (*snp_it)->GetPosteriors()[0];
			post[1] = (*snp_it)->GetPosteriors()[1];
			post[2] = (*snp_it)->GetPosteriors()[2];
#ifdef DEBUG
cout << "Zero known_snp_count for " << (*snp_it)->GetPos() << endl;
#endif
		}
*/
		if(post[1] > 1.0)
			cout << "Posterior het for somatic " << (*snp_it)->GetPos() << " > 1.0" << endl;
    		(*snp_it)->add_somatic_posteriors(post);

#ifdef DEBUG
cout << "Posterior het prob for som " << (*snp_it)->GetPos() << " = " << post[1] << endl;
#endif
	}
}

// Calculate posterior probability
int somaticHaplotypeProbability(vector<SNP*>::iterator snp_it, double probs[3], int *known_hap_count)
{
        char ref = (*snp_it)->GetRef();
        char alt = (*snp_it)->GetAlt();
        int real_count = 0, ref_ct = 0, alt_ct = 0;
	int ref1_ct = 0, alt1_ct = 0, ref2_ct = 0, alt2_ct = 0, err_ct=0;
	int contig_max = 0, contig_new = 1, max_stretch = 0;
	int count = 0, gap = 0, stretch = 0, contiguous = 1, contig_start = 0;
	int mapq1 = 0, proximal_ins_sum = 0, proximal_del_sum = 0;
	int clus_flag = 0, clus_in_flag = 1, clus_del_flag = 1;
	int max_reads = (*snp_it)->GetReadCount();
	double probs1[3], probs2[3];
        double prevalence;

	ref_ct = (*snp_it)->GetRefCount();
	alt_ct = (*snp_it)->GetAltCount();
	
	// dynamically calculate the prevalence
	prevalence = (double)(alt_ct)/(double)(ref_ct+alt_ct);
	*known_hap_count = 0;
        probs[0] = 1.0; probs[1] = 1.0; probs[2] = 1.0;
        probs1[0] = 1.0; probs1[1] = 1.0; probs1[2] = 1.0;
        probs2[0] = 1.0; probs2[1] = 1.0; probs2[2] = 1.0;

#ifdef DEBUG
cout << " Count on snp " << (*snp_it)->GetPos() << " is " << max_reads << endl;
#endif

	// First for loop looks for maximum contiguous region
        for(count=0; count<max_reads; count++) {
                READ *rd = (*snp_it)->GetRead(count);
                READ *pd = (*snp_it)->GetRead(count-1);

                if(rd->GetKnownOverlapCount()>0&&rd->GetDiscordance()==0) { // No gap exists in hap assignment
			(*known_hap_count)++;
			stretch++;
			if(contig_new==1) {
				contig_start = count;
				contig_new = 0;
			}
                } else { // Gap exists in hap assignment. Not contiguous.
			contiguous = 0;
			contig_new = 1;
			if(stretch>max_stretch) { // Identify longest stretch of contiguous assignments
				max_stretch = stretch;
				contig_max = contig_start;
			}
			stretch = 0;
		}
		if(contiguous==1) {
			max_stretch=stretch;
			contig_max=0;
		} else if(stretch>max_stretch) {
			max_stretch = stretch;
			contig_max = contig_start;
		}
#ifdef DEBUG
cout << "SNP=" << (*snp_it)->GetPos() << " Gap=" << gap << " contiguous=" << contiguous << " max_reads=" << max_reads << " max_stretch=" << max_stretch << " contig_max=" << contig_max << endl;
#endif
	}

	if(contiguous==0)
		gap+=2;
	if(max_stretch<(int)(9*max_reads/10)) // If max stretch of contiguous assigned reads is at least 90% of total number, consider reads for genotype score assignments.
		gap++;

#ifdef DEBUG
cout << " Gap=" << gap << " contiguous=" << contiguous << " max_reads=" << max_reads << " max_stretch=" << max_stretch << " contig_max=" << contig_max << " prevalence = " << prevalence << endl;
#endif

	// Second for loop performs actual genotype score calculation
        for(count=0; count<max_reads; count++) {
		bool proximal_ins = 0, proximal_del = 0;
                char all, pall = 'N';
		int flag = 0;
                double qual, qualscore;
                READ *rd = (*snp_it)->GetRead(count);
                READ *pd = (*snp_it)->GetRead(count-1);

                int hap = rd->GetHap();
                double haprob = rd->GetHapProb();
		double hethap = sqrt((1.0-haprob)*haprob);
		int snp_count = rd->GetSnpCount();

                if(snp_count==0) {
#ifdef DEBUG
cout << "For som " << (*snp_it)->GetPos() << " skipping read " << rd->GetPos() << " for happrob calculation" << endl;
#endif
			haprob = 0.5;
                }

                // Obtain allele, base qual and #proximal insertions/deletions on current read
		long snp_position = (*snp_it)->GetPos();
                for(int s_pos=0; s_pos<snp_count; s_pos++) {
                        if(rd->GetSnp(s_pos)->GetPos() == snp_position) {
				proximal_ins = rd->GetProximalInsert(s_pos);
				proximal_del = rd->GetProximalDelete(s_pos);
                                all = rd->GetAllele(s_pos);
                                qual = rd->GetQualScore(s_pos)-33;
                                qualscore = pow(10.0,-(qual/10.0));
                                break;
                        }
                }
		if(hap==1) { // See comments later
                	if(all == ref) {
                        	ref1_ct++;
                	} else if(all == alt) {
                                alt1_ct++;
			} else 
				err_ct++;
                } else if(hap==2) {
                	if(all == ref) {
                        	ref2_ct++;
                	} else if(all == alt) {
                                alt2_ct++;
			} else 
				err_ct++;
		}
                real_count++;

		if(qual<5) // No base with quality < 5 should be considered
			continue;
		if(all==alt && qual>=20) // At least one alternate allele with base quality >20
			mapq1++;
		if(proximal_ins==1) { // Read contains insertion in vicinity
			proximal_ins_sum+= 1;
		}
		if(proximal_del==1) { // Read contains deletion in vicinity
			proximal_del_sum+= 1;
		}
		if(all==alt) {
			clus_in_flag--;
			clus_del_flag--;
			if(proximal_ins==1)
				clus_in_flag++;
			if(proximal_del==1)
				clus_del_flag++;
	    	}
		proximal_ins = 0;
		proximal_del = 0;

		// Work on contiguous stretch only
		if(count>=contig_max&&count<contig_max+max_stretch) {
			// Since we are haplotype-aware, two separate set of calculations are performed
			// First set assumes somatic mutation lies on haplotype 1
			// Second set assumes somatic mutation lies on haplotype 2
			// In the end we consider that haplotype to contain somatic mutation, which has higher het probability assigned
			if(hap==1) {
                		if(all == ref) {
                	        	probs1[0] *= ((haprob)*(1.0-qualscore));
                	        	probs1[1] *= ((1.0-prevalence-qualscore/4.0)*haprob);
                	        	probs1[2] *= ((1.0-haprob) * (qualscore/3.0));
                	        	probs2[0] *= ((haprob)*(1.0-qualscore));
                	        	probs2[1] *= ((1.0-prevalence-qualscore/4.0)*(1.0-haprob));
                	        	probs2[2] *= ((1.0-haprob) * (qualscore/3.0));
                		} else if(all == alt) {
                	        	probs1[0] *= (haprob*(1.0*qualscore/3.0));
                	        	probs1[1] *= ((1.0*qualscore/3.0)*(1.0-haprob));
                	        	probs1[2] *= ((1.0-haprob)*(1.0*qualscore/3.0));
                	        	probs2[0] *= ((1.0-haprob) * (qualscore/3.0));
                	        	probs2[1] *= (prevalence*haprob);
                	        	probs2[2] *= ((haprob)*(1.0-qualscore));
				}
                	} else if(hap==2) {
                		if(all == ref) {
                	        	probs1[0] *= ((haprob)*(1.0-qualscore));
                	        	probs1[1] *= ((1.0-prevalence-qualscore/4.0)*(1.0-haprob));
                	        	probs1[2] *= ((1.0-haprob) * (qualscore/3.0));
                	        	probs2[0] *= ((haprob)*(1.0-qualscore));
                	        	probs2[1] *= ((1.0-prevalence-qualscore/4.0)*haprob);
                	        	probs2[2] *= ((1.0-haprob) * (qualscore/3.0));
                		} else if(all == alt) {
                	        	probs1[0] *= ((1.0-haprob) * (qualscore/3.0));
                	        	probs1[1] *= (prevalence*haprob);
                	        	probs1[2] *= ((haprob)*(1.0-qualscore));
                	        	probs2[0] *= (haprob*(1.0*qualscore/3.0));
                	        	probs2[1] *= ((1.0*qualscore/3.0)*(1.0-haprob));
                	        	probs2[2] *= ((1.0-haprob)*(1.0*qualscore/3.0));
				}
			}
		} else if(rd->GetKnownOverlapCount()>=0) { // Insufficient contiguous regions. Haplotype unaware
                	if(all == ref) {
                        	probs1[0] *= ((1.0-qualscore));
                        	probs1[1] *= (1.0-prevalence-qualscore/3.0);
                        	probs1[2] *= ((qualscore/3.0));
                        	probs2[0] *= ((1.0-qualscore));
                        	probs2[1] *= (1.0-prevalence-qualscore/3.0);
                        	probs2[2] *= ((qualscore/3.0));
                	} else if(all == alt) {
                        	probs1[0] *= ((qualscore/3.0));
                        	probs1[1] *= (prevalence-qualscore/3.0);
                        	probs1[2] *= ((1.0-qualscore));
                        	probs2[0] *= ((qualscore/3.0));
                        	probs2[1] *= (prevalence-qualscore/3.0);
                        	probs2[2] *= ((1.0-qualscore));
			} else {
                        	probs1[0] *= ((2.0*qualscore/3.0));
                        	probs1[1] *= (qualscore/2.0);
                        	probs1[2] *= ((2.0*qualscore/3.0));
                        	probs2[0] *= ((2.0*qualscore/3.0));
                        	probs2[1] *= (qualscore/2.0);
                        	probs2[2] *= ((2.0*qualscore/3.0));
			}
	    }

#ifdef DEBUG
cout << "For som " << (*snp_it)->GetPos() << " bearing allele " << all << ref << alt << " on read " << rd->GetPos() << ", qualscore: " << qualscore << " has known_hap_count " << rd->GetKnownOverlapCount() << " and flag " << flag << " bearing haplotype " << hap << " with haprob " << haprob << " for happrob calculation " << probs1[0] << ":" << probs1[1] << ":" << probs1[2]  << " and " << probs2[0] << ":" << probs2[1] << ":" << probs2[2] << " and proxdel " << proximal_del << " and proxins " << proximal_ins << " and clus_in " << clus_in_flag << " and clus_del " << clus_del_flag << " and clus " << clus_flag << endl;
#endif
        }

	// Pick one of two sets of scores
	double norm1 = probs1[0] + probs1[1] + probs1[2];
	double norm2 = probs2[0] + probs2[1] + probs2[2];
	if(probs1[1]/norm1>probs2[1]/norm1) {
		probs[0]=probs1[0];probs[1]=probs1[1];probs[2]=probs1[2];
	} else if(probs2[1]/norm2>probs1[1]/norm1) {
		probs[0]=probs2[0];probs[1]=probs2[1];probs[2]=probs2[2];
	} else if(probs1[0]/norm1<probs2[0]/norm2) {
		probs[0]=probs1[0];probs[1]=probs1[1];probs[2]=probs1[2];
	} else if(probs2[0]/norm2<probs1[0]/norm1) {
		probs[0]=probs2[0];probs[1]=probs2[1];probs[2]=probs2[2];
	} else {
		probs[0]=probs1[0];probs[1]=probs1[1];probs[2]=probs1[2];
	}

	PROXIMAL_READ_CT = 2 + (int)ceil((double)(real_count)/10.0);
	if(clus_in_flag>=1||clus_del_flag>=1)
		clus_flag = 1;
	// FILTERS
	if(ref_ct<3 || alt_ct<1 || (gap!=3&&alt1_ct>=2&&alt2_ct>=2) || mapq1==0 || proximal_ins_sum>=PROXIMAL_READ_CT&&alt_ct<5 || proximal_del_sum>=PROXIMAL_READ_CT&&alt_ct<5 || clus_flag==1) { // || (gap==3&&alt_ct>=ref_ct))
		gap=5;
		if(mapq1==0)
			gap = 6;
		if(proximal_ins_sum>=PROXIMAL_READ_CT&&alt_ct<5 || proximal_del_sum>=PROXIMAL_READ_CT&&alt_ct<5)
			gap = 7;
		if(clus_flag==1)
			gap = 8;
		probs[0] = 1.0; probs[1] = probs[2] = 0.0;
	}
        if(probs[0]<0.0||probs[1]<0.0||probs[2]<0.0) {
                cout << "Negative: " << probs[0] << ":" << probs[1] << ":" << probs[2] << endl;
                probs[0] = 1.0; probs[1] = probs[2] = 0.0;
		gap=5;
	}

#ifdef DEBUG
cout << "Accepted mutation " << (*snp_it)->GetPos() << " Prox=" << PROXIMAL_READ_CT << " Gap=" << gap << " known=" << (*snp_it)->GetKnown() << " real_count=" << real_count << " ref_ct=" << ref_ct << " ref1_ct=" << ref1_ct << " ref2_ct=" << ref2_ct << " alt_ct=" << alt_ct << " alt1_ct=" << alt1_ct << " alt2_ct=" << alt2_ct << " proxdelsum=" << proximal_del_sum << " proxinssum=" << proximal_ins_sum << " clus_in=" << clus_in_flag << " clus_del=" << clus_del_flag << " clus=" << clus_flag << endl;
#endif
	return gap;
}

int GetMismatchSum(READ *read)
{
	int read_snp_count = read->GetSnpCount();
	SNP **snp_list = read->GetSnpList();
	int it=0;
	double sum=0.0;

	while(it<read_snp_count) {
		SNP *snp = snp_list[it];
		if(snp->GetKnown()==0) {
			sum += snp->GetQualScore();
		}
		it++;
	}
	return (int)sum;
}

// Obtain list of overlapping known (dbsnp) sites between adjacent reads
void GetKnownSnpList(SNP** known_snp_list, int *known_snp_count, int *known_index, int t)
{
	READ *prev_read = reads_list[t-1];
	int prev_read_snp_count = prev_read->GetSnpCount();
	SNP **prev_snp_list = prev_read->GetSnpList();
	READ *curr_read = reads_list[t];
	int curr_read_snp_count = curr_read->GetSnpCount();
	SNP **curr_snp_list = curr_read->GetSnpList();
	int it1 = 0, it2 = 0, it3 = 0, it4 = 0;

	while(it1<prev_read_snp_count && it2<curr_read_snp_count) {
		SNP* snp1 = prev_snp_list[it1];
		SNP* snp2 = curr_snp_list[it2];
		long snp1_pos = snp1->GetPos();
		long snp2_pos = snp2->GetPos();

		if(snp1_pos == snp2_pos) {
			if(snp1->GetKnown()==2) {
				known_snp_list[*known_snp_count] = snp1;
				known_index[2*(*known_snp_count)] = it1;
				known_index[(2*(*known_snp_count))+1] = it2;
				known_snp_list[*known_snp_count]->IncrKnownOverlapCount();
				(*known_snp_count)++;
			}
			it1++;it2++;
		} else if(snp1_pos < snp2_pos) {
			it1++;
		} else if(snp1_pos > snp2_pos) {
			it2++;
		}
	}
}

// Compute naive bayes probability
double compute_new_emission(SNP **reads_snp_list, int count, int t, int *index, int hap)
{
	double err_rate = errate;
	double known_snp_rate = 0.1;
	double gen_prior[3];
	double snp_rate = known_snp_rate;
	double gen_lik[3] = {0.9, 0.01, 0.09};
	double obs_lik[3] = {0.0, 0.0, 1.0};
	double qual_lik[3] = {1.0, 1.0, 1.0};
	double prob = 0.0;
	double lik = 0.0;
	char dummy = 'N';

	char ref = reads_snp_list[count]->GetRef();
	char alt = reads_snp_list[count]->GetAlt();
	int pos = reads_snp_list[count]->GetPos();
	char all1 = GetPosAllele(t-1,pos);
	char all2 = GetPosAllele(t,pos);
	int quals = GetPosQual(t,pos);
	double qual = pow(10.0,-quals/10);
	int obt = (ref==all1) ? ((ref==all2) ? 1 : 2) : ((ref==all2) ? 3 : 4);

	if(all1==dummy||all2==dummy)
		return 0.0;

	qual_lik[0] = qual/2.0;
	qual_lik[1] = 1.0 - qual;
	qual_lik[2] = qual/2.0;
	gen_prior[0] = 1 - snp_rate - snp_rate*snp_rate;
	gen_prior[1] = snp_rate;
	gen_prior[2] = snp_rate*snp_rate;
	double good = (1-err_rate)*(1-err_rate);
	double bad = (1-err_rate)*err_rate;
	double ugly = err_rate*err_rate;

#ifdef DEBUG
cout << "\t" << all1 << "\t" << all2;
#endif
	for(int j=0; j<3; j++) {
		lik += qual_lik[j] * gen_prior[j];
	}

	for(int i=0; i<3; i++) {
		gen_lik[i] = (qual_lik[i] * gen_prior[i])/lik;

		switch(i) {
		case 0: // genotype = A/A
			if(obt==1&&hap==1)
				obs_lik[i] = good;
			else if(obt==1&&hap==2)
				obs_lik[i] = good;
			else if(obt==2&&hap==1)
				obs_lik[i] = bad;
			else if(obt==2&&hap==2)
				obs_lik[i] = bad;
			else if(obt==3&&hap==1)
				obs_lik[i] = bad;
			else if(obt==3&&hap==2)
				obs_lik[i] = bad;
			else if(obt==4&&hap==1)
				obs_lik[i] = ugly;
			else if(obt==4&&hap==2)
				obs_lik[i] = ugly;
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		case 1: // genotype = A/B
			if(obt==1&&hap==1)
				obs_lik[i] = good + ugly;
			else if(obt==1&&hap==2)
				obs_lik[i] = 2*bad;
			else if(obt==2&&hap==1)
				obs_lik[i] = 2*bad;
			else if(obt==2&&hap==2)
				obs_lik[i] = ugly + good;
			else if(obt==3&&hap==1)
				obs_lik[i] = 2*bad;
			else if(obt==3&&hap==2)
				obs_lik[i] = good + ugly;
			else if(obt==4&&hap==1)
				obs_lik[i] = ugly + good;
			else if(obt==4&&hap==2)
				obs_lik[i] = 2*bad;
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		case 2: // genotype = B/B
			if(obt==1&&hap==1)
				obs_lik[i] = ugly;
			else if(obt==1&&hap==2)
				obs_lik[i] = ugly;
			else if(obt==2&&hap==1)
				obs_lik[i] = bad;
			else if(obt==2&&hap==2)
				obs_lik[i] = bad;
			else if(obt==3&&hap==1)
				obs_lik[i] = bad;
			else if(obt==3&&hap==2)
				obs_lik[i] = bad;
			else if(obt==4&&hap==1)
				obs_lik[i] = good;
			else if(obt==4&&hap==2)
				obs_lik[i] = good;
			else
				cout << "Invalid observation and/or haplotype " << obt << ", " << hap << endl;
		break;
		}

#ifdef DEBUG
cout << "\t" << gen_lik[i] << " " << obs_lik[i] << " " << gen_lik[i] * obs_lik[i];
#endif
		prob += gen_lik[i] * obs_lik[i];
	}
	reads_snp_list[count]->addEmission(hap,prob);
#ifdef DEBUG
cout << "\t" << prob << "\t" << log(prob) << endl;
#endif
	return prob;
}

