#include<stdio.h>
#include <assert.h>
#include <math.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<string>
#include <iomanip>
#include <fstream>
#include<map>
#include<vector>
#include<sstream>
#include<algorithm>

using namespace std;

#include "read.h"
#include "snp.h"
#include "hmm.h"

/*
#include "ActionClassifier.h"
#include "BayesianClassifier.h"
#include "Domain.h"
*/

#define MAX_READS 10000
#define MAX_READ_LEN 4000
#define FULL
//#define DEBUG
//#define FULLDEBUG

int region = 21;
int coverage = 10;
double errate = 0.02;
int EM=1;
int flag_iter = 1;

vector<SNP*> snp_list;
vector<READ*> reads_list;
vector<string> known_snps;
vector<long int> known_soms;

struct gt_st {
	int known;
	double gt[3];
};

map<long, gt_st> gt_map;

struct read_span {
	long start;
	long end;
};

read_span *read_info = new read_span[1000];

vector<string> &split(const string &s, char delim, vector<string> &elems) {
	stringstream ss(s);
	string item;
	while(getline(ss, item, delim))
		elems.push_back(item);
	return elems;
}

// Get a seg fault due to mem corruption error here.
// Fixed by setting \0 at end of line copied from file
vector<string> split(const string &s, char delim,int flag) {
	vector<string> elems;
	return split(s, delim, elems);
}

int rev_s_len(const char *seq)
{
	int i = strlen(seq)-1;
	int len = 0;
	int j = 1;

	if(seq[i] != 'S')
		return 0;
	else {
		i--;
		while(seq[i]>='0'&&seq[i]<='9') {
			len = 10*j*(seq[i]-'0') + len;
			i--;
			j *= 10;
		}
		if(len==0) {
			cout << "Invalid CIGAR. Ends with " << seq[strlen(seq)-1] << endl;
			exit(1);
		}
		return len/10;
	}
}

int s_len(const char* seq)
{
	int i = 0;
	int len = 0;

	while(seq[i]>='0'&&seq[i]<='9') {
		len = 10*len + seq[i]-'0';
		i++;
	}
	if(len==0) {
		cout << "Invalid cigar. Beginning with " << seq[0] << endl;
		exit(1);
	}
	if(seq[i]=='S')
		return len;
	else
		return 0;
}

int get_deletion_count(string cigar, int limit)
{
	//char *c_cigar = new char[cigar.length()+1];
	char c_cigar[1000];
	strcpy(c_cigar,cigar.c_str());
	c_cigar[cigar.length()] = '\0';
	int it = 0;
	int len = 0;
	int dlen = 0;
	int mlen = 0;
	int ilen = 0;
	int slen = 0;
	int dmlim = 0;

	// slen serves no great purpose. Not utilized here.
	while(it<strlen(c_cigar)) {
		if(c_cigar[it]>='0'&&c_cigar[it]<='9') {
			len = len*10 + c_cigar[it] - '0';
		} else if(c_cigar[it]=='S') {
			slen += len;
			len = 0;
#ifdef FULLDEBUG
if(limit<10000) cout << "Slen = " << slen << endl;
#endif
		} else if (c_cigar[it]=='M') {
			mlen += len;
			len = 0;
#ifdef FULLDEBUG
if(limit<10000) cout << "Mlen = " << mlen << endl;
#endif
			if(mlen+dlen-ilen==limit) {
				dmlim = 0;
#ifdef FULLDEBUG
if(limit<10000) cout << "Limit= " << limit << endl << "Mlen+Dlen = " << mlen+dlen << endl;
#endif
				break;
			} else if(mlen+dlen-ilen > limit) {
				dmlim = 0;
#ifdef FULLDEBUG
if(limit<10000) cout << "Limit= " << limit << endl << "Mlen+Dlen = " << mlen+dlen << endl;
#endif
				break;
			}
		} else if(c_cigar[it]=='D') {
			dlen += len;
			len = 0;
#ifdef FULLDEBUG
if(limit<10000) cout << "Dlen = " << dlen << endl;
#endif
			if(mlen+dlen-ilen==limit) {
				dmlim = 1;
#ifdef FULLDEBUG
if(limit<10000) cout << "Limit " << limit << endl << "Dlen+Mlen = " << dlen+mlen << endl;
#endif
				break;
			} else if(mlen+dlen-ilen > limit) {
				dmlim = 1;
#ifdef FULLDEBUG
if(limit<10000) cout << "Limit " << limit << endl << "Dlen+Mlen = " << dlen+mlen << endl;
#endif
				break;
			}
		} else if(c_cigar[it]=='I') {
			ilen += len;
			len = 0;
#ifdef FULLDEBUG
if(limit<10000) cout << "Ilen = " << ilen << endl;
#endif
		} else {
			cout << "Invalid character in cigar string " << c_cigar[it] << endl;
		}
		it++;
	}

//	delete [] c_cigar;
	if(dmlim)
		return -1;
	else
		return dlen - ilen;
}

int get_read_span(const char *file_name)
{
	int iter = 0;
	int span_iter = 0;
	string line;
	char str[20000];
	FILE *file = fopen(file_name, "r");

	if(file==NULL) {
		printf("Cannot open file %s\n",file_name);
		exit(1);
	}

	while(true) {
		fgets(str,sizeof(str),file);
		str[strlen(str)-1] = '\0';
		if(feof(file)) {
			printf("Empty file %s\n",file_name);
			exit(1);
		}
		if(str[0]!='@')
			break;
	}

	iter++;
	line = str;
	vector<string> read_line = split(line, '\t',0);
	read_info[span_iter].start = atol(read_line[3].c_str());

	while(true) {
		fgets(str,sizeof(str),file);
		str[strlen(str)-1] = '\0';
		if(feof(file)) {
			read_info[span_iter].end = atol(read_line[3].c_str());
			return span_iter+1;
		}
		iter++;
		line = str;
		read_line = split(line, '\t',0);
		if(!(iter%MAX_READS)) {
			read_info[span_iter].end = atol(read_line[3].c_str());
			span_iter++;
			read_info[span_iter].start = atol(read_line[3].c_str());
		}
	}
	fclose(file);
}

void read_known_snp_file(const char *file)
{
	string line;
	FILE *snp_file = fopen(file, "r");
	if(snp_file==NULL) {
		printf("Cannot open snp file %s\n",file);
		exit(1);
	}

	while(true) {
		char str[50];
		fgets(str,sizeof(str),snp_file);
		str[strlen(str)-1] = '\0';
		if(feof(snp_file))
			break;
		line = str;
		vector<string> snp_line = split(line, ' ',0);

		// Consider hets only. Others are not useful for haplotype calling
		if(snp_line[4]==snp_line[5])
			continue;
		string snpit = snp_line[1] + snp_line[4] + snp_line[5] + snp_line[8];
		known_snps.insert(known_snps.end(),snpit);
	}
}

void read_known_som_file(const char *file)
{
	long int pos;
	FILE *som_file = fopen(file, "r");
	if(som_file==NULL) {
		printf("Cannot open som file %s\n",file);
		exit(1);
	}

	while(true) {
		char str[50];
		fgets(str,sizeof(str),som_file);
		str[strlen(str)-1] = '\0';
		if(feof(som_file))
			break;
		pos = atol(str);

		known_soms.insert(known_soms.end(), pos);
	}
	std::sort(known_soms.begin(), known_soms.end());
}

void read_snp_file(FILE *snp_file, long snp_end)
{
	string line;

	vector<long int>::iterator som_it;
	vector<string>::iterator vec_start = known_snps.begin();
	while(true) {
		char str[1000];
		fpos_t prev_pos;
		fgetpos(snp_file,&prev_pos);
		fgets(str,sizeof(str),snp_file);
		str[strlen(str)-1] = '\0';
		if(feof(snp_file))
			break;
		if(str[0]=='#')
			continue;
		line = str;
		vector<string> vcf_line = split(line, '\t',0);

		// Exclude non-biallelic snps and non-hets
		if(vcf_line[4].find(',')==string::npos) {
		//if(vcf_line[7].find("DP=2;")==string::npos && vcf_line[7].find("DP=3;")==string::npos) { // && vcf_line[7].find("DP=4;")==string::npos && vcf_line[7].find("DP=5;")==string::npos) {
			long pos = atol(vcf_line[1].c_str());
			if(pos>snp_end) {
				fsetpos(snp_file,&prev_pos);
				break;
			}

			string info = vcf_line[9];
			vector<string> format = split(info, ':',0);
			vector<string> gl3 = split(format[1], ',',0);

			// known = 0(error), 1(known), 2(novel), 3(null ref) - discarding value 3
			int known = 0, a1 = 0, a2 = 0, som = 0;
			for(vector<string>::iterator known_snpit = vec_start; known_snpit != known_snps.end(); known_snpit++) {
				long int full_info = atol((*known_snpit).c_str());
				if( full_info/1000 > pos) {
					break;
				}
				if(full_info/1000 == pos) {
					known = full_info%10 == 0 ? 2 : 1;
					a1 = (full_info/10)%10;
					a2 = (full_info/100)%10;
					vec_start = known_snpit;
					break;
				}
			}
			// Changing quality score
			som_it = find(known_soms.begin(), known_soms.end(), pos);
			som = (som_it==known_soms.end()) ? 0 : 1;
			SNP *snp = new SNP(vcf_line[0], pos, vcf_line[3].c_str()[0], vcf_line[4].c_str()[0], gl3, known, som, atof(vcf_line[5].c_str()));
			snp_list.insert(snp_list.end(), snp);
#ifdef FULLDEBUG
cout << "Insert.. " << pos << endl;
#endif
		//}
		}
	}
	cout << "SNP MAP size = " << snp_list.size() << endl;
}

int read_sam_files(FILE *fq_file)
{
	int i;
	int iter = 0;
	string line;

	vector<SNP*>::iterator snp_begin = snp_list.begin();
	while(iter<MAX_READS) {
		char str[20000];
		fgets(str,sizeof(str),fq_file);
		str[strlen(str)-1] = '\0';
		if(feof(fq_file))
			break;
		char star = str[0];

		if(star != '@') { // for each new read in the file
			line = str;
			vector<string> read_line = split(line, '\t', 0);
			string hapl_s = read_line[0];
			vector<string> hapl_v = split(hapl_s,':',0);
			int hapl_i = atoi(hapl_v[2].c_str());
			long start = atol(read_line[3].c_str());
#ifdef DEBUG
cout << "Read: " << start << " has haplotype: " << hapl_i << endl;
#endif
			string cigar = read_line[5];
			int slen = s_len(cigar.c_str());
			int rev_slen = rev_s_len(cigar.c_str());
			int Ndels = get_deletion_count(cigar,10000);
if(Ndels==-1) {
	cout <<"Unexpected -1 for get_deletion_count" << endl;
}
			int length = (int)(read_line[9].length()) - slen - rev_slen + Ndels;
#ifdef FULLDEBUG
cout << slen << "\t" << rev_slen << "\t" << Ndels << "\t" << endl;
#endif
#ifdef FULLDEBUG
cout << "Start " << start << ", Length = " << length << endl;
#endif
			READ *read = new READ(start,length,hapl_i); // create new read object

#ifdef FULLDEBUG
cout << start << "\t" << (*snp_begin)->GetPos() << endl;
#endif
			for(vector<SNP*>::iterator snp_it = snp_begin; snp_it != snp_list.end(); snp_it++) {
				long snp_pos = (*snp_it)->GetPos();
				if(snp_pos >= start) {
					if(snp_pos <= start+length-1) {
						int snp_exists = 0;
						int type = 2;
						char snp_ref = (*snp_it)->GetRef();
						char alt_ref = (*snp_it)->GetAlt();
#ifdef FULLDEBUG
cout << snp_pos << "\t" << start << endl;
#endif
						int ndels = get_deletion_count(cigar,snp_pos-start+1);
#ifdef FULLDEBUG
cout << "Ndels = " << ndels << endl;
#endif
						if(ndels==-1)
							continue;

						char allele = read_line[9][snp_pos - start + slen - ndels];
						int qual = read_line[10][snp_pos - start + slen - ndels];

						if(allele==snp_ref)
							type = 0;
						else if(allele==alt_ref)
							type = 1;

						read->addsnp(*snp_it, allele, qual);
						(*snp_it)->append(type,read);
					} else {
						break;
					}
				} else {
					snp_begin = snp_it;
				}
			}
			reads_list.insert(reads_list.end(), read); // add it to the vector
			iter++;
		}
	}

	cout << "number of reads = " << reads_list.size() << endl;
	cout << "number of snps = " << snp_list.size() << endl;
	return iter;
	//return reads_list.size();
}

double runNB(ofstream &stateOutput, ofstream &distanceOutput, double initProb, long snp_start, long snp_end, int nbSymbols)
{
	int i;
	int nbDimensions = 1, nbStates = 2;
	int isGaussian = 0;
	int *listNbSymbols = new int[(nbDimensions)+1];;
	double **transitionMatrix, **emissionMatrix;
	double transProb = 0.5, emissionProb;

	for (i=1;i<=nbDimensions;i++){
		listNbSymbols[i] = nbSymbols;
	}
#ifdef DEBUG
cout << "Nbsymbols = " << nbSymbols << endl;
#endif
	//NaiveBayes(snp_start,snp_end, nbSymbols);
	GibbsSampling(snp_start,snp_end, nbSymbols);

	for(vector<SNP*>::iterator snp_it = snp_list.begin(); snp_it != snp_list.end(); snp_it++) { // check for snps in vector
		long pos = (*snp_it)->GetPos();
		if(pos<=snp_start)
			continue;
		else if(pos>snp_end)
			break;

		double gt[3], gt2[3], normal, st[3];
		gt[0] = (*snp_it)->GetPosteriors()[0];
		gt[1] = (*snp_it)->GetPosteriors()[1];
		gt[2] = (*snp_it)->GetPosteriors()[2];
		gt2[0] = (*snp_it)->GetGenLik()[0];
		gt2[1] = (*snp_it)->GetGenLik()[1];
		gt2[2] = (*snp_it)->GetGenLik()[2];
		st[0] = (*snp_it)->GetSomaticPosteriors()[0];
		st[1] = (*snp_it)->GetSomaticPosteriors()[1];
		st[2] = (*snp_it)->GetSomaticPosteriors()[2];

		int gp = gt[0]>gt[1] ? (gt[0]>gt[2] ? 0:2) : (gt[1]>gt[2] ? 1:2);

		// Assigning priors to likelihood for calculating posteriors.
		double known_prior[3], novel_prior[3], prior[3], sam_post[3], hmm_post[3], som_post[3];

		known_prior[1] = 0.1;
		known_prior[2] = known_prior[1]*known_prior[1];
		known_prior[0] = 1 - known_prior[1] - known_prior[2];
		novel_prior[1] = 0.001;
		novel_prior[2] = novel_prior[1]*novel_prior[1];
		novel_prior[0] = 1 - novel_prior[1] - novel_prior[2];

		if((*snp_it)->GetKnown()==1) {
			prior[0] = known_prior[0];
			prior[1] = known_prior[1];
			prior[2] = known_prior[2];
		} else {
			prior[0] = novel_prior[0];
			prior[1] = novel_prior[1];
			prior[2] = novel_prior[2];
		}
		prior[0] = 1;
		prior[1] = 1;
		prior[2] = 1;
		sam_post[0] = gt2[0]*prior[0];
		sam_post[1] = gt2[1]*prior[1];
		sam_post[2] = gt2[2]*prior[2];
		hmm_post[0] = gt[0]*prior[0];
		hmm_post[1] = gt[1]*prior[1];
		hmm_post[2] = gt[2]*prior[2];
		som_post[0] = st[0]*prior[0];
		som_post[1] = st[1]*prior[1];
		som_post[2] = st[2]*prior[2];

		normal = sam_post[0] + sam_post[1] + sam_post[2];
		sam_post[0] /= normal;
		sam_post[1] /= normal;
		sam_post[2] /= normal;
		normal = hmm_post[0] + hmm_post[1] + hmm_post[2];
		hmm_post[0] /= normal;
		hmm_post[1] /= normal;
		hmm_post[2] /= normal;
		normal = som_post[0] + som_post[1] + som_post[2];
		som_post[0] /= normal;
		som_post[1] /= normal;
		som_post[2] /= normal;

		stateOutput << (*snp_it)->GetChr() << "\t" << (*snp_it)->GetPos() << "\t" << (*snp_it)->GetKnown() << "\t" << gp << "\t" << hmm_post[0] << "\t" << hmm_post[1] << "\t" << hmm_post[2] << "\t" << sam_post[0] << "\t" << sam_post[1] << "\t" << sam_post[2] << "\t" << sam_post[1]/sam_post[0] << "\t" << hmm_post[1]/hmm_post[0] << "\t" << (*snp_it)->GetSom() << "\t" << som_post[1] << endl;
	}
	for(vector<READ*>::iterator read_it = reads_list.begin()+reads_list.size()-(long)nbSymbols+1; read_it != reads_list.end(); read_it++) {
		distanceOutput << (*read_it)->GetHap() << "\t" << (*read_it)->GetPos() << "\t" << (*read_it)->GetSnpCount() << "\t" << (*read_it)->GetKnownCount() << endl;
	}

	int highHap = (*reads_list.begin())->GetHap();
	double newProb = (*reads_list.begin())->GetHapProb();
	if(highHap==1)
		newProb = newProb;
	else if(highHap==2)
		newProb = 1.0 - newProb;
	else
		cout << "Incorrect Haplotype for last read: " << newProb << endl;

	return newProb;
}

void delete_objects(long snp_start, long snp_end)
{
	for(vector<SNP*>::iterator snp_it = snp_list.begin(); snp_it != snp_list.end();) {
#ifdef FULLDEBUG
cout << (*snp_it)->GetPos() << endl;
#endif
		if((*snp_it)->GetPos() < snp_start) {
			cout << "Unexpected snp below lower limit. Had to be deleted earlier: " << (*snp_it)->GetPos() << endl;
		} else if((*snp_it)->GetPos() >= snp_end) {
			break;
		} else {
			cout << "Deleting snp object: " << (*snp_it)->GetPos() << endl;
			delete (*snp_it);
			snp_list.erase(snp_list.begin());
		}
	}

	for(vector<READ*>::iterator read_it = reads_list.begin(); read_it != reads_list.end();) {
#ifdef FULLDEBUG
cout << (*read_it)->GetPos() << endl;
#endif
		if((*read_it)->GetPos() < snp_start) {
			cout << "Deleting read object: " << (*read_it)->GetPos() << endl;
			delete (*read_it);
			reads_list.erase(reads_list.begin());
		} else if((*read_it)->GetPos()+(*read_it)->GetLen() >= snp_end) {
			break;
		} else {
			cout << "Deleting read object: " << (*read_it)->GetPos() << endl;
			delete (*read_it);
			reads_list.erase(reads_list.begin());
		}
	}
}

int main(int argc, char **argv)
{
	int i;
	const char *ext = ".sam";
cout << argc << endl;

#ifdef FULL
	const char *snpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/vcfs/0.02_0.02_10_4000_21.vcf";
	const char *knownsnps="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/snp_21.list";
	const char *knownsoms="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/som_0.02_0.02_10_4000_21.list";
	const char *file_base = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/reads/0.02_0.02_10_4000_21.sorted";
	const char *base_name="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output/0.02_0.02_10_4000_1_21";
	char knownsnpfile[200];
	const char *knownsomfile;
	const char *base;

	if(argc>1) {
		file_base=argv[1];
		snpfile=argv[2];
		region=atoi(argv[3]);
		sprintf(knownsnpfile,"/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/snp_%d.list",region);
		knownsomfile=argv[4];
		base=argv[5];
		coverage=atol(argv[6]);
		errate = atof(argv[7]);
	} else {
		strcpy(knownsnpfile,knownsnps);
		knownsomfile=knownsoms;
		base=base_name;
	}
#else
	const char *snpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/short_read_1_21.1.1.vcf";
	const char *file_base = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/short_read_1.sorted";
	const char *knownsnpfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/snp_21.list";
	const char *knownsomfile="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/snps/som_0.02_0.02_10_4000_21.list";
	const char *base="/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/HMM/output/output_short_1";
#endif

	char file_name[200];
	sprintf(file_name, "%s%s", file_base, ext);

	read_known_snp_file(knownsnpfile);
	int span = get_read_span(file_name);
	read_known_som_file(knownsomfile);

	char stateOutputName[256];
	char distanceOutputName[256];

	sprintf(stateOutputName,"%s%s", base, ".sta");
	sprintf(distanceOutputName,"%s%s", base, ".obs");

	int iter=0;
	double prob1 = 0.99;
	FILE *snp_file = fopen(snpfile, "r");
	if(snp_file==NULL) {
		printf("Cannot open snp file %s\n",snpfile);
		exit(1);
	}
	FILE *fq_file = fopen(file_name, "r");
	if(fq_file==NULL) {
		printf("Cannot open fq file %s\n",file_name);
		exit(1);
	}
	ofstream stateOutput(stateOutputName);
	ofstream distanceOutput(distanceOutputName);

	//while(iter<span) {
	while(iter<1) {
		read_snp_file(snp_file, read_info[iter].end + MAX_READ_LEN); // Reads specified range or from pos 1
		int read_ct = read_sam_files(fq_file); // Reads 10000 at a time
		if(iter==0) {
			distanceOutput << "1\t" << reads_list[0]->GetPos() << "\t" << reads_list[0]->GetSnpCount() << "\t" << reads_list[0]->GetKnownCount() << endl;
		} else {
			read_ct++;
		}
		// Runs haplotype calling on [start,end] start position range of reads
		// Runs genotype calling on (start,end]
		prob1 = runNB(stateOutput, distanceOutput, prob1, read_info[iter].start, read_info[iter].end, read_ct);
		// Delete reads [start,end.end<end]
		// Delete snps [start,end)
		delete_objects(read_info[iter].start, read_info[iter].end);
		iter++;
	}

		distanceOutput.close();
		stateOutput.close();
		fclose(fq_file);
		fclose(snp_file);
	return 0;
}

