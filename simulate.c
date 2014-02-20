// Usage information
// > ./a.out -snp_rate=<double>, -err_rate=<double>, -coverage=<int>, -max_read_len=<int>, -min_read_len=<int>, -ref_file=<string> -snp_file=<string> -indel_file=<string> -ncells=<int> -out_base=<string> -region=<int>
//
// All parameters are optional. Check is performed for duplicate parametes, but not for inconsistency of values.

#include<stdio.h>
#include <stdlib.h>
#include<string.h>

#define FULL
//#define DEBUG

#define NUM_SOM 10000
#define SOM_MIN 10
#define SOM_RANGE 40

int reads_hap[111] = {1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1,2,1,1,1,1,2,2,2,1,2,1,2,1,1,2,2,1,2,1,2,2,1,2,1,1,2,1,2,1,2,1,2,2,1};
//long int nbp[23] = {0, 24724971, 24295114, 19950182, 19127306, 18085786, 17089999, 15882142, 14627482, 14027325, 13537473, 13445238, 13234953, 11414298, 10636858, 10033891, 8882725, 7877474, 7611715, 6381165, 6243596, 4694432, 4969143};
//long int nbp[23] = {5000000, 247249719, 242951149, 199501827, 191273063, 180857866, 170899992, 158821424, 146274826, 140273252, 135374737, 134452384, 132349534, 114142980, 106368585, 100338915, 88827254, 78774742, 76117153, 63811651, 62435964, 46944323, 49691432};
long int nbp[23] = {5000000,249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566};

struct snp_struct {
	int position;
	char all[2];
};

struct indel_struct {
	int position;
	char reflen;
	char all[2][100];
};

struct read_struct {
	char *str;
	long start;
	int len;
	double del_errors;
	double snp_errors;
};

struct som_struct {
	int position;
	int prevelance;
	char alt;
};

// By default, it will run only for chr21 if FULL is defined, or a small portion of it if FULL is not defined.
#ifdef FULL
struct snp_struct snps[400000];
struct indel_struct indels[40000];
struct som_struct sites[1000];
//long int N_BP = 46944323; // chr21
long int N_BP = 48129895; // chr21
int reads_max = 150000;
int region = 21;
#else
struct snp_struct snps[400000];
struct indel_struct indels[40000];
struct som_struct site[1000];
long int N_BP = 2000000;
int reads_max = 10000;
int region = 0;
#endif

int CLT = 30;

double def_snp_rate = 0.02;
double def_err_rate = 0.02;
double snp_rate;
double err_rate;

char A[3] = "TCG";
char T[3] = "CGA";
char C[3] = "GAT";
char G[3] = "ATC";

int coverage;
int prevalence = 25;
int def_coverage = 30;
int read_len = 4000;
int max_ref_len = 250000000;
char *snp_file, *indel_file, *som_file;
FILE *fq_file, *log_file;
const char *file_base = "short_read_1";
//const char *ref_file = "/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/reference/bcm_hg18.fasta";
const char *ref_file = "reference/hg19.fasta";

int read_cmp(const void *r1, const void *r2)
{
	return (((struct read_struct*)r1)->start < ((struct read_struct*)r2)->start ? -1 : (((struct read_struct*)r1)->start == ((struct read_struct*)r2)->start ? 0 : 1) );
}

void simulate(int chrnum);

int main(int argc, char **argv)
{
	snp_rate = def_snp_rate;
	err_rate = def_err_rate;
	coverage = def_coverage;
	snp_file = (char*)malloc(sizeof(char)*200);
	indel_file = (char*)malloc(sizeof(char)*200);
	som_file = (char*)malloc(sizeof(char)*200);

#ifdef DEBUG
printf("RAND_MAX = %ld\n", RAND_MAX);
#endif
	// override parameters
	if(argc>1) {
		int param = 1;
		int snp_ck = 0, err_ck = 0, len_ck = 0, cov_ck = 0, ref_ck = 0, snf_ck = 0, ind_ck = 0, som_ck = 0, cel_ck = 0, out_ck = 0, reg_ck = 0, prev_ck = 0;
		for(param=1;param<argc;param++) {
			char *token[2];
			token[0] = strtok(argv[param], "=");
			if(token[0]==NULL) {
				printf("Ignoring incorrect parameter #%d\n", param);
			} else if((token[1] = strtok(NULL, "=")) == NULL) {
				printf("Ignoring incorrect parameter #%d\n", param);
			} else {
				if(strcmp(token[0],"-snp_rate")==0) {
					if(snp_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						snp_ck++;
						snp_rate = atof((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-err_rate")==0) {
					if(err_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						err_ck++;
						err_rate = atof((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-read_len")==0) {
					if(len_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						len_ck++;
						read_len = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-coverage")==0) {
					if(cov_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						cov_ck++;
						coverage = atoi((const char *)token[1]);
						if(coverage<2) {
							printf ("Insufficient coverage %dX provided. Ignoring..\n",coverage);
							coverage = def_coverage;
						}
					}
				} else if(strcmp(token[0],"-ref_file")==0) {
					if(ref_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						ref_ck++;
						ref_file = token[1];
					}
				} else if(strcmp(token[0],"-snp_file")==0) {
					if(snf_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						snf_ck++;
						snp_file = token[1];
					}
				} else if(strcmp(token[0],"-indel_file")==0) {
					if(ind_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						ind_ck++;
						indel_file = token[1];
					}
				} else if(strcmp(token[0],"-som_file")==0) {
					if(som_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						som_ck++;
						som_file = token[1];
					}
				} else if(strcmp(token[0],"-out_base")==0) {
					if(out_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						out_ck++;
						file_base = token[1];
					}
				} else if(strcmp(token[0],"-region")==0) {
					if(reg_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						reg_ck++;
						region = atoi((const char *)token[1]);
					}
				} else if(strcmp(token[0],"-freq")==0) {
					if(prev_ck>0) {
						printf("Duplicate parameter %s=%s, ignoring..\n", token[0],token[1]);
					} else {
						prev_ck++;
						prevalence = atoi((const char *)token[1]);
					}
				} else {
					printf("Invalid parameter %s=%s. Ignoring..\n", token[0],token[1]);
				}
			}
		}
		if(snp_rate>err_rate) {
			printf("SNP rate (default 0.02) should be lower than error rate (default 0.02). Using default values now...\n");
//				snp_rate = def_snp_rate;
//				err_rate = def_err_rate;
		}
	}

	int chrnum;
	char file_name[200];
	char *ext = ".fq";
	char *log = ".log";

	sprintf(file_name, "%s%s", file_base, ext);
	fq_file = fopen(file_name, "w");

	if (fq_file == NULL)
		printf("Error opening output file [%s]\n", (const char *) file_name);

	simulate(region);
	fclose(fq_file);
}

void simulate(int chrnum)
{
	int c=0;
	long int i=0;
	N_BP = nbp[chrnum];
	int som_ct = (NUM_SOM * N_BP) / 3000000000;
	sprintf(snp_file, "snps/snp_%d.list",chrnum);
	sprintf(indel_file, "snps/indel_%d.list",chrnum);
	long int n_reads = (nbp[chrnum] * coverage)/read_len;
//n_reads = 50000;
	struct read_struct *reads = (struct read_struct*)malloc(sizeof(struct read_struct)*n_reads);

	for(c=0;c<n_reads;c++) {
		reads[c].str = NULL;
	}

	int ctr=0;
	char line[100];
	char *ref_s = (char*) malloc(sizeof(char)*max_ref_len);
	
	// Read ref fasta here
	FILE *ref = fopen(ref_file,"rt");
	if(ref == NULL) { printf("Can't open ref file \"%s\"\n", ref_file); exit(1); }

	int read_flag = 0;
	if(chrnum==0)
		chrnum=21;
	while(1) {
		fgets(line,sizeof(line),ref);
		if(feof(ref))
			break;
		int c_ptr = 0;
		int line_len = strlen(line) - 1;
		if(line[0]=='>') {
			if(read_flag==1) {
				break;
			} else {
				if((line_len==5&&line[4]-48==chrnum)||(line_len==6&&(line[4]-48)*10+line[5]-48==chrnum))
					read_flag = 1;
					continue;
			}
		}
		if(read_flag==0)
			continue;
		for(c_ptr = 0;c_ptr<line_len;c_ptr++) {
			ref_s[ctr+c_ptr] = line[c_ptr];
		}
		ctr += line_len;
	}
	ref_s[ctr] = '\0';
	fclose(ref);

	// Read snp file here
	FILE *snp_f = fopen(snp_file,"rt");
	if(snp_f == NULL) { printf("Can't open snp file \"%s\"\n", snp_file); exit(1); }

	int sctr = 0;
	while(1) {
		int pos=0,hap=0;
		char ref,alt;

		fgets(line,sizeof(line),snp_f);
		if(feof(snp_f))
			break;
		strtok(line,"\t");
		pos = atoi(strtok(NULL,"\t"));
		strtok(NULL,"\t");
		//strtok(NULL," ");
		//strtok(NULL," ");
		//strtok(NULL," ");
		ref = strtok(NULL,"\t")[0];
		alt = strtok(NULL,"\t")[0];

		// Consider hets only. Everything else is useless for haplotype calling
		if(ref==alt)
			continue;
		snps[sctr].position = pos;
		snps[sctr].all[0] = ref;
		snps[sctr].all[1] = alt;
		sctr++;
	}
	fclose(snp_f);

	// Read indel file here
	FILE *indel_f = fopen(indel_file,"rt");
	if(indel_f == NULL) { printf("Can't open indel file \"%s\"\n", indel_f); exit(1); }

	int ictr = 0;
	while(1) {
		int pos=0,hap=0,sit=0;
		char sptr[100];

		fgets(line,sizeof(line),indel_f);
		line[strlen(line)-1] = '\0';
		if(feof(indel_f))
			break;
		strtok(line,"\t");
		indels[ictr].position = atoi(strtok(NULL,"\t"));
		strtok(NULL,"\t");
		memset(sptr,'\0',100);
		strcpy(sptr,strtok(NULL,"\t"));
		indels[ictr].reflen = strlen(sptr);
		//strtok(NULL," ");
		//strtok(NULL," ");
		memset(indels[ictr].all[0],'\0',100);
		memset(indels[ictr].all[1],'\0',100);
		//strcpy(indels[ictr].all[0],strtok(NULL," "));
		strcpy(indels[ictr].all[0],sptr);
		strcpy(indels[ictr].all[1],strtok(NULL,"\t"));
		ictr++;
	}
	fclose(indel_f);

	srand ( (unsigned)time ( NULL ) );
	for(i=0;i<n_reads;i++) {
#ifdef FULL
		reads[i].start = (int)((rand()*(N_BP-read_len))/RAND_MAX) + 1;
#else
		int read_position_start = 9880000;
		reads[i].start = (int)((rand()*N_BP)/RAND_MAX) + read_position_start;
#endif
	}
	for(i=0;i<n_reads;i++) {
		reads[i].len = read_len;
		reads[i].str = (char*)malloc(sizeof(char)*(reads[i].len+1));
	}
	for(i=0;i<n_reads;i++) {
		int clt;
		double sum=0, ssum=0;
		for(clt=0;clt<CLT;clt++) {
			sum += (((rand()*(err_rate*1000))/RAND_MAX) + (err_rate*1000/2));
			ssum += (((rand()*(snp_rate*1000))/RAND_MAX) + (snp_rate*1000/2));
		}
		reads[i].del_errors = sum/(CLT*10);
		reads[i].snp_errors = ssum/(CLT*10);
		sum = 0; ssum = 0;
	}

	char buffer[20];
	FILE *som_f = fopen(som_file,"w");
	if(som_f == NULL) { printf("Can't open som file \"%s\"\n", som_file); exit(1); }

	for(i=0;i<som_ct;i++) {
		sites[i].position = (int)((rand()*N_BP)/RAND_MAX) + 1;
		sites[i].prevelance = prevalence;

		if(ref_s[sites[i].position-1] == 'A' || ref_s[sites[i].position-1] == 'a')
			sites[i].alt = A[(int)((long int)3*rand()/RAND_MAX)];
		else if(ref_s[sites[i].position-1] == 'C' || ref_s[sites[i].position-1] == 'c')
			sites[i].alt = C[(int)((long int)3*rand()/RAND_MAX)];
		else if(ref_s[sites[i].position-1] == 'G' || ref_s[sites[i].position-1] == 'g')
			sites[i].alt = G[(int)((long int)3*rand()/RAND_MAX)];
		else if(ref_s[sites[i].position-1] == 'T' || ref_s[sites[i].position-1] == 't')
			sites[i].alt = T[(int)((long int)3*rand()/RAND_MAX)];
		else if(ref_s[sites[i].position-1] == 'N' || ref_s[sites[i].position-1] == 'n') {
			i--;
			continue;
		} else
			printf("Encountered unexpected base %c at position %d\n", ref_s[sites[i].position], sites[i].position);

		memset(buffer,'\0',20);
		sprintf(buffer, "%d", sites[i].position);
		fputs(buffer,som_f);
		fputs("\n",som_f);

#ifdef DEBUG
printf("Created mutation %d at position %d with prevelance %d and allele %c for original %c\n", i, sites[i].position, sites[i].prevelance, sites[i].alt, ref_s[sites[i].position-1]);
#endif
	}
	fclose(som_f);

#ifdef DEBUG
printf("Before sorting\n");
for(i=0; i<n_reads; i++) {
	printf("Read %d starts at %d,%d,%f,%f\n", i, reads[i].start,reads[i].len,reads[i].del_errors,reads[i].snp_errors);
}
#endif
	qsort(reads, n_reads, sizeof(struct read_struct), read_cmp);

#ifdef DEBUG
printf("After sorting\n");
for(i=0; i<n_reads; i++) {
	printf("Read %d starts at %d,%d,%d,%f,%f\n", i, reads[i].start,reads[i].len,reads[i].start+reads[i].len,reads[i].del_errors,reads[i].snp_errors);
}
#endif

	int snp_ct_start = 0, indel_ct_start = 0;
	char *qualstr = (char *)malloc(sizeof(char) * read_len);

	for(i=0;i<n_reads;i++) {
		int j, k = 0;
		int qsum;

		int j_st = 0, k_st = 0;
		int snp_flag = 0, indel_flag = 0;
		for(j=0;j<reads[i].len;j++) {
			int homon_hash = 100;
			int snp_ct = 0, common_snp = 0, som = 0, som_snp = 0;
			int ctj_start = snp_ct_start > j_st ? snp_ct_start : j_st;
			for(snp_ct = ctj_start; snp_ct <= ctr; snp_ct++) {
				common_snp = 0;
				if(snps[snp_ct].position < reads[i].start) {
					if(snp_flag==0) {
						snp_ct_start = snp_ct+1;
						snp_flag = 1;
					}
				} else if(snps[snp_ct].position == reads[i].start+j+1) {
					j_st = snp_ct+1;
					reads[i].str[k] = snps[snp_ct].all[reads_hap[i%100]-1];
					qualstr[k]='I';
					qsum = rand()%10 + 31;
					qsum += 33;
#ifdef FULLDEBUG
printf("True\t%d\n",qsum);
#endif
					qualstr[k]=(int)qsum;
					common_snp = 1;
#ifdef DEBUG
printf("Replacing common snp %c with %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),snps[snp_ct].all[reads_hap[i%100]-1],reads[i].start+1,k+1,j+1,reads[i].start+j+1);
#endif
					break;
				} else {
					break;
				}
			}

			int indel_ct = 0, common_indel = 0;
			int ctk_start = indel_ct_start > k_st ? indel_ct_start : k_st;
			for(indel_ct = ctk_start; indel_ct <= ctr; indel_ct++) {
				common_indel = 0;
				if(indels[indel_ct].position < reads[i].start) {
					if(indel_flag==0) {
						indel_ct_start = indel_ct+1;
						indel_flag = 1;
					}
				} else if(indels[indel_ct].position == reads[i].start+j+1) {
					k_st = indel_ct+1;
					int inct = 0;
					while(reads[i].str[k] = (indels[indel_ct].all[reads_hap[i%100]-1])[inct]) {
						qualstr[k]='I';
					qsum = rand()%10 + 31;
					qsum += 33;
#ifdef FULLDEBUG
printf("True\t%d\n",qsum);
#endif
					qualstr[k]=(int)qsum;
						k++;inct++;
					}
					k--;
					common_indel = 1;
					//j += indels[indel_ct].reflen - strlen(indels[indel_ct].all[reads_hap[i%100]-1]);
					j += (indels[indel_ct].reflen - 1);
#ifdef DEBUG
printf("Indeling site %c with %s on read %d up to position %d,%d,%d\n",*(ref_s+reads[i].start+j),indels[indel_ct].all[reads_hap[i%100]-1],reads[i].start+1,k,j+1,reads[i].start+j+1);
#endif
					break;
				} else {
					break;
				}
			}
			if(common_indel==1) continue;

			for(som=0;som<som_ct;som++) {
				som_snp = 0;
				if(sites[som].position == reads[i].start+j+1) {
					if(som%2+reads_hap[i%100]%2==1) {
						long int rnd = (long int)100*rand()/RAND_MAX;
						if(rnd <= sites[som].prevelance) {
							reads[i].str[k] = sites[som].alt;
							qualstr[k]='I';
							qsum = rand()%10 + 31;
							qsum += 33;
							qualstr[k]=(int)qsum;
							som_snp = 1;
#ifdef DEBUG
printf("Adding somatic mutation at ref %c to %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),sites[som].alt,reads[i].start+1,k+1,j+1,reads[i].start+j+1);
#endif
						}
					}
					break;
				}
			}

			double check = (100*(double)rand())/RAND_MAX;
			if(check > reads[i].del_errors+reads[i].snp_errors) {
				if(common_snp==0&&som_snp==0) {
					reads[i].str[k] = *(ref_s+reads[i].start+j);
					qualstr[k]='I';
					qsum = rand()%10 + 31;
					qsum += 33;
#ifdef FULLDEBUG
printf("True\t%d\n",qsum);
#endif
					qualstr[k]=(int)qsum;
#ifdef FULLDEBUG
printf("inserting %c on read %d at position %d,%d,%d\n",reads[i].str[k],reads[i].start+1,k+1,j+1,reads[i].start+j+1);
#endif
				}
				k++;
			} else if(check <= reads[i].snp_errors) {
				if(*(ref_s+reads[i].start+j) == 'A'||*(ref_s+reads[i].start+j) == 'a') {
					reads[i].str[k] = A[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+reads[i].start+j) == 'T'||*(ref_s+reads[i].start+j) == 't') {
					reads[i].str[k] = T[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+reads[i].start+j) == 'C'||*(ref_s+reads[i].start+j) == 'c') {
					reads[i].str[k] = C[(int)((long int)3*rand()/RAND_MAX)];
				} else if(*(ref_s+reads[i].start+j) == 'G'||*(ref_s+reads[i].start+j) == 'g') {
					reads[i].str[k] = G[(int)((long int)3*rand()/RAND_MAX)];
				} else { reads[i].str[k] = *(ref_s+reads[i].start+j); }
				qualstr[k]='!';
				qsum = rand()%10 + 21;
				qsum += 33;
#ifdef FULLDEBUG
printf("False\t%d\n",qsum);
#endif
				qualstr[k]=(int)qsum;
				k++;
#ifdef DEBUG
printf("replacing %c with %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),reads[i].str[k-1],reads[i].start+1,k,j+1,reads[i].start+j+1);
#endif
			} else {
#ifdef DEBUG
printf("Deleting %c on read %d at position %d,%d,%d\n",*(ref_s+reads[i].start+j),reads[i].start+1,k+1,j+1,reads[i].start+j+1);
#endif
			}
		}
		reads[i].str[k] = '\0';
		qualstr[k] = '\0';

		fprintf(fq_file, "@Start_POS:%d:%d\n", reads[i].start+1,reads_hap[i%100]);
		fprintf(fq_file, "%s\n", reads[i].str);
		fprintf(fq_file, "+\n");
		fprintf(fq_file, "%s\n", qualstr);
	} // end n_reads
	for(i=0;i<n_reads;i++) {
		free(reads[i].str);
	}
	free(qualstr);
	free(ref_s);
}

