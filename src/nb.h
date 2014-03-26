void GibbsSampling(int snp_start, int snp_end, int nbSymbols);
void NaiveBayes(int snp_start, int snp_end, int nbSymbols);
void GetCommonSnpList(SNP** reads_snp_list, SNP** known_snp_list, int *common_snp_count, int *known_snp_count, int *common_index, int *known_index, int t);
void GetKnownSnpList(SNP** known_snp_list, int *known_snp_count, int *known_index, int t);
int GetMismatchSum(READ *read);
void GetAllSnps(SNP** known_snp_list, int *known_snp_count, int *known_index, READ *rd);
double compute_new_emission(SNP**reads_snp_list, int count, int t, int *index, int hap);
double compute_new_emission(SNP**reads_snp_list, int count, char all1, char all2, int *index, int hap);
double compute_full_emission(SNP**reads_snp_list, int count, int *index, int hap);
double compute_full_emission(SNP**snp_list, int count, int t, int *index, int hap);
void haplotypeProbability(vector<SNP*>::iterator snp_it, double happ[3]);
int somaticHaplotypeProbability(vector<SNP*>::iterator snp_it, double happ[3], int *known_hap_count);
void genotypeProbability(vector<SNP*>::iterator snp_it, double genp[2]);
void UpdateGenotypes(long start, long end);
int UpdateBaumWelchGenotypePosteriors(long start, long end);
void FindSomaticMutations(long start, long end);
char GetPosAllele(int t, long pos);
int GetPosQual(int t, long pos);

//===============================================================================
//===============================================================================
