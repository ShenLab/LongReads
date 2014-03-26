
#include<vector>

class READ;

class SNP {

public:

  SNP(void);
  //SNP(long snp_pos, char snp_ref, char snp_alt, int type, READ *read, vector<string>, bool known);
  SNP(string chr, long snp_pos, char snp_ref, char snp_alt, int known, double qual);

  ~SNP(void);

  char *GetChr();
  long GetPos();
  void append(int type, READ *read);
  READ *GetRead(int pos);
  char GetRef();
  char GetAlt();
  int GetReadCount();
  int GetRefCount();
  int GetAltCount();
  int GetErrCount();
  double* GetGenLik();
  int GetKnown();
  int GetSom();
  double GetQualScore();
  int GetOverlapCount();
  int GetKnownOverlapCount();
  int GetEmissionState();
  void add_posteriors(double posterior[3]);
  void add_somatic_posteriors(double posterior[3]);
  void assign_genotype(int gt, double genp);
  void addEmission(int haplotype, double probability);
  void IncrOverlapCount();
  void IncrKnownOverlapCount();
  double* GetPosteriors();
  double* GetSomaticPosteriors();

private:
  char chr[10];
  int known;
  char ref;
  char alt;
  int refcount;
  int altcount;
  int errcount;
  int count;
  int overlap_count;
  int known_overlap_count;
  long position;
  double qual;
  double likelihood_ratio;
  double genprob;
  READ **reads;
  double *gl;
  double *posterior;
  double *somatic_posterior;
  double emission[3];

void add_read(int type, READ *read);

};

//===============================================================================
//===============================================================================
