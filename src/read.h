
class SNP;

class READ {

public:
  READ(void);
  READ(const READ&);
  READ(long start, int length);
  READ (int i);
  READ (double i);

  ~READ(void);

  void addsnp(SNP *snp, char allele, int qual, bool inp, bool delp, int known);
  void assignHaplotype(int haplotype, double prob, int flag);
  long GetPos(void);
  int GetLen(void);
  int GetHap(void);
  int GetDiscordance(void);
  double GetHapProb(void);
  void AddKnownCount(int);
  SNP *GetSnp(int pos);
  char GetAllele(int pos);
  char GetKnownAllele(int pos);
  bool GetProximalInsert(int pos);
  bool GetProximalDelete(int pos);
  int GetQualScore(int pos);
  int GetKnownQualScore(int pos);
  int GetSnpCount(void);
  int GetKnownCount(void);
  int GetKnownOverlapCount(void);
  SNP **GetSnpList();
  SNP **GetKnownList();

  READ* operator=(int i);
  READ* operator+=(double i);
  READ* operator-=(double i);
  READ* operator*=(double i);
  READ* operator/=(double i);
  READ* operator+=(READ* read);
  READ* operator-=(READ* read);
  READ* operator*=(READ* read);
  READ* operator/=(READ* read);
  operator char();
  operator int();
  operator double();
  double operator*(READ* read);
  double operator-(double i);
  double operator+(double i);

private:
  int length;
  int snp_count;
  int known_count;
  int known_overlap_count;
  int hap;
  int discordance;
  long start;
  double hprob;
  SNP **snps;
  SNP **known_snps;
  bool *inr;
  bool *delr;
  char *alleles;
  int *qualstring;
  char *known_alleles;
  int *known_qualstring;
};

