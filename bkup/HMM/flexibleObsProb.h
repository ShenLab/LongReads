/* ************************************************************************ *
 * ************************************************************************ *

   File: vectorObsProb.h
   The class CFlexibleObsProb defines operations 
   for the observation distribution B
   when observations are vectors with independent components

  * ************************************************************************ *

  * ************************************************************************ *

   Authors: Daniel DeMenthon & Marc Vuilleumier
   Date:  2-27-99 

 * ************************************************************************ *

   Modification Log:
	

 * ************************************************************************ *
   Log for new ideas:
 * ************************************************************************ *
               Language and Media Processing
               Center for Automation Research
               University of Maryland
               College Park, MD  20742
 * ************************************************************************ *
 * ************************************************************************ */
 
//===============================================================================

	
//===============================================================================

class CFlexibleObsProb: public CObsProb{
	
  public:
        CFlexibleObsProb(int *listNbSymbols, int nbStates, int nbComponents, int isGaussian);         	
	CFlexibleObsProb(void);
	~CFlexibleObsProb(void);
	void Start(void);
	void StartIter(void);
	void BWSum(double *gamma, CObs *obs);
	void BWSum(double *gamma, int symbol);
	void SKMSum(int state, CObs *obs);
	double EndIter();
	double EndLogIter();
	void UnLog();
	void End();
	void InitStateProb(double);
        void SetVal(int state, int symbol, double value);
	CObs* PickObservation(int state);
        CObs** MapStateToObs(void);
	void Print(ostream &outFile);
	double at(int state, CObs *obs);
	double at(int state, int symbol);
	void addEmission(int state, int symbol, double emission);
	int GetM(void){cout<<"Wrong class used"<<endl; return 0;};
#if 1
	CObs* ReadObsFrom(ifstream &inFile);
	void ReadFileHeader(ifstream &inFile);
	void PrintFileHeader(ostream &outFile);
#endif
   
  private:
	void InitParameters(void);
	
  private:
  	int mDimension;
//	CGaussianObsProb **mComponentProb;// array of ptrs to gaussian probs, size = mDimension
	CObsProb **mComponentProb;// array of ptrs to discrete or gaussian probs, size = mDimension
};

//===============================================================================
//===============================================================================

//===============================================================================
