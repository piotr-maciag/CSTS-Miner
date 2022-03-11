#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string>
#include <functional>
#include <ctime>
#include <chrono>

using namespace std;
using namespace std::chrono;

struct STPoint
{
	int eventID;
	string eventType;
	double spatialX;
	double spatialY;
	double temporal;
};


struct Sequence
{
	//int seqID;
	vector<string> sequence;
	vector<vector<STPoint *>> tailEventSet;
	double PI = 1000.0;

	Sequence* firstParent = NULL;
	Sequence* secondParent = NULL;

	vector<Sequence*> children;

	vector<Sequence*> Cmax;
	vector<Sequence*> RC;

	//bool isAppClosed = true;
	bool isClosed = true;
};


extern int sizeP;
extern STPoint* dataP;
extern fstream uchwyt;

extern vector<vector <STPoint>> dataset;
extern vector<vector <STPoint>> sortedDataset;

extern vector<vector<Sequence *>> SequencesSetSPTree;

extern double theta;
extern double approxMargin;

extern double phase1Time;
extern double phase2Time;

extern double R;
extern double T;

void LoadDataset(string);
int CountInstances(string);
void InsertInstance(STPoint);
void TransformData();
void PrintDataset();
void PrintSortedDataset();

void SortDataset();

vector<STPoint> ForwardSweep(vector<STPoint>);
void PrintSequences();
int CountPatternsCSTS_Miner();
int CountPatterns();
int CountPatternsClosed();
vector<Sequence *> GenAndVerify(vector<Sequence *> SeqVec);
vector<Sequence *> GenAndVerifyCSTPM(vector<Sequence *> SeqVec);
void ExpandSequence(Sequence * Seq);
//void CalculateClosureSet(Sequence * seq);
void CSTS_Miner();
void STBFM() ;
void CST_SPMiner() ;
void STS_Miner();
void CSTPM();
void VerifyL1Parent(Sequence * seq, Sequence * seqPar);
void VerifySupersequence(Sequence * seq, Sequence * seqPar);

void ClearSequencesSet();
void ClearDataset();
void ClearStructures();
void WriteStatistics(fstream &fileRCsizes, fstream &fileCmaxsizes);

void InsertIntoSequSet(Sequence seq);




