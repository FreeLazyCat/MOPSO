#ifndef MOPSO_H
#define MOPSO_H

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <list>
#include <cstdlib>
#include <ctime>
#include <pthread.h>

using namespace std;

struct structBest;
struct Particle;
struct Rep;
struct Atom;
struct POINT;
struct VECTOR;
struct ANGLE;
struct LAMBDA;

typedef list<Rep> myRep;            //   a wonderful list, the data struct of the repositpry, which
typedef const char CONCH;
//store the best solutions that have been found (pateto front)


// fileDisposal
    // Input Data

int getLines(char *address);
void getAtom(Atom *atom, int cntLines, int &atomNum, char *address);        // find the word "ATOM" in the file

char * inputSeq();
void inputParticle(Particle &particle, int k, char *seq);
void inputParticles(Particle *particle);
void inputParOrigin(Particle &particle, int fileNum, int cntLines);
void inputParPhi(Particle &particle, int fileNum);
    // a function that will open file and throw error if it can't open file
void openFile(const char fileDir[10000], ifstream &infile);
    // a function that will remove file and throw error if it can't remove file
void removeFile(const char * str);

    // for the object function
void printPdb(Particle * particle);
    // for debug
        // print rep in several files
void printRep(myRep & rep, int len);
        // print the particles' status in one file
void printParticleCost(Particle * particle, int iterator);
    // print answer rep
    void outputAnswer(myRep &rep);
    //  print run time
    void printTime(const int choice, const int loopTimes);
    //  create a new fold for the answer
void createNewFold();

//  multiply better input
void multiplyBetterInput();
void multiplyParticle(Particle particle, int repNum);

void runScwrl(int index);                   //  run scwrl program after print .pdb

// get the cost for the object function 1 and 2
double getCost0(int i);
double getCost1(int i);
double getCost2(int i);
void rewriteTemplate(int i, const char *templateFile, const char *tmpIFile, const char *tmpI);
double getEnergyValue(const char *fileName);
void printCurrentId(int id);

// initialize particle
void initializeParticles(Particle *particle);
void initializeParticle(Particle &particle);

// put input particles into rep
void updateRep(Particle *particle, myRep &rep, int loopTimes);

//  MOPSOFunction
    //  apply the formulation of PSO
void PSOAdaptionForPhi(Particle &particle, myRep &rep, int it);
void PSOAdaptionForXYZ(Particle &particle, myRep &rep, int it);

    //  decide if particle[i] is dominated, for all i
void decideDominated(Particle * particle);

    //  compare the new position with the pBest and update it
void updatePBest(Particle * particle);

    // calculate f(x)
void getAllParticleCostForTest(Particle * particle);
void getAllParticleCost(Particle * particle);

void *getAParticleCost(void* particle);

    //if the result is true, then cost1 is dominated by cost2
bool isDominated(double *cost1, double *cost2);

    //    add the new non-dominated particles into rep / eliminate some in the rep that is not so good now
void putNewParticleIntoRep(Particle * particle, myRep & rep, int iterator);

    //decide whether particle[i] can be put into rep or not, as well as eliminate the rep particles
    //which are dominated by particle[i]
    // if return true then can put
bool canPutIntoRep(double *Cost, myRep & rep);

void sieveRep(myRep & rep);         //  if there are too many paticles in rep, restrict it in a set value


//  MOPSOAidFunction
    //  cat c1-ss-c2  e.g. :  "mj" , 520 , "forever"  ->  mj520forever        used to conbine file name
char *catStrIntStr(const char *c1, int ss, const char *c2);
char *catStrIntStr(const char *c1, const char *c2, int ss, const char *c3);

    //  cat c1-c2
char *catStrStr(CONCH *c1, CONCH *c2);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7);
char *catStrStr(CONCH *c1, CONCH *c2, CONCH *c3, CONCH *c4, CONCH *c5, CONCH *c6, CONCH *c7, CONCH *c8);

    //  copy p2 to p1
void cpyDoubleArray(double *& p1, double *& p2, int n);

    //  check if V is out of range
void checkV(double *& tmp, int n);

//  check if position is out of range
void checkOrigin(double *&tmp, double *&tmpV, int n);

    //  check if position is out of range
void checkP(double *& tmp, double *& tmpV, int n);

    //  get the right rep_h for this update time
int rouletteWheel(myRep & rep);

    //  get the w
double getInertiaWeight(double it,double MaxIt);

    //  for test object function
double dis1(double x, double y, double z);
double dis2(double x, double y, double z);

    //  get the AA that pr represent
char * getRelativeAA(char pr);

void strReplace(char *src, const char *shortStr, const int start, const int length);

int strToInt(char *str);

int getInputParameter(char **argv, int &nPop, int &MaxIt);

int getMax(int a, int b);

bool doubleEqual(double a, double b);

//  convert rotation to coordinary
void convertRotationToCoordinary(Particle &particle);


//  Matrix Function
    //  ans = a*b
void matrixProduct(double a[4][4], double b[4][4], double ans[4][4]);

    //  set a[4][4] to be full of 0
void initialize_0(double a[4][4]);

    //  dan wei(chinese) matrix
void initialize_1(double a[4][4]);

    //  get the R in the paper
void getR(int it, double * a, double theta, double R[4][4]);

    //  get the new coordinate for itAtom
void getNewCoordinary(int itAtom, double R[4][4], double * a);

void cpyOriginTOAddO_origin(int it, double * x, double * xAddO);

void getNewCoordinaryForO(int it, double R[4][4], double * xAddO, int numAtom);

    //  newR = preR
void cpyMatrix(double preR[4][4], double newR[4][4]);

    //   dan wei hua (chinese)
void unitization(double v[3]);

    //set Q[3][3] the paper used
void getQ(double Q[4][4], double v[3], double theta);


//  get the best reps
void sortRepByAngle(myRep &rep, ANGLE *angle);
void sortRepByLambda(myRep &rep);
void getAngle(int i, const POINT *point, ANGLE *angle);
double vectorDot(VECTOR u, VECTOR v);
double cross(VECTOR u, VECTOR v);
double length(VECTOR u);
bool pdPlusPi(POINT first, POINT second, POINT third);
void angleQsort(int l, int r, ANGLE *a);
void lambdaQsort(int l, int r, LAMBDA *a);
void pointQsort(int l, int r, POINT *a);

// setTime
void setTime();
void getEndTime();

// dispose the input parameters and files previously
void preDisposeInputParametersAndFiles(char **argv);
void disposePDB();
void createSeqTxt();
void createParticleTxt(int k);
void createPhi(int k);
bool isImportantAtom(const Atom atom);
char getAbbreviation(char *str);
void getArgv();

// check similarity
void checkAllParticleSimilarity(Particle *particle);
void *checkOneParticleSimilarity(void *p);
void checkOneParticleSimilarity(Particle &particle);
void runTM_score(int index);                //  run TM_score to get the similarity of 2 particle
double getTM_score(int index);
bool particleSimilarity(double TM_scrore);      // if return 1 then two particles are similar
void becomeInitialParticle(Particle &particle);

//free up space
void freeUpSpace(Particle *& particle);

// apply variables
void applyVariable();

// code for Debug
void printIn(int option);
void printOut(int option);

//  const statement
extern const char rootAddress[];
extern char *inputAddress;
extern char *logAddress;
extern const char *energyFileAddress;
extern const char *tempFileAddress;
extern const char *defaultFileAddress;
extern const char *QUACKoutFileAddress;
extern const char *charmmFileAddress;
extern const char *refine_1Address;
extern const char *TM_scoreAddress;
extern char *answerAddress;
extern const char *draftAddress;
extern const char *scoreAddress;
extern const char *databaseAddress;
extern const char *mybinAddress;
extern const char *strideAddress;

extern const int nVar;
extern const double VarMin;
extern const double VarMax ;
extern const int VarSize[];
extern const double VelMax;             //  =20 without rama_map


// MOPSO Settings
extern int nPop;                  // Population Size
extern const int nRep;                // Repository Size
extern int  MaxIt;          // Maximum Number of Iterations
extern const double Criterion;
extern const int lambdaLoopTimes;
extern const int tidSize;
extern const int answerRepNumber;
extern time_t starTime;
extern time_t endTime;
extern char **argv;

extern const double phi1;
extern const double phi2;
extern const double phi;
extern const double chi ;    // 0.73
extern const int bufferLen;

extern const double wMin;                        //  Inertia Weight
extern const double wMax;
extern const double wDamp;                       //  Inertia Weight Damping Ratio
extern const double c1;                 //  Personal Learning Coefficient
extern const double c2;                 //  Global Learning Coefficient

extern const double Alpha;       //Grid Inflation Parameter
extern const int nGrid;               //Number of Grids per each Dimension
extern const int Beta;                   //Leader Selection Pressure Parameter
extern const int Gamma;             // Extra (to be deleted) Repository Member Selection Pressure

extern const int numObjective;      //  Multiple Objectives
extern const double TM_scoreThreshold; // the threshold of TM-score
extern const double PI;
extern const int multiplyNumber;
extern const double INF;

// my tools
extern Particle *particle;
extern myRep rep;
extern ANGLE *angle;
extern int *sortAns;

//  struct definition
// pBest
struct structBest{
    double *Position;
    double *addO_origin;
    double *Cost;
};

struct Particle{
    double * origin;                          //    the coordinate(x,y,z) of every atom (not include atom O),
    //  three indexes represent an Atom's position
    double * addO_origin;                   // the coordinate(x,y,z) of every atom (include the atom O)
    double * Position;                      //  the angles
    double * old_position;                  //  used in the calculation of angle change
    double * Velocity;                      //  will be used in the MOPSO
    double * Cost ;                         //  will be used int the MOPSO, store every f(x) the particle have
    int * GridIndex, * GridSubIndex;
    bool dominated;                     //   represent if the particle is dominated by some other,
    //  which means it is certainly not the best answer
    int numAA, numAtom;         //  the number of AA and the number of atoms
    int sizeOfOrigin, sizeOfAddO_origin, sizeOfPosition, sizeOfVelocity;    //  the size of arrays
    int index;
    char *seq;
    structBest Best;                    //  pBest
};

struct Rep{
    double *Cost;
    double *addO_origin;
    double *Position;
    int sizeOfAddO_origin;
    int iterator;
};

struct Atom {
    char name[10], group[10];
    int number;
    double x, y, z;
};

struct POINT{
    double x, y;
    int id;
};

struct VECTOR{
    double x, y;
    VECTOR(double x1, double y1){
        x = x1;
        y = y1;
    }
};

struct ANGLE{
    double value;
    int id;
};

struct LAMBDA{
    double value;
    int id;
};

#endif // MOPSO_H
