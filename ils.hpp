/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#ifndef _ILS_H_
#define _ILS_H

#include "solution.hpp"
#define Infinity 99999999
using namespace std;

class ILS {
private:
  // Solution solution;
  map<float, Solution> elite;
  int N;
  int D;
  int We;
  int W;
  int WD;
  int S;
  int B;
  int ND;
  float M;
  float T;
  int V;
  float costUAV;
  float **t;
  int *mMax;
  int *mMin;
  float **C;
  int **b;
  int seed;
  vector<vector<int> > beta;
  vector<vector<int> > bi;
  unsigned decipresi;
  double timeCplex;

  // Parameters
  unsigned E;  // size of elite solutions
  unsigned Kt; // iterations number
  unsigned Km; // number that every Km iterations calls to the MILP solver
  unsigned Kg; // number that every Kg iterations calls to the global reset
  unsigned Kl; // number that every Kl iterations calls to the local reset
  float p1;    // probability for perturbation 1
  float p2;    // probability for perturbation 2
  float p3;    // probability for perturbation 3
  float p4;    // probability for perturbation 4
  float p5;    // probability for perturbation 5
  float wp;     // probability to fix variables to zero in the MILP solver
  float Tc;    // the time limit to solve by Cplex

  // auxiliar function
  int rnd(unsigned low, unsigned high);
  void print_sol(Solution &s);
  void computeRoute(Solution &s, vector<float> &c);
  float computeCost(Solution &s);
  void coveringAreaByNode(Solution &s, int i);
  void copyS(Solution &s, Solution &sbest);
  void updateElite(Solution &s, float cost);

  // operator
  bool initialSolution(Solution &s);
  void perturbation2(Solution &s);
  void perturbation3(Solution &s);
  void perturbation4(Solution &s);
  void perturbation5(Solution &s);
  bool twoOpt(vector <int> &s);
  bool removeNode(Solution &s);
  bool LS_swap(Solution &s);
  bool LS_swap_outnodes(Solution &s);
  float solver(vector<int> &y, float currentCost, float timeLimit, Solution &s, bool flag, bool initSol, bool original, int tries, bool fixed);
  float coveringSolver(Solution &s, float add);
public:
  default_random_engine generator;
  uniform_real_distribution<double> randValue;
  string nameAlgorithm;
  // results data
  float initialCost;
  float costFinal;
  float timeF;
  // function
  ILS(int seed, int N, int D, int We, int W, int S, int B, int **bo, int ND, float T, int V, float M,
      vector<vector<int> > &b, vector<vector<int> > &bi,  float **t, int *mMax, int *mMin, float **C, float costUAV, float *parameters);
  virtual ~ILS();
  void run(int timeLimitAlgorithm, int variant);
  void runMH(int timeLimitAlgorithm);
  void printOutput(string Instance);
};
#endif /* _ILS_H_ */
