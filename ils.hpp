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
  void perturbation1(Solution &s);
  void perturbation2(Solution &s);
  void perturbation3(Solution &s);
  void perturbation4(Solution &s);
  bool twoOpt(vector <int> &s);
  bool removeNode(Solution &s);
  bool LS_swap(Solution &s);
  bool LS_swap_outnodes(Solution &s);
  float solver(vector<int> &y, float currentCost, float timeLimit, Solution &s, bool flag, bool initSol, bool original, int tries, bool fixed);
  float coveringSolver(Solution &s, float add);
public:
  default_random_engine generator;
  uniform_real_distribution<double> randValue;
  // function
  ILS(int seed, int N, int D, int We, int W, int S, int B, int **bo, int ND, float T, int V, float M,
      vector<vector<int> > &b, vector<vector<int> > &bi,  float **t, int *mMax, int *mMin, float **C, float costUAV);
  virtual ~ILS();
  void run(float *result);
  void runMH(float *result);

};
#endif /* _ILS_H_ */
