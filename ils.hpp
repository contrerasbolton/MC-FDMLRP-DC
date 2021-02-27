/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#ifndef _ILS_H_
#define _ILS_H

#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <algorithm>
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN

class ILS {
private:
  struct Solution {
    int *Sd;
    vector<vector<int> > Ss;
    vector<int> noInSs;
    int *coveredArea;
    int feasible;
    float Lcost;
    float Rcost;
    float cost;
  };
  Solution solution;
  int N;
  int D;
  int K;
  int S;
  int B;
  int ND;
  float T;
  int V;
  float **t;
  int *mMax;
  int *mMin;
  int **C;
  int seed;
  vector<vector<vector<int> > > alpha;
  vector<vector<int> > beta;
  vector<vector<int> > bi;
  // function
  int rnd(unsigned low, unsigned high);
  void print_sol(Solution &s);
  float computeCost(Solution &s);
  void coveringAreaByDepot(Solution &s, int d);
  void coveringAreaByNode(Solution &s, int i);
  void initialSolution(Solution &s);
  void perturbation1(Solution &s);
  bool twoOpt(vector <int> &s);
  bool removeNode(Solution &s);
public:
  default_random_engine generator;
  uniform_real_distribution<double> randValue;
  // function
  ILS(int seed, int N, int D, int K, int S, int B, int ND, float T, int V,
      vector<vector<vector<int> > > &a, vector<vector<int> > &b, vector<vector<int> > &bi,  float **t, int *mMax, int *mMin, int **C);
  virtual ~ILS();
  void run();

};
#endif /* _ILS_H_ */
