/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#ifndef _SOLUTION_H_
#define _SOLUTION_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <iomanip>
#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

#define Infinity 99999999
using namespace std;
ILOSTLBEGIN

class Solution {
public:
  Solution();
  Solution(const Solution& receive);
  Solution(int N, int V, int WD);

  vector<int> Sd;
  vector<vector<int> > Ss;
  vector<bool> noInSs;
  vector<int> coveredArea;
  int feasible;
  float Lcost;
  float Rcost;
  float cost;
};

#endif /* _SOLUTION_H_ */
