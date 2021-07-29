/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "solution.hpp"

Solution::Solution()
{
}

Solution::Solution(const Solution& receive)
{
  this->feasible = receive.feasible;
  this->Lcost = receive.Lcost;
  this->Rcost = receive.Rcost;
  this->cost = receive.cost;

  this->Sd = receive.Sd;
  this->coveredArea = receive.coveredArea;
  for(unsigned i = 0; i < receive.Ss.size(); i++)
    {
      this->Ss.push_back(vector<int>());
      for(unsigned j = 0; j < receive.Ss[i].size(); j++)
	this->Ss[i].push_back(receive.Ss[i][j]);
    }

  this->noInSs = receive.noInSs;
}

Solution::Solution(int N, int V, int WD)
{
  // initialize depot
  Sd = vector<int>(WD);

  // initialize areas to cover
  coveredArea = vector<int>(N, 0);
  for(auto i = 0; i < N; i++)
    coveredArea.push_back(0);

  // initialize routes
  for(auto d = 0; d < WD * 2; d++)
    Ss.push_back(vector<int>());

  // initialize noInSs
  noInSs = vector<bool>(V);
  //noInSs.reserve(V);
}
