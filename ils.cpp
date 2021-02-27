/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "ils.hpp"

ILS::ILS(int seed, int N, int D, int K, int S, int B, int ND, float T, int V,
         vector<vector<vector<int> > > &a, vector<vector<int> > &b, vector<vector<int> > &bi, float **t, int *mMax, int *mMin, int **C)
{
  this->seed = seed;
  this->N = N;
  this->D = D;
  this->K = K;
  this->S = S;
  this->B = B;
  this->ND = ND;
  this->T = T;
  this->V = V;
  this->alpha = a;
  this->beta = b;
  this->bi = bi;
  this->t = t;
  this->mMax = mMax;
  this->mMin = mMin;
  this->C = C;

  // Initialize
  generator.seed(seed);

  // solution.Sd = new int[D];
  // for(auto i = 0; i < D; i++)
  //   solution.Sd[i] = 1;
}

ILS::~ILS()
{

}

void ILS::initialSolution(Solution &s)
{
  // s.Sd[1] = 0;
  // nodes {D, D + 1, ..., D + S - 1, D + S}
  vector<int> nodes;
  for(auto i = 0; i < S; i++)
    nodes.push_back(i + D);

  // random nodes
  shuffle(nodes.begin(), nodes.end(), default_random_engine(seed));

  // reset routes
  for(auto c = 0; c < D * 2; c++)
    s.Ss[c].clear();

  // for(unsigned i = 0; i < nodes.size(); i++)
  //   cout << nodes[i] << " ";
  // cout << endl;

  float currentCost = 0.0;
  int u = 0, v = 0, next = 0;

  int lastPosition = 0;
  for(auto c = 0; c < D * 2; )
    {
      //cout << c << endl;
      int d = c / 2;
      vector<int> aux;
      unsigned i = next;
      float routeCost = 0.0;
      if(s.Sd[d] == 1)
        {
          aux.push_back(d);
          u = nodes[i];
          routeCost += t[d][u];
          //cout << "("<< d << " " << u << ") ";
          for(; i < nodes.size() - 1; i++)
            {
              u = nodes[i];
              v = nodes[i + 1];
              currentCost = t[u][v];
              if(routeCost + currentCost <= T)
                {
                  routeCost += t[u][v];
                  coveringAreaByNode(s, u);
                  //cout << "("<< u << " " << v << ") ";
                  aux.push_back(u);
                }
              else // returning to depot
                {
                  u = nodes[i];
                  //cout << "("<< u << " " << d << ") " << endl;
                  next = i + 1;
                  coveringAreaByNode(s, u);
                  routeCost += t[u][d];
                  //s.Rcost += routeCost;
                  aux.push_back(u);
                  aux.push_back(d);
                  s.Ss[c] = aux;
                  c++;
                  break;
                }
            }
          // last route finished
          if(i == nodes.size() - 1)
            {
              u = nodes[i];
              //cout << "("<< u << " " << d << ") " << endl;
              coveringAreaByNode(s, u);
              routeCost += t[u][d];
              aux.push_back(u);
              aux.push_back(d);
              s.Ss[c] = aux;
              //s.Rcost += routeCost;
              break;
            }
        }
      else
        c += 2;
      lastPosition = i;
    }
  //cout << "last " << lastPosition << " " << S << endl;

  if(lastPosition < S)
    {
      s.noInSs.clear();
      for(auto i = lastPosition; i < S; i++)
        s.noInSs.push_back(nodes[i]);
    }

}

void ILS::run()
{
  cout << "ILS" << endl;
  Solution s;

  // depot
  s.Sd = new int[D];
  for(auto i = 0; i < D; i++)
    s.Sd[i] = rnd(0, 2);

  // areas to cover
  s.coveredArea = new int [N];
  for(auto i = 0; i < N; i++)
    s.coveredArea[i] = 0;

  // initialize routes
  for(auto d = 0; d < D * 2; d++)
    s.Ss.push_back(vector<int>());

   // build initial solution
   initialSolution(s);

  float s_cost = computeCost(s);
  cout << "cost " << s_cost << endl;
  print_sol(s);

  Solution sBest;
  float s_bestCost = s_cost;
  for(auto iter = 0; iter < 1000; iter++)
    {
      if(iter % 5 == 0)
        initialSolution(s);

      // 2-opt for each route
      for(auto c = 0; c < D * 2; c++)
        if(s.Ss[c].size())
          twoOpt(s.Ss[2]);

      // s_cost = computeCost(s);
      // cout << "cost " << s_cost << endl;
      // print_sol(s);

      removeNode(s);
      s_cost = computeCost(s);
      cout << "iter " << iter << ": "<< s_cost << endl;
      if(s_bestCost > s_cost)
        {
          s_bestCost = s_cost;
          sBest = s;
        }
    }
  s_cost = computeCost(sBest);
  cout << "cost " << s_bestCost << endl;
  print_sol(sBest);
}

void ILS::coveringAreaByDepot(Solution &s, int d)
{
  int cover = -1;
  for(unsigned k = 0; k < alpha[d].size(); k++)
    {
      for(unsigned i = 0; i < alpha[d][k].size(); i++)
        {
          cover = alpha[d][k][i];
          s.coveredArea[cover]++;
        }
    }
}

void ILS::coveringAreaByNode(Solution &s, int i)
{
  int cover = -1;
  for(unsigned k = 0; k < beta[i].size(); k++)
    {
      cover = beta[i][k];
      s.coveredArea[cover]++;
    }
}

void ILS::perturbation1(Solution &s)
{
  // int count2 = 0;
  // for(auto d = 0; d < D; d++)
  //   {
  //     if(s1.Sd[i])
  //       count1++;
  //     else if(s1.Sd[i] == 2)
  //       count2++;
  //   }
  // if(count1
}

bool ILS::removeNode(Solution &s)
{
  int node = 0, area = 0;
  int route = 0;
  do {
    route = rnd(0, s.Ss.size() - 1);
  } while(s.Ss[route].size() == 0);

  for(unsigned i = 1; i < s.Ss[route].size() - 1; i++)
    {
      node = s.Ss[route][i];
      bool flag = true;
      for(unsigned j = 0; j < beta[node].size(); j++)
        {
          area = beta[node][j];
          if(s.coveredArea[area] < 2)
            {
              flag = false;
              break;
            }
        }
      if(flag)
        {
          s.Ss[route].erase(s.Ss[route].begin() + i); // remove node
          s.noInSs.push_back(node); // add to list noInSs
          for(unsigned j = 0; j < beta[node].size(); j++) // decreasing covered areas
            s.coveredArea[beta[node][j]]--;
          i--;
        }
    }
  // if all the nodes are removed so remove this route
  if(s.Ss[route].size() == 2)
    {
      s.Ss[route].pop_back();
      s.Ss[route].pop_back();
      // if(route % 2 == 0)
      //   {
      //     if(s.Ss[route + 1].size() == 0)
    }
  return true;
}

bool ILS::twoOpt(vector <int> &s)
{
  //cout << "size " << s.size() << endl;
  if(s.size() < 3)
    return false;
  float diff = 0;
  bool improve = false;
  int depot = s[0];

  // for(unsigned i = 0; i < s.size(); i++)
  //   cout << s[i] << " ";
  // cout << endl;

  for(unsigned i = 0; i < s.size(); i++)
    {
      for(unsigned j = i + 2; j < s.size() - 1; j++)
        {
          // newEdges ((i, j) + (i + 1, j + 1)) - currentEdges ((i, i + 1) + (j, j + 1))
          diff = t[s[i]][s[j]] + t[s[i + 1]][s[j + 1]] - t[s[i]][s[i + 1]] - t[s[j]][s[j + 1]];
           if(diff < 0)
            {
              reverse(s.begin() + i + 1, s.begin() + j + 1);
              improve = true;
            }
        }
    }
  if(depot != s[0])
    {
      cout << "must fix the depot at begin" << endl;
      exit(0);
    }

  if(improve)
    return true;
  else
     return false;
}

float ILS::computeCost(Solution &s)
{
  // reset covered Area
  for(auto i = 0; i < N; i++)
    s.coveredArea[i] = 0;

  // compute location cost
  s.Lcost = 0.0;
  for(auto i = 0; i < D; i++)
    {
      if(s.Sd[i] != 0)
        {
          s.Lcost += C[i][s.Sd[i] - 1];
          coveringAreaByDepot(s, i);
        }
    }
  //cout << s.Lcost << endl;

  // compute drones' routes
  int u = 0, v = 0;
  float routeCost = 0.0;
  s.Rcost = 0.0;
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size())
        {
          for(unsigned j = 0; j < s.Ss[i].size() - 2; j++)
            {
              u = s.Ss[i][j];
              v = s.Ss[i][j + 1];
              //cout << u << " " << v << endl;
              routeCost += t[u][v];
              coveringAreaByNode(s, v);
            }
          u = s.Ss[i][s.Ss[i].size() - 2];
          v = s.Ss[i][s.Ss[i].size() - 1];
          routeCost += t[u][v];
        }
    }

  s.Rcost = routeCost;
  //cout << s.Rcost << endl;
  // final cost
  s.cost = s.Lcost + s.Rcost;
  // check if the solution is feasible
  s.feasible = 1;
  for(auto i = 0; i < N; i++)
    if(!s.coveredArea[i])
      s.feasible = 0;
  //cout << s.cost << endl;

  return s.cost;
}

int ILS::rnd(unsigned low, unsigned high)
{
  return low + generator() % ((high + 1) - low);
}

void ILS::print_sol(Solution &s)
{
  cout << "Sd:" << endl;
  for(auto i = 0; i < D; i++)
    cout << s.Sd[i] << " ";
  cout << endl;
  cout << "Ss:" << endl;
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(i % 2 == 0)
        cout << "Depot " << i / 2 << ": "<< endl;

      cout << "\tRoute_" << i % 2 << ": ";
      for(unsigned j = 0; j < s.Ss[i].size(); j++)
        {
          cout << s.Ss[i][j] << " ";
        }
      cout << endl;
    }

  cout << "noInSs = ";
  for(unsigned i = 0; i < s.noInSs.size(); i++)
    cout << s.noInSs[i] << " ";
  cout << endl;

  cout << "cA = ";
  for(auto i = 0; i < N; i++)
      cout << s.coveredArea[i] << " ";
  cout << endl;

  cout << "Feasible = " << s.feasible << endl;

}
