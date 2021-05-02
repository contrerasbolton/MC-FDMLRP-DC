/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "ils.hpp"

unsigned sizeElite = 10;

ILS::ILS(int seed, int N, int D, int K, int S, int B, int **beta, int ND, float T, int V, float M,
         vector<vector<vector<int> > > &a, vector<vector<int> > &b, vector<vector<int> > &bi, float **t, int *mMax, int *mMin, float **C)
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
  this->M = M;
  this->b = beta;
  // Initialize
  generator.seed(seed);
  randValue = uniform_real_distribution<double>(0.0, 1.0);
  // solution.Sd = new int[D];
  // for(auto i = 0; i < D; i++)
  //   solution.Sd[i] = 1;
}

ILS::~ILS()
{

}

ILOMIPINFOCALLBACK5(timeLimitCallback,
                    IloCplex, cplex,
                    IloBool,  aborted,
                    IloNum,   timeStart,
                    IloNum,   timeLimit,
                    IloNum,   UB)
{

  //out << "oli " << UB << " " << getBestObjValue() << " " << getIncumbentObjValue() << endl;
  if(getBestObjValue() > UB)
    {
      getEnv().out() << endl
                     << " it is not explored because LB > current UB." << endl;
      aborted = IloTrue;
      abort();
    }
  if ( !aborted  &&  hasIncumbent() ) {
    IloNum gap = 100.0 * getMIPRelativeGap();
    IloNum timeUsed = cplex.getCplexTime() - timeStart;



    // if ( timeUsed > 1 )
    //    getEnv().out() << timeUsed << endl;
    if ( timeUsed > timeLimit // && gap < acceptableGap
         ) {
      getEnv().out() << endl
                     << "Good enough solution at "
                     << timeUsed << " sec., gap = "
                     << gap << "%, quitting." << endl;
      aborted = IloTrue;
      abort();
    }
  }
}

bool ILS::initialSolution(Solution &s)
{
  // cout << "Start Initial Solution" << endl;
  // reset covered Area
  for(auto i = 0; i < N; i++)
    s.coveredArea[i] = 0;

  // reset depot
  int maxK[2] = {0, 0};
  for(auto i = 0; i < D; i++)
    s.Sd[i] = -1;

  // assures the minimum K = 1
  int d = 0;

  // assures the minimum K = 0
  for(auto i = 0; i < mMin[0]; i++)
    {
      do {
        d = rnd(0, D - 1);
      } while(s.Sd[d] != -1);
      s.Sd[d] = 0;
    }
  maxK[0]++;
  // assures the minimum K = 1
  for(auto i = 0; i < mMin[1]; i++)
    {
      do {
        d = rnd(0, D - 1);
      } while(s.Sd[d] != -1);
      s.Sd[d] = 1;
    }
  maxK[1]++;

  for(auto i = 0; i < D; i++)
    {
      if(s.Sd[i] == -1)
        {
          s.Sd[i] = rnd(-1, K - 1);
          if(maxK[s.Sd[i]] <= mMax[0] && maxK[s.Sd[i]] <= mMax[1])
            {
              if(s.Sd[i] != -1)
                maxK[s.Sd[i]]++;
            }
          else if(maxK[s.Sd[i]] >= mMax[0] && maxK[s.Sd[i]] <= mMax[1])
            {
              s.Sd[i] = rnd(-1, K - 2);
              if(s.Sd[i] == 0)
                maxK[s.Sd[i]]++;
            }
          else if(s.Sd[i] >= mMax[1] && maxK[s.Sd[i]] <= mMax[0])
            {
              s.Sd[i] = rnd(-1, K - 2);
              maxK[s.Sd[i]]++;
            }
          else
            s.Sd[i] = -1;
        }
    }

  // cout << "Sd" << endl;
  for(auto i = 0; i < D; i++)
    {
      // cout << s.Sd[i] << " ";
      if(s.Sd[i] > -1)
        coveringAreaByDepot(s, i, s.Sd[i]);
    }
  // reset routes
  for(auto c = 0; c < D * 2; c++)
    s.Ss[c].clear();

  // cout << "only depot" << endl;

  // for(auto i = 0; i < N; i++)
  //   cout << s.coveredArea[i] << " ";
  // cout << endl;

  vector<int> nodes;
  //ector<bool> mark(V, false);
  for(auto i = 0; i < V; i++)
    s.noInSs[i] = false;

  // Nodes subset according to uncovered areas
  for(auto i = 0; i < N; i++)
    {
      if(!s.coveredArea[i])
        {
          for(unsigned j = 0; j < bi[i].size(); j++)
            {
              int node = bi[i][j];
              if(!s.noInSs[node])
                {
                  s.noInSs[node] = true;
                  nodes.push_back(node);
                }
            }
        }
    }

  // random nodes
  shuffle(nodes.begin(), nodes.end(), default_random_engine(seed));

  // for(unsigned i = 0; i < nodes.size(); i++)
  //   cout << nodes[i] << " ";
  // cout << endl;
  // exit(0);
  // float currentCost = 0.0;
  int u = 0, v = 0;
  for(auto c = 0; c < D * 2; )
    {
      d = c / 2;
      vector<int> aux;
      float routeCost = 0.0;
      //cout << "Test in Sd: " << d << ": " << c << " = " << s.Sd[d] << endl;
      if(s.Sd[d] == 0 && nodes.size())
        {
          aux.push_back(d);
          u = nodes.back();
          routeCost += t[d][u];
          //cout << "("<< d << " " << u << ")* ";
          while(nodes.size())
            {
              u = nodes.back();
              v = nodes[nodes.size() - 2];
              if(routeCost + t[u][v] <= T)
                {
                  routeCost += t[u][v];
                  coveringAreaByNode(s, u);
                  //cout << "("<< u << " " << v << ") ";
                  aux.push_back(u);
                  nodes.pop_back();

                  if(!nodes.size()) // there no nodes
                    {
                      //cout << "what" << endl;
                      routeCost += t[u][d];
                      aux.push_back(d);
                      s.Ss[c] = aux;
                      break;
                    }
                }
              else // returning to depot
                {
                  u = aux.back();
                  //cout << "var " << routeCost + t[u][d] << endl;
                  if(routeCost + t[u][d] <= T) // returning from last u
                    {
                      //cout << "("<< u << " " << d << ")_ " << endl;

                      aux.push_back(d);
                      routeCost += t[u][d];
                      s.Ss[c] = aux;
                    }
                  else //cannot returning from last u
                    {
                      if(aux.size() > 2) // at least a depot and two node, so remove the last node
                        {
                          //cout << "("<< u << " " << d << ")_ " << endl;
                          coveringAreaByNode(s, u);
                          aux.push_back(d);
                          routeCost += t[u][d];
                          s.Ss[c] = aux;
                        }
                      else // in this case the route is removed
                        {
                          //cout << "remove" << endl;
                          if(u >= D)
                            nodes.push_back(u);
                          // remove route
                          aux.pop_back();
                          aux.pop_back();
                        }
                    }
                  c++;
                  break;
                }
            }
        }
      else
        c += 2;
    }
  //cout << endl;

  // cout << "Ss:" << endl;
  // for(unsigned i = 0; i < s.Ss.size(); i++)
  //   {
  //      if(i % 2 == 0)
  //        cout << "Depot " << i / 2 << ": "<< endl;

  //     cout << "\tRoute_" << i % 2 << ": ";
  //     if(s.Ss[i].size())
  //       {
  //         float temp = 0.0;
  //         for(unsigned j = 0; j < s.Ss[i].size() - 1; j++)
  //           {
  //              cout << s.Ss[i][j] << " ";
  //             //temp += t[s.Ss[i][j]][s.Ss[i][j + 1]];
  //           }
  //         cout << s.Ss[i][s.Ss[i].size() - 1] << " ";
  //         // cout << ": " << temp << endl;
  //       }
  //     cout << endl;
  //   }
  // cout << "quedan:" << endl;
  // for(unsigned i = 0; i < nodes.size(); i++)
  //   cout << nodes[i] << " ";
  // // exit(0);

  if(nodes.size())
    {
      //cout << "quedan (" << nodes.size() << ")"<< endl;
      // check if the solution is feasible
      for(auto ii = 0; ii < N; ii++)
        {
          if(!s.coveredArea[ii])
            {
              s.feasible = 0;
              // for(unsigned i = 0; i < nodes.size(); i++)
              //   cout << nodes[i] << " ";
              // cout << endl;
              // for(auto i = 0; i < N; i++)
              //   cout << s.coveredArea[i] << " ";
              // cout << endl;

              for(auto i = nodes.size() - 1; nodes.size(); i--)
                {
                  int u, v, w;
                  float minimum = 9999999.9;
                  int here = -1, route = 0;
                  w = nodes.back();
                  //cout << "arreglar " << w << endl;
                  for(auto c = 0; c < D * 2; c++)
                    {
                      d = c / 2;
                      // cout << "r" << c << endl;
                      if(s.Sd[d] == 0)
                        {
                          if(s.Ss[c].size() > 0)
                            {
                              for(unsigned r = 0; r < s.Ss[c].size() - 1; r++)
                                {
                                  u = s.Ss[c][r];
                                  v = s.Ss[c][r + 1];

                                  // cout << u << " " << v << "-> ";
                                  float newCost = t[u][w] + t[w][v] - t[u][v];
                                  // cout << newCost << " " << endl;
                                  if(newCost < minimum)
                                    {
                                      //cout << u << " " << v << "-> " << newCost << " " << endl;
                                      newCost = minimum;
                                      here = r;
                                      route = c;
                                    }

                                }
                            }
                          else
                            {
                              //cout << "no había y se agrega: " << d << " " << w << " " << d<< endl;
                              s.Ss[c].push_back(d);
                              s.Ss[c].push_back(w);
                              coveringAreaByNode(s, w);
                              s.Ss[c].push_back(d);
                              nodes.pop_back();
                              w = nodes.back();
                              here = -1;
                              break;
                            }
                        }
                    }
                  if(here != -1)
                    {
                      //cout << "agrego " << w << endl;

                      s.Ss[route].insert(s.Ss[route].begin() + here + 1, w);
                      coveringAreaByNode(s, w);
                      nodes.pop_back();
                      ii = 0;
                      break;
                    }
                }
              // test = false;
            }
        }
    }

  for(unsigned i = 0; i < nodes.size(); i++)
    s.noInSs[nodes[i]] = false;

  // for(auto i = D; i < V; i++)
  //   {
  //     s.noInSs[i] = false;
  //     if(mark[i])
  //       s.noInSs[i]
  //   }

  // if(lastPosition < S)
  //   {
  //     s.noInSs.clear();
  //     for(auto i = lastPosition; i < S; i++)
  //       s.noInSs.push_back(nodes[i]);
  //   }
  // cout << "aqui en initial solution" << endl;

  // for(auto i = 0; i < N; i++)
  //   cout << s.coveredArea[i] << " ";
  // cout << endl;

  //cout << computeCost(s) << endl;
  computeCost(s);
  // print_sol(s);
  // cout << "End initial Solution" << endl;

  // exit(0);
  return true;
}

void ILS::copyS(Solution &give, Solution &receive)
{
  for(auto i = 0; i < D; i++)
    receive.Sd[i] = give.Sd[i];

  for(auto i = 0; i < N; i++)
    receive.coveredArea[i] = give.coveredArea[i];

  for(auto i = 0; i < V; i++)
    receive.noInSs[i] = give.noInSs[i];

  for(unsigned i = 0; i < give.Ss.size(); i++)
    {
      receive.Ss[i].clear();
      for(unsigned j = 0; j < give.Ss[i].size(); j++)
        receive.Ss[i].push_back(give.Ss[i][j]);
    }
}

float ILS::solver(vector<int> &y, float currentUB, float timeLimit, Solution &s, bool flag, bool initSol)
{
  if(initSol)
    {
      vector<int> v(2, 0);

      for(auto i = 0; i < D; i++)
        {
          y[i] = s.Sd[i];
          if(y[i] != -1)
            v[y[i]]++;
        }

      int d = rnd(0, D - 1);

      if(v[0] < mMax[0])
        {
          if(y[d] != 0)
            {
              if(y[d] != -1)
                {
                  if(v[1] > mMin[1])
                    {
                      y[d] = 0;
                    }
                }
            }
          else if(y[d] == -1)
            {
              y[d] = 0;
            }
          else if(y[d] == 0)
            {
              if(v[0] > mMin[0])
                {
                  if(v[1] < mMax[1])
                    y[d] = 1;
                  else
                    y[d] = -1;
                }
            }
        }
      // cout << "******new**********************" << endl;
      // for(auto i = 0; i < D; i++)
      //   {
      //    cout << y[i] << " ";
      //   }
    }
  IloEnv env;
  try
    {
      IloModel model1(env);
      stringstream name;

      vector<int> L(N, 0);
      // Definition of the decision variables
      // decision variables P
      IloBoolVarArray P(env, N);
      for(auto i = 0; i < N; ++i)
        {
          name << "P_" << i;
          P[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
          name.str("");
        }

      // decision variables v
      IloBoolVarArray v(env, S);
      for(auto i = 0; i < S; i++)
        {
          name << "v_" << i;
          v[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
          name.str("");
        }

      // decision variables u
      IloNumVarArray u(env, V);
      for(auto i = 0; i < V; i++)
        {
          name << "u_" << i;
          u[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, name.str().c_str());
          name.str("");
        }

      // decision variables w
      IloNumVarArray w(env, V);
      for(auto i = 0; i < V; i++)
        {
          name << "w_" << i;
          w[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Int, name.str().c_str());
          name.str("");
        }

      // decision variables x
      IloArray<IloNumVarArray> x(env, V);
      for(auto i = 0; i < V; i++)
        {
          x[i] = IloNumVarArray(env, V);
          for(auto j = 0; j < V; j++)
            {
              name << "x_" << i << "_" << j;
              x[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
              name.str("");
            }
        }
      for(auto d = 0; d < D; d++)
        {
          for(auto k = 0; k < K; k++)
            {
              if(y[d] == k)
                {
                  for(unsigned i = 0; i < alpha[d][k].size(); i++)
                    {
                      int cover = alpha[d][k][i];
                      L[cover] = 1;
                    }
                }
            }
        }

      //cout << "sum " << M << endl;
      M = 50;
      // cout << "L: ";
      // for(auto i = 0; i < N; i++)
      //   cout << L[i] << " ";
      // cout << endl;
      // exit(0);

      // Second stage
      IloModel model2(env);
      IloCplex cplex2(model2);
      // // Function objetive: 1) \Sum_{i in N} P_i
      // IloExpr fo2(env);
      // for(auto i = 0; i < N; i++)
      //   fo2 += P[i];
      // model2.add(IloMaximize(env, fo2));
      // Another function objective: min route cost
      IloExpr fo2(env);
      for(auto i = 0; i < V; i++)
        for(auto j = 0; j < V; j++)
          if(i != j)
            fo2 += t[i][j] * x[i][j];

      for(auto d = 0; d < D; d++)
        if(y[d] > -1)
          fo2 += C[d][y[d]];
      model2.add(IloMinimize(env, fo2));

      // for(auto i = 0; i < V; i++)
      //   {
      //     model2.add(x[0][i] == 0);
      //     model2.add(x[i][0] == 0);
      //     model2.add(x[2][i] == 0);
      //     model2.add(x[i][2] == 0);

      //   }

      // original function objective is became to constraint
      // Constraints 3) \Sum_{i in N} L_i + L_i == N
      IloRangeArray constraint3(env, 1);
      IloExpr exp(env);
      for(auto i = 0; i < N; i++)
        exp += L[i] + P[i];
      name << "constraint3";
      constraint3[0] = IloRange(env, N, exp, N, name.str().c_str());
      name.str("");
      model2.add(constraint3);

      // Drone routing constraints
      // Constraints 7) \Sum_{k in K} y_{dk} <= M_k, for all d in D
      IloRangeArray constraint7(env, D);
      for(auto d = 0; d < D; d++)
      	{
      	  IloExpr exp(env);
          for(auto j = 0; j < V; j++)
            exp += x[d][j];
          //exp -= ND * y[d][0];
          if(y[d] == 0)
            exp -= ND * 1;
      	  name << "constraint7_" << d;
      	  constraint7[d] = IloRange(env, -IloInfinity, exp, 0, name.str().c_str());
      	  name.str("");
      	}
      model2.add(constraint7);

      // Constraints 8) \Sum_{j in V} x_{dj} - \Sum_{j in V} x_{jd} = 0, for all d in D
      IloRangeArray constraint8(env, D);
      for(auto d = 0; d < D; d++)
        {
          IloExpr sum1(env);
          for(auto j = 0; j < V; j++)
            sum1 += x[j][d];

          IloExpr sum2(env);
          for(auto j = 0; j < V; j++)
            sum2 += x[d][j];

          name << "constraint8_" << d;
      	  constraint8[d] = IloRange(env, 0, sum2 - sum1, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint8);

      // Constraints 9) \Sum_{j in V} x_{ij} = v_i, for all i in S
      IloRangeArray constraint9(env, S);
      for(auto i = D, c = 0; i < V; i++, c++)
        {
          IloExpr sum(env);
          for(auto j = 0; j < V; j++)
            sum += x[i][j];
          sum -= v[c];
          name << "constraint9_" << i;
      	  constraint9[c] = IloRange(env, 0, sum, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint9);

      // Constraints 10) \Sum_{j in V} x_{ji} = v_i, for all i in S
      IloRangeArray constraint10(env, S);
      for(auto i = D, c = 0; i < V; i++, c++)
        {
          IloExpr sum(env);
          for(auto j = 0; j < V; j++)
            sum += x[j][i];
          sum -= v[c];
          name << "constraint10_" << i;
      	  constraint10[c] = IloRange(env, 0, sum, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint10);

      // decision variables g
      IloArray<IloNumVarArray> g(env, V);

      for(auto i = 0; i < V; i++)
        {
          g[i] = IloNumVarArray(env, V);
          for(auto j = 0; j < V; j++)
            if(i != j)
              {
                name << "g_" << i << "_" << j;
                g[i][j] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, name.str().c_str());
                name.str("");
              }
        }

      // Constraints 11) \Sum_{j in V} g_{ij} - \Sum_{j in V} g_{ji} == \Sum_{j in V} t_{ij} * x_{ij}, for all i in S,x i != j
      IloRangeArray  constraint11(env, S);
      for(auto i = D, c = 0; i < V; i++, c++)
        {
          IloExpr sum1(env);
          for(auto j = 0; j < V; j++)
            if(i != j)
              sum1 += g[i][j];

          IloExpr sum2(env);
          for(auto j = 0; j < V; j++)
            if(i != j)
              sum2 += g[j][i];

          IloExpr sum3(env);
          for(auto j = 0; j < V; j++)
            if(i != j)
              sum3 += t[i][j] * x[i][j];

          name << "constraint11_" << i;
          constraint11[c] = IloRange(env, 0, sum1 - sum2 - sum3, 0, name.str().c_str());
          name.str("");
        }
      model2.add(constraint11);

      IloArray<IloRangeArray> constraint12a(env, V);
      for(auto i = D; i < V; i++)
        {
          constraint12a[i] = IloRangeArray(env, V);
          for(auto j = D; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint12_UB" << i << "_" << j;
                  constraint12a[i][j] = IloRange(env, -IloInfinity, g[i][j] - T * x[i][j], 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint12a[i][j]);
                }
            }
        }
      IloArray<IloRangeArray> constraint12b(env, V);
      for(auto i = D; i < V; i++)
        {
          constraint12b[i] = IloRangeArray(env, V);
          for(auto j = D; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint12_LB_" << i << "_" << j;
                  constraint12b[i][j] = IloRange(env, 0, g[i][j] - x[i][j] * t[i][j], IloInfinity, name.str().c_str());
                  name.str("");
                  model2.add(constraint12b[i][j]);
                }
            }
        }

      IloArray<IloRangeArray> constraint13(env, D);
      for(auto i = 0; i < D; i++)
        {
          constraint13[i] = IloRangeArray(env, V);
          for(auto j = 0; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint13_LB_" << i << "_" << j;
                  constraint13[i][j] = IloRange(env, 0, g[i][j] - x[i][j] * t[i][j], 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint13[i][j]);
                }
            }
        }

      IloArray<IloRangeArray> constraint13b(env, V);
      for(auto i = 0; i < V; i++)
        {
          constraint13b[i] = IloRangeArray(env, D);
          for(auto j = 0; j < D; j++)
            {
              if(i != j)
                {
                  name << "constraint13_UB_" << i << "_" << j;
                  constraint13b[i][j] = IloRange(env, -IloInfinity, g[i][j] - x[i][j] * T, 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint13b[i][j]);
                }
            }
        }

      // Constraints 14) \Sum_{j in S} beta_{ji} * v_j >= P_i
      IloRangeArray constraint14(env, N);
      for(auto i = 0; i < N; i++)
        {
          IloExpr exp(env);
          for(auto j = 0; j < S; j++)
            exp += b[j][i] * v[j];
          exp -= P[i];
          name << "constraint14_" << i;
          constraint14[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model2.add(constraint14);

      // Constraints 15) u_i = d, for all i in D
      IloRangeArray constraint15(env, D);
      for(auto d = 0; d < D; d++)
        {
          name << "constraint15_" << d;
      	  constraint15[d] = IloRange(env, d, w[d], d, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint15);

      // Constraints 16) w_i +- w_j < = (D - 1)  * (1 - x_{ji}), for all i in S, j in V, i != j
      IloArray<IloRangeArray> constraint16(env, V);
      for(auto i = 0; i < V; i++)
        {
          constraint16[i] = IloRangeArray(env, V);
          for(auto j = 0; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint16_" << i << " " << j;
                  constraint16[i][j] = IloRange(env, -IloInfinity, w[i] - w[j] - (D - 1) * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint16[i][j]);
                }
            }
        }

      // Constraints 17) w_j +- w_i < = (D - 1)  * (1 - x_{ji}), for all i in S, j in V, i != j
      IloArray<IloRangeArray> constraint17(env, V);
      for(auto i = 0; i < V; i++)
        {
          constraint17[i] = IloRangeArray(env, V);
          for(auto j = 0; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint17_" << i << " " << j;
                  constraint17[i][j] = IloRange(env, -IloInfinity, w[j] - w[i] - (D - 1) * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint17[i][j]);
                }
            }
        }

      // Constraints 18) x_{ij} = 0, for all i in D, j in D
      IloArray<IloRangeArray> constraint18(env, D);
      for(auto i = 0; i < D; i++)
        {
          constraint18[i] = IloRangeArray(env, D);
          for(auto j = 0; j < D; j++)
            {
              if(i != j)
                {
                  name << "constraint18_" << i << " " << j;
                  constraint18[i][j] = IloRange(env, 0, x[i][j], 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint18[i][j]);
                }
            }

        }
      // Constraints 19) L_i + P_i <= 1, for all i in N
      IloRangeArray constraint19(env, N);
      for(auto i = 0; i < N; i++)
        {
          name << "constraint19_" << i;
          int Li = L[i];
      	  constraint19[i] = IloRange(env, -IloInfinity, Li + P[i], 1, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint19);

      for(auto i = 0; i < V; i++)
        {
          model2.add(x[i][i] == 0);
        }

      if(initSol)
        {
          IloNumVarArray startVar(env);
          IloNumArray startVal(env);
          IloNumArray2 xSol(env, V);
          for(auto i = 0; i < V; ++i)
            {
              xSol[i] = IloNumArray(env, V);
              // for(auto j = 0; j < V; j++)
              //   xSol[i][j] = 0;
            }
          // Read initial solution
          int u = 0, v = 0;
          for(auto d = 0; d < D * 2; d++)
            {
              if(s.Ss[d].size())
                {
                  for(unsigned i = 0; i < s.Ss[d].size() - 1; i++)
                    {
                      u = s.Ss[d][i];
                      v = s.Ss[d][i + 1];
                      xSol[u][v] = 1.0;
                      //cout << u << " " << v << endl;
                    }
                }
            }

          for(auto u = 0; u < V; u++)
            {
              for(auto v = 0; v < V; v++)
                {
                  startVar.add(x[u][v]);
                  startVal.add(xSol[u][v]);
                }
              xSol[u].end(); // free
            }
          xSol.end(); // free

          cplex2.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto, "MIPStart");
        }


      // Cplex Parameters
      //timeLimit = 10;
      // cplex2.exportModel("model2.lp");
      cplex2.setParam(IloCplex::Param::Threads, 1);
      cplex2.setParam(IloCplex::Param::TimeLimit, timeLimit);
      cplex2.setParam(IloCplex::Param::Emphasis::MIP, 1);
      cplex2.setParam(IloCplex::Param::MIP::Limits::Solutions, 3);
      cplex2.use(timeLimitCallback(env, cplex2, IloFalse, cplex2.getCplexTime(), timeLimit, currentUB));
      cplex2.setOut(env.getNullStream());
      cplex2.setParam(IloCplex::Param::MIP::Display, 0);

      // Solve
      if(!cplex2.solve())
        {
          cout << "can't solve'" << endl;
          return 99999999;
        }

      // float routingCost = cplex2.getObjValue();
      // float routingTime = cplex2.getTime();
      // cout << "UB   = " << routingCost << endl;
      // cout << "LB " << cplex2.getBestObjValue() << endl;
      // cout << "time = " << routingTime << endl;

      if(flag)
        {
          // reset routes
          //s.noInSs.clear();
          cout << "check the change of s.noInSs" << endl;
          exit(0);

          vector<bool> check(V, true);

          for(auto c = 0; c < D * 2; c++)
            s.Ss[c].clear();
          int route = 0;
          for(unsigned i = 0; i < y.size(); i++)
            s.Sd[i] = y[i];

          for(auto i = 0; i < D; i++)
            {
              for(auto j = D; j < V; j++)
                {
                  if(cplex2.getValue(x[i][j]) > 0.9)
                    {
                      route = i * 2;
                      int u = i, v = j;
                      if(s.Ss[route].size())
                        route++;
                      while(v < V)
                        {
                          if(cplex2.getValue(x[u][v]) > 0.9)
                            {
                              //cout << "(" << u << ", " << v << ") ";
                              s.Ss[route].push_back(u);
                              check[u] = 0;
                              if(v < D)
                                {
                                  //cout << "(" << u << ", " << v << ") " << endl;
                                  //cout << "*" << endl;

                                  s.Ss[route].push_back(v);
                                  //i = route;

                                  // route++;
                                  break;
                                }
                              u = v;
                              v = -1;
                            }
                          v++;
                        }
                    }
                }
            }
          // for(auto i = D; i < V; i++)
          //   {
          //     if(check[i])
          //       s.noInSs.push_back(i);
          //   }

        }

      // getchar();
      // cout << "P: ";
      // for(auto i = 0; i < N; i++)
      //   cout << cplex2.getValue(P[i]) << " ";
      // cout << endl;

      // cout << "L: ";
      // for(auto i = 0; i < N; i++)
      //   cout << L[i] << " ";
      // cout << endl;

      // cout << "x: " << endl;
      float routeCost = 0.0;
      for(auto i = 0; i < V; i++)
        {
          for(auto j = 0; j < V; j++)
            {
              if(cplex2.getValue(x[i][j]) > 0.9)
                {
                  routeCost += t[i][j];
                  // cout << "x[" << i << "][" << j << "] = "<< cplex2.getValue(x[i][j]) << " -> " << t[i][j] << endl;
                }
            }
        }

      // cout << "v: ";
      // for(auto i = 0; i < S; i++)
      //   cout << cplex2.getValue(v[i]) << " ";
      // cout << endl;

      float B2 = 0.0;
      float finalCost = cplex2.getObjValue();
      for(auto i = 0; i < D; i++)
        B2 += C[i][y[i]];


      // cout << "w: ";
      // for(auto i = 0; i < V; i++)
      //   cout << cplex2.getValue(w[i]) << " ";
      // cout << endl;

      // cout << "Depot: ";
      // for(auto i = 0; i < D; i++)
      //   cout << y[i] << " ";
      // cout << endl;

      // cout << "Budget = " << B2 << endl;
      // cout << "Route cost  = " << routeCost << endl;
      // cout << "Total cost" << B2 + routeCost << " " << finalCost <<endl;

      // if(finalCost < M)
      //   {
      //     // cout << "M = " << M << endl;
      //     M = finalCost;
      //     // cout << "updated M = " << M << endl;
      //   }

      return finalCost;
      // cout << "Total time  = " << locationTime + routingTime << endl;
    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
  return -1;
}

void ILS::runMH()
{
  int iterMax = 10;
  float timelimit = 20;
  float currentCost = 0;
  Solution s;
  vector<int> initial(D, 0);
  // depot
  s.Sd = new int[D];

  // areas to cover
  s.coveredArea = new int [N];
  for(auto i = 0; i < N; i++)
    s.coveredArea[i] = 0;

  // initialize routes
  for(auto d = 0; d < D * 2; d++)
    s.Ss.push_back(vector<int>());

  float min = solver(initial, 9999999, 20, s, true, false);

  auto start = chrono::steady_clock::now();
  for(auto iter = 0; iter < iterMax; iter++)
    {
      vector<int> Sd(D, -1);
      int d = 0;
      // test[i] = rnd(0, 1);
      // assures the minimum K = 0
      for(auto i = 0; i < mMin[0]; i++)
        {
          do {
            d = rnd(0, D - 1);
          } while(Sd[d] != -1);
          Sd[d] = 0;
        }
      // assures the minimum K = 1
      for(auto i = 0; i < mMin[1]; i++)
        {
          do {
            d = rnd(0, D - 1);
          } while(Sd[d] != -1);
          Sd[d] = 1;
        }
      for(auto i = 0; i < D; i++)
        if(Sd[i] == -1)
          Sd[i] = rnd(-1, K - 1);

      currentCost = solver(Sd, currentCost, timelimit, s, false, false);
      if(currentCost < min)
        min = currentCost;
    }

  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  auto totalTime  = chrono::duration <double, std::ratio<1>> (diff).count();
  cout << "Cost = " << min << endl;
  cout << "Time = " << totalTime << endl;
}

void ILS::updateElite(Solution &s1, float cost)
{
  if(s1.feasible)
    {
      if(elite.size() > 0)
        {
          map<float, Solution>::iterator rit = --elite.end();
          //cout << "dasddasdasdadgjhdgjkasdgasdhgasdgasjdgaskdjghasdhasgdjkg " <<  rit->first << " " << cost << endl;
          if(elite.size() < sizeElite || rit->first > cost)
            {
              if(elite.find(cost) == elite.end())
                {
                  // cout << ":O" << endl;
                  Solution s;
                  s.Sd = new int[D];

                  // areas to cover
                  s.coveredArea = new int[N];
                  for(auto i = 0; i < N; i++)
                    {
                      s.coveredArea[i] = 0;
                    }

                  // initialize routes
                  for(auto d = 0; d < D * 2; d++)
                    {
                      s.Ss.push_back(vector<int>());
                    }

                  // initialize noInSs
                  s.noInSs = new bool[V];
                  //computeCost(s1);
                  copyS(s1, s);
                  elite.insert(make_pair(cost, s));
                  // cout << population.size() << endl;
                  while(elite.size() > sizeElite)
                    elite.erase(--elite.end());
                }
            }
        }
      else
        {
          // cout << "|:O " << elite.size() << endl;
          Solution s;
          s.Sd = new int[D];

          // areas to cover
          s.coveredArea = new int[N];
          for(auto i = 0; i < N; i++)
            {
              s.coveredArea[i] = 0;
            }

          // initialize routes
          for(auto d = 0; d < D * 2; d++)
            {
              s.Ss.push_back(vector<int>());
            }

          // initialize noInSs
          s.noInSs = new bool[V];
          //computeCost(s1);
          copyS(s1, s);
          //computeCost(s1);
          // print_sol(s);
          elite.insert(make_pair(cost, s));
        }
    }
}

void ILS::run(float *result)
{
  cout << "ILS" << endl;
  auto start = chrono::steady_clock::now();

  Solution s;
  Solution sBest;
  Solution sBestBest;

  // depot
  s.Sd = new int[D];
  sBest.Sd = new int[D];
  sBestBest.Sd = new int[D];

  // areas to cover
  s.coveredArea = new int[N];
  sBest.coveredArea = new int[N];
  sBestBest.coveredArea = new int[N];
  for(auto i = 0; i < N; i++)
    {
      s.coveredArea[i] = 0;
      sBest.coveredArea[i] = 0;
      sBestBest.coveredArea[i] = 0;
    }

  // initialize routes
  for(auto d = 0; d < D * 2; d++)
    {
      s.Ss.push_back(vector<int>());
      sBest.Ss.push_back(vector<int>());
      sBestBest.Ss.push_back(vector<int>());
    }

  // initialize noInSs
  s.noInSs = new bool[V];
  sBest.noInSs = new bool[V];
  sBestBest.noInSs = new bool[V];

  // build initial solution
  int nTimes = 10000;
  for(auto times = 0; times < nTimes; times++)
    {
      cout << times << endl;
      initialSolution(s);
      // cout << "here0?" << endl;
      computeCost(s);
      if(!s.feasible)
        {
          float alpha = randValue(generator);
          if(alpha < 0.25)
            {
              //cout << "pert1" << endl;
              perturbation1(s);
            }
          else if(alpha < 0.50)
            {
              //cout << "pert2" << endl;
              perturbation2(s);
            }
          else if(alpha < 0.75)
            {
              //cout << "pert3" << endl;
              perturbation3(s);
            }
          else
            {
              //cout << "pert4" << endl;
              perturbation4(s);
            }
          computeCost(s);

          for(auto i = 0; i < 10; i++)
            {
              // cout << "2opt" << endl;
              for(auto c = 0; c < D * 2; c++)
                if(s.Ss[c].size())
                  twoOpt(s.Ss[c]);
              // computeCost(s);
              // cout << "remove" << endl;
              removeNode(s);
              // computeCost(s);

              // cout << "swap" << endl;
              LS_swap(s);
              computeCost(s);
            }
          if(s.feasible)
            break;
        }

      for(auto i = 0; i < 10; i++)
        {
          for(auto c = 0; c < D * 2; c++)
            if(s.Ss[c].size())
              twoOpt(s.Ss[c]);
          removeNode(s);
          LS_swap(s);
          computeCost(s);
        }

      if(s.feasible)
        break;
      // if(times > 9 && s.feasible == 0)
      //   {
      //     vector<int> initial(D);
      //     //cout << "*********************************" << endl;
      //     for(auto i = 0; i < D; i++)
      //       {
      //         //cout << s.Sd[i] << " ";
      //         initial[i] = s.Sd[i];
      //       }
      //     // cout << endl;
      //     // cout << "*********************************" << endl;
      //     solver(initial, 999999, 10, s, true, false);
      //     computeCost(s);
      //     // print_sol(s);
      //     // exit(0);
      //     // getchar();
      //     if(s.feasible)
      //       break;
      //   }
    }
  // cout << "dsliii" << endl;

  // exit(0);

  for(auto t = 0; t < 3; t++)
    {
      for(auto c = 0; c < D * 2; c++)
        if(s.Ss[c].size())
          twoOpt(s.Ss[c]);
      removeNode(s);
      LS_swap(s);
    }

  float s_cost = computeCost(s);
  cout << "cost " << s_cost << " " << s.feasible << endl;
  // print_sol(s);

  vector<float> solutions;
  float s_bestCost = 99999999;
  float s_bestBestCost = 99999999;
  if(s.feasible)
    {
      s_bestBestCost = s_bestCost = s_cost;
      copyS(s, sBest);
      copyS(sBest, sBestBest);
      solutions.push_back(s_cost);
    }

  float alpha = 0;
  for(auto iter = 1; iter < 1000000; iter++)
    {
      // cout << "perturbation1 1" << endl;
      if(iter % 10000 == 0)
        cout << "iter " << iter << endl;

      if(iter % 20000 == 0)
        {
          // cout << "reset" << endl;
          s_bestCost = 9999999;

          if(elite.size())
            {
              map<float, Solution>::iterator it = elite.begin();
              int position = rnd(0, elite.size() - 1);
              // cout << "rnd " << position << endl;
              advance(it, position);
              copyS(it->second, s);
              s_cost = it->first;
            }
          else
            {
              initialSolution(s);
              s_cost = computeCost(s);
            }
          // 2-opt for each route
          for(auto c = 0; c < D * 2; c++)
            {
              if(s.Ss[c].size())
                twoOpt(s.Ss[c]);
          }
          s_cost = computeCost(s);
          removeNode(s);
          s_cost = computeCost(s);
          LS_swap(s);
          s_cost = computeCost(s);
        }
      if(iter % 2000 == 0)
        {
          initialSolution(s);
          s_cost = computeCost(s);
          // cout << "here?" << endl;
        }
      float test = randValue(generator);
      if(test < 0.3)
        {
          alpha = randValue(generator);
          // cout << "alpha " << alpha << endl;

          if(alpha < 0.25)
            {
              // cout << "pert1" << endl;
              perturbation1(s);
            }
          else if(alpha < 0.50)
            {
              // cout << "pert2" << endl;
              perturbation2(s);
            }

          else if(alpha < 0.75)
            {
              // cout << "pert3" << endl;
              perturbation3(s);
            }
          else
            perturbation4(s);
          // cout << "here3" << endl;

          computeCost(s);
        }
      else if(test < 0.6 && elite.size())
        {

          // cout << "qua? " << elite.size() << endl;
          // print_sol(s);
          // cout << "dsadasdadasd" << endl;

          map<float, Solution>::iterator it = elite.begin();
          int position = rnd(0, elite.size() - 1);

          // cout << "rnd " << position << endl;

          // if(position != 0)
          advance(it, position);
          copyS(it->second, s);
          // print_sol(s);
          s_cost = computeCost(s);
          //s_cost = it->first;
        }

      // cout << "2-opt" << endl;
      // 2-opt for each route
      for(auto c = 0; c < D * 2; c++)
        {
          //cout << "2-opt " << c << endl;
          if(s.Ss[c].size())
            twoOpt(s.Ss[c]);
        }

      // for(auto i = 0; i < N; i++)
      //   cout << s.coveredArea[i] << " ";
      // cout << endl;
      s_cost = computeCost(s);
      // cout << "cost " << s_cost << endl;
      // print_sol(s);

      // cout << "removing nodes " << endl;
      removeNode(s);
      s_cost = computeCost(s);
      // print_sol(s);

      // LS: swap
      // cout << "swap" << endl;
      LS_swap(s);
      s_cost = computeCost(s);
      //cout << "iter " << iter << ": "<< s_cost << endl;
      if((s_bestCost > s_cost) && s.feasible)
        {
          s_bestCost = s_cost;
          copyS(s, sBest);
          // cout << "qui" << endl;
          if(s_bestBestCost > s_cost)
            {
              s_bestBestCost = s_cost;
              copyS(sBest, sBestBest);
              cout << iter << " "<< s_cost << endl;
              solutions.push_back(s_cost);
            }

          //getchar();
        }
      // cout << "iteration ended" << endl;
      // for(auto i = 0; i < N; i++)
      //   cout << s.coveredArea[i] << " ";
      // cout << endl;
      //exit(0);

    }
  cout << "Best solution:" << endl;
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  auto totalTime  = chrono::duration <double, std::ratio<1>> (diff).count();

  s_bestBestCost = computeCost(sBestBest);
  cout << "Best cost " << s_bestBestCost << endl;
  print_sol(sBestBest);
  cout << "Best cost " << s_bestBestCost << endl;
  cout << "Best solution found:" << endl;

  if(solutions.size())
    {
      for(unsigned i = 0; i < solutions.size(); i++)
        cout << i+1 << " "<< solutions[i] << endl;
      cout << "Time = " << totalTime << endl;
    }
  else
    {
      cout << "A solution feasible are not found" << endl;
      return;
    }

  // cout << "Sd: ";
  // for(auto i = 0; i < D; i++)
  //   cout << s.Sd[i] << " ";
  cout << "Elite solution:" << endl;
  for(std::map<float, Solution>::iterator it = elite.begin(); it != elite.end(); ++it)
    cout << it->first << endl;


  result[0] = solutions[0];
  result[1] = s_bestBestCost;
  result[2] = totalTime;

  delete s.Sd;
  delete sBest.Sd;
  delete sBestBest.Sd;
  delete s.coveredArea;
  delete sBest.coveredArea;
  delete sBestBest.coveredArea;
  delete s.noInSs;
  delete sBest.noInSs;
  delete sBestBest.noInSs;
}

void ILS::coveringAreaByDepot(Solution &s, int d, int k)
{
  int cover = -1;
  // cout << "considerado " << d << " cubre a:" << endl;
  for(unsigned i = 0; i < alpha[d][k].size(); i++)
    {
      cover = alpha[d][k][i];
      // cout << cover << " ";
      s.coveredArea[cover]++;
    }

  // cout << endl;
}

void ILS::coveringAreaByNode(Solution &s, int i)
{
  int cover = -1;
  // cout << "considerado " << i << " cubre a:" << endl;
  for(unsigned k = 0; k < beta[i].size(); k++)
    {
      cover = beta[i][k];
      // cout << cover << " ";
      s.coveredArea[cover]++;
    }
  // cout << endl;
}

void ILS::perturbation1(Solution &s)
{
  // cout << "dhajksdashkajdhaskjdhsakdjhaskjdhasdkjhasjkdhasjkdhaskdhs" << endl;

  //print_sol(s);
  vector<int> candidate;

  for(auto i = D; i < V; i++)
    {
      if(!s.noInSs[i])
        candidate.push_back(i);
    }

  if(candidate.size())
    {
      int route = 0;
      do {
        route = rnd(0, s.Ss.size() - 1);
      } while(s.Ss[route].size() == 0);

      //cout << "route " << route << endl;
      // cout << s.noInSs.size() << endl;

      int i = rnd(0, candidate.size() - 1);

      //cout << "noInSs[i] = " << s.noInSs[i] << endl;

      int r = rnd(1, s.Ss[route].size() - 2);

      //cout << "r = " << r << endl;

      s.Ss[route].insert(s.Ss[route].begin() + r, candidate[i]);
      s.noInSs[candidate[i]] = true;
    }
  // print_sol(s);
  // cout << "dhajksdashkajdhaskjdhsakdjhaskjdhasdkjhasjkdhasjkdhaskdhs" << endl;

}

void ILS::perturbation2(Solution &s)
{
  // cout << "perturbation2 " << computeCost(s) << endl;
  int route = 0;
  do {
    route = rnd(0, s.Ss.size() - 1);
  } while(s.Ss[route].size() == 0);

  if(s.Ss[route].size() > 3)
    {
      //cout << "route " << route << endl;
      // cout << s.noInSs.size() << endl;

      int i = rnd(1, s.Ss[route].size() - 2);
      //cout << "Ss[i] = " << s.Ss[route][i] << endl;

      int j = rnd(1, s.Ss[route].size() - 2);
      //cout << "Ss[j] = " << s.Ss[route][j] << endl;

      int aux = s.Ss[route][i];
      s.Ss[route][i] = s.Ss[route][j];
      s.Ss[route][j] = aux;
      //cout << "done" << endl;
    }
}

void ILS::perturbation3(Solution &s)
{
  // cout << "perturbation3 " << computeCost(s) << endl;
  int count = 0;
  vector<int> r;

  for(auto i = 0; i < D * 2 && count < 2; i++)
    {
      if(s.Ss[i].size() > 2)
        {
          r.push_back(i);
          count++;
        }
    }

  if(count > 1)
    {
      int numRnd = rnd(0, r.size() - 1);
      int route1 = r[numRnd];

      int route2 = 0;
      do {
        numRnd = rnd(0, r.size() - 1);
        route2 = r[numRnd];
      } while(route1 == route2);
      // cout << "route 1 " << route1 << endl;
      // cout << "route 2 " << route2 << endl;
      int i = rnd(1, s.Ss[route1].size() - 2);
      // cout << "Ss[r1][i] = " << s.Ss[route1][i] << endl;

      int j = rnd(1, s.Ss[route2].size() - 2);
      // cout << "Ss[r2][j] = " << s.Ss[route2][j] << endl;

      int aux = s.Ss[route1][i];
      s.Ss[route1][i] = s.Ss[route2][j];
      s.Ss[route2][j] = aux;
      //cout << "done" << endl;
    }
}

void ILS::perturbation4(Solution &s)
{
  vector<int> index;
  int d = 0;
  // cout << "perturbation4 " << computeCost(s) << endl;
  // print_sol(s);
  int minK[2] = {0, 0};
  for(unsigned i = 0; i < s.Ss.size(); i += 2)
    {
      d = i / 2;
      if(s.Sd[d] == 0 && s.Ss[i].size())
        {
          minK[0]++;
          if(s.Ss[i].size() < 5 && s.Ss[i + 1].size() < 5)
            index.push_back(i);
        }
    }
  if(!index.size() || minK[0] < 2)
    return;

  // for(unsigned i = 0; i < index.size(); i++)
  //   cout << "route " << index[i]<< endl;

  int j = rnd(0, index.size() - 1);
  int route = index[j];
  // cout << "chosen route " << route << endl;

  d = route / 2;

  // cout << "lalal" << s.Ss[route].size() << endl;

  vector<int> nodes;
  for(unsigned i = 1; i < s.Ss[route].size() - 1; i++)
    nodes.push_back(s.Ss[route][i]);

  // cout << "var" << endl;

  s.Ss[route].clear();
  if(s.Ss[route + 1].size())
    {
      for(unsigned i = 1; i < s.Ss[route + 1].size() - 1; i++)
        nodes.push_back(s.Ss[route + 1][i]);
      s.Ss[route + 1].clear();
    }

  for(unsigned i = 0; i < alpha[d][0].size(); i++)
    {
      s.coveredArea[alpha[d][0][i]]--;
      // nodes.push_back(alpha[d][0][i]);
      // for(unsigned j = 0; j < s.noInSs.size(); j++)
      //   {
      //     if(nodes.back() == s.noInSs[j])
      //     {
      //       s.noInSs.erase(s.noInSs.begin() + j);
      //       break;
      //     }
      //   }
    }

  shuffle(nodes.begin(), nodes.end(), default_random_engine(seed));
  s.Sd[d] = -1;
  int w = 0, u = 0, v = 0;
  for(unsigned i = 0; i < nodes.size(); i++)
    {
      w = nodes[i];
      //cout << "node " << w  << endl;

      float minimum = 9999999.9;
      int here = -1;
      for(auto c = 0; c < D * 2; c++)
        {
          d = c / 2;
          // cout << "r" << c << endl;
          if(s.Sd[d] == 0)
            {
              if(s.Ss[c].size() > 0)
                {
                  for(unsigned r = 0; r < s.Ss[c].size() - 1; r++)
                    {
                      u = s.Ss[c][r];
                      v = s.Ss[c][r + 1];
                      // cout << u << " " << v << "-> ";
                      float newCost = t[u][w] + t[w][v] - t[u][v];
                      // cout << newCost << " " << endl;
                      if(newCost < minimum)
                        {
                          //cout << u << " " << v << "-> " << newCost << " " << endl;
                          newCost = minimum;
                          here = r;
                          route = c;
                        }

                    }
                }
            }
        }
      if(here != -1)
        {
          //cout << "agrego " << w << endl;
          s.Ss[route].insert(s.Ss[route].begin() + here + 1, w);
          coveringAreaByNode(s, w);
        }
    }
  // computeCost(s);
  // print_sol(s);

  // cout << "cA = ";

  // for(auto i = 0; i < N; i++)
  //   {
  //     cout << s.coveredArea[i] << " ";
  //   }
  // cout << endl;

  vector<int> missing;
  vector<bool> mark(V, false);
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size())
        {
          for(unsigned j = 0; j < s.Ss[i].size() - 1; j++)
            mark[s.Ss[i][j]] = true;
        }
    }

  // Nodes subset according to uncovered areas
  for(auto i = 0; i < N; i++)
    {
      if(!s.coveredArea[i])
        {
          for(unsigned j = 0; j < bi[i].size(); j++)
            {
              int node = bi[i][j];
              if(!mark[node])
                {
                  mark[node] = true;
                  missing.push_back(node);
                }
            }
        }
    }
  // cout << "missing: " << endl;

  // for(unsigned i = 0; i < missing.size(); i++)
  //   {
  //     cout << missing[i] << " ";

  //   }
  // cout << endl;

  if(missing.size())
    {
      //cout << "quedan (" << missing.size() << ")"<< endl;
      // check if the solution is feasible
      for(auto ii = 0; ii < N; ii++)
        {
          if(!s.coveredArea[ii])
            {
              for(auto i = missing.size() - 1; missing.size(); i--)
                {
                  int u, v, w;
                  float minimum = 9999999.9;
                  int here = -1, route = 0;
                  w = missing.back();
                  //cout << "arreglar " << w << endl;
                  for(auto c = 0; c < D * 2; c++)
                    {
                      d = c / 2;
                      // cout << "r" << c << endl;
                      if(s.Sd[d] == 0)
                        {
                          if(s.Ss[c].size() > 0)
                            {
                              for(unsigned r = 0; r < s.Ss[c].size() - 1; r++)
                                {
                                  u = s.Ss[c][r];
                                  v = s.Ss[c][r + 1];

                                  // cout << u << " " << v << "-> ";
                                  float newCost = t[u][w] + t[w][v] - t[u][v];
                                  // cout << newCost << " " << endl;
                                  if(newCost < minimum)
                                    {
                                      //cout << u << " " << v << "-> " << newCost << " " << endl;
                                      newCost = minimum;
                                      here = r;
                                      route = c;
                                    }

                                }
                            }
                          // else
                          //   {
                          //     //cout << "no había y se agrega: " << d << " " << w << " " << d<< endl;
                          //     s.Ss[c].push_back(d);
                          //     s.Ss[c].push_back(w);
                          //     coveringAreaByNode(s, w);
                          //     s.Ss[c].push_back(d);
                          //     missing.pop_back();
                          //     w = missing.back();
                          //     here = -1;
                          //     break;
                          //   }
                        }
                    }
                  if(here != -1)
                    {
                      //cout << "agrego " << w << endl;

                      s.Ss[route].insert(s.Ss[route].begin() + here + 1, w);
                      coveringAreaByNode(s, w);
                      //s.noInSs.erase(remove(s.noInSs.begin(), s.noInSs.end(), w), s.noInSs.end());
                      s.noInSs[w] = true;
                      missing.pop_back();
                      ii = 0;
                      break;
                    }
                }
              // test = false;
            }
        }
    }

  // computeCost(s);
  // print_sol(s);

  //exit(0);

}
void ILS::computeRoute(Solution &s, vector<float> &c)
{
  c = vector<float>(D * 2, 0);
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size())
        {
          for(unsigned j = 1; j < s.Ss[i].size(); j++)
            c[i] += t[s.Ss[i][j - 1]][s.Ss[i][j]];
        }
    }
}

bool ILS::LS_swap(Solution &s)
{
  vector<float> Cost;
  computeRoute(s, Cost);
  // print_sol(s);
  bool flag = false;
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      //cout << "que ondi " << i << " really?! " <<  s.Ss[i].size() << endl;
      if(s.Ss[i].size())
        {
          for(unsigned j = 1; j < s.Ss[i].size() - 1; j++)
            {
              for(unsigned route = 0; route < s.Ss.size(); route++)
                {
                  if(s.Ss[route].size())
                    {
                      for(unsigned l = 1; l < s.Ss[route].size() - 1; l++)
                        {
                          if(i == route && j == l)
                            continue;
                          else
                            {
                              //cout << i << " "<< j << " - " << route << " "<< l << endl;
                              // d(s[i - 1], s[l]) + d(s[l], s[i + 1]) - d(s[i - 1], i) - d(s[i], s[i + 1])
                              int opt = 0;
                              int a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;
                              if(s.Ss[route][l] != s.Ss[i][j + 1])
                                {
                                  //cout << "norep" << endl;
                                  e = t[s.Ss[i][j - 1]][s.Ss[route][l]];
                                  f = t[s.Ss[route][l]][s.Ss[i][j + 1]];

                                  a = t[s.Ss[i][j - 1]][s.Ss[i][j]];
                                  b = t[s.Ss[i][j]][s.Ss[i][j + 1]];

                                  g = t[s.Ss[route][l - 1]][s.Ss[i][j]];
                                  h = t[s.Ss[i][j]][s.Ss[route][l + 1]];

                                  c = t[s.Ss[route][l - 1]][s.Ss[route][l]];
                                  d = t[s.Ss[route][l]][s.Ss[route][l + 1]];
                                  if(!s.feasible)
                                    {
                                      //cout << "dsadasd "  << e - a +  h - d << endl;

                                      if(i == route && (e - a +  h - d < 0))
                                        opt = 1;
                                      if(i != route && (e + f - a - b + Cost[i] < T) && (g + h - c - d + Cost[route] < T))
                                        opt = 1;
                                    }
                                  else
                                    {
                                      //cout << "lalalala "  << e + f - a - b + g + h - c - d << endl;

                                      if(s.Ss[route].size() > 3)
                                        {
                                          if((e + f - a - b + g + h - c - d < 0) && (e + f - a - b + Cost[i] < T) && (g + h - c - d + Cost[route] < T))
                                            opt = 1;
                                        }
                                      else
                                        {
                                          if((e + f - a - b + g + h - c - d < 0) && (e + f - a - b + Cost[i] < T) && (g + h - c - d + Cost[route] < T))
                                            opt = 1;

                                        }
                                    }
                                }
                              else
                                {
                                  //cout << "rep" << endl;
                                  a = t[s.Ss[i][j - 1]][s.Ss[i][j]];
                                  b = t[s.Ss[route][l]][s.Ss[route][l + 1]];
                                  c = t[s.Ss[i][j - 1]][s.Ss[route][l]];
                                  d = t[s.Ss[i][j]][s.Ss[route][l + 1]];

                                  if(- a - b + c + d < 0)
                                    opt = 1;
                                }

                              if(opt)
                                {
                                  //cout << "cambio " << s.Ss[i][j] << " <-> " << s.Ss[route][l] << endl;
                                  int aux = s.Ss[i][j];
                                  s.Ss[i][j] = s.Ss[route][l];
                                  s.Ss[route][l] = aux;
                                  computeCost(s);
                                  // print_sol(s);
                                  computeRoute(s, Cost);
                                  //exit(0);

                                  // return true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  // print_sol(s);
  // cout << "end swap" << endl;

  return flag;
}

bool ILS::removeNode(Solution &s)
{
  if(!s.feasible)
    return false;
  int node = 0, area = 0;
  int route = 0;

  // for(auto ii = 0; ii < N; ii++)
  //   cout << s.coveredArea[ii] << " ";
  // cout << endl;
  // print_sol(s);
  do {
    route = rnd(0, s.Ss.size() - 1);
  } while(s.Ss[route].size() == 0);
  // for(auto route = 0; route < s.Ss.size(); route++)
  //   if(s.Ss[route].size())
  {

    // cout << "recomputelalala" << endl;

    // for(auto ii = 0; ii < N; ii++)
    //   cout << s.coveredArea[ii] << " ";
    // cout << endl;
    for(unsigned i = 1; i < s.Ss[route].size() - 1; i++)
      {
        if(s.Ss[route].size() > 3)
          {
            node = s.Ss[route][i];
            bool flag = true;
            //cout << node << endl;
            for(unsigned j = 0; j < beta[node].size(); j++)
              {
                area = beta[node][j];
                // cout << area << " ";
                if(s.coveredArea[area] < 2)
                  {
                    flag = false;
                    break;
                  }
              }
            // cout << endl;
            if(flag)
              {
                // cout << "se va " << endl;
                for(unsigned j = 0; j < beta[node].size(); j++) // decreasing covered areas
                  {
                    // cout << beta[node][j] << " " << endl;;
                    s.coveredArea[beta[node][j]] -= 1;

                    // for(auto ii = 0; ii < N; ii++)
                    //   cout << s.coveredArea[ii] << " ";
                    // cout << endl;
                  }
                // cout << "que " << i << endl;
                s.Ss[route].erase(s.Ss[route].begin() + i); // remove node
                s.noInSs[node] = false; // add to list noInSs
                i--;
              }
            // cout << endl;

            // for(auto ii = 0; ii < N; ii++)
            //   cout << s.coveredArea[ii] << " ";
            // cout << endl;
            // computeCost(s);
            // cout << "recompute" << endl;

            // for(auto ii = 0; ii < N; ii++)
            //   cout << s.coveredArea[ii] << " ";
            // cout << endl;
          }
      }
  }
  //print_sol(s);
  // if all the nodes are removed so remove this route
  // if(s.Ss[route].size() == 2)
  //   {
  //     s.Ss[route].pop_back();
  //     s.Ss[route].pop_back();
  //     if(route % 2 == 0 && s.Ss[route + 1].size() == 0)
  //       {
  //         s.Sd[(route / 2)] = 0;
  //         cout << "lo hice" << endl;
  //         exit(0);
  //       }
  //     else if(route % 2 == 1 && s.Ss[route - 1].size() == 0)
  //       s.Sd[(route / 2)] = 0;
  //   }

  // for(unsigned i = 0; i < s.Ss.size(); i++)
  //   {
  //     if(i % 2 == 0 && s.Sd[i / 2] == 1 && s.Ss[i].size() == 0 && s.Ss[i + 1].size() == 0)
  //       {
  //         cout << "chau" << endl;

  //         s.Sd[(route / 2)] = 0;
  //         exit(0);

  //       }
  //   }


  // for(auto i = 0; i < N; i++)
  //   cout << s.coveredArea[i] << " ";
  // cout << endl;

  for(auto i = 0; i < N; i++)
    if(!s.coveredArea[i])
      {
        s.feasible = 0;
        cout << "infeasible?" << endl;
        return false;
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
      if(s.Sd[i] > -1)
        {
          s.Lcost += C[i][s.Sd[i]];
          coveringAreaByDepot(s, i, s.Sd[i]);
        }
    }

  // cout << "\nLcost: "<< s.Lcost << endl;
  for(auto i = D; i < V; i++)
    s.noInSs[i] = false;

  // compute drones' routes
  int u = 0, v = 0;
  float routeCost = 0.0;
  float routeCostReal = 0.0;
  s.Rcost = 0.0;
  s.feasible = 1;
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size() > 0 && s.Sd[i / 2] != 0)
        {
          cout << "error this route is not watchtower" << endl;
          print_sol(s);
          exit(0);
        }

      if(s.Ss[i].size())
        {
          routeCost = 0;
          routeCostReal = C[0][3];

          for(unsigned j = 0; j < s.Ss[i].size() - 2; j++)
            {
              u = s.Ss[i][j];
              v = s.Ss[i][j + 1];
              // cout << u << " " << v << " " << t[u][v]<< endl;
              routeCost += t[u][v];
              routeCostReal +=  0.02 * t[u][v];
              // cout << "covering node " << v << endl;
              coveringAreaByNode(s, v);
              s.noInSs[u] = 1;
              s.noInSs[v] = 1;
            }
          u = s.Ss[i][s.Ss[i].size() - 2];
          v = s.Ss[i][s.Ss[i].size() - 1];
          s.noInSs[u] = 1;
          s.noInSs[v] = 1;

          // Remove it after!
          if(v != s.Ss[i][0])
            {
              print_sol(s);
              cout << "error!!!!!!!!!!!!!!!!!!!!!" << endl;
              exit(0);
            }

          // cout << u << " " << v << endl;
          routeCost += t[u][v];
          routeCostReal +=  0.02 * t[u][v];
          if(routeCost > T)
            s.feasible = 0;
          s.Rcost += routeCostReal;
        }
    }


  // for(auto i = D; i < V; i++)
  //   {
  //     if(check[i] != s.noInSs[i])
  //       {
  //         cout << "it is repetead " << i << endl;
  //         print_sol(s);
  //         exit(0);
  //       }
  //   }

  // cout << "Rcost: "<< s.Rcost << endl;
  // final cost
  s.cost = s.Lcost + s.Rcost;
  // check if the solution is feasible
  for(auto i = 0; i < N; i++)
    if(!s.coveredArea[i])
      s.feasible = 0;
  // cout << "Total " << s.cost  << " Feasible "<< s.feasible << endl;

  // for(auto i = 0; i < N; i++)
  //   cout << s.coveredArea[i] << " ";
  // cout << endl;
  if(routeCost > 10000)
    {
      print_sol(s);
      cout << "error big number !!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(0);
    }
  // cout << "oli" << endl;
  //print_sol(s);
  updateElite(s, s.cost);
  return s.cost;
}

int ILS::rnd(unsigned low, unsigned high)
{
  return low + generator() % ((high + 1) - low);
}

void ILS::print_sol(Solution &s)
{
  float cost = 0.0;
  cout << "Sd:" << endl;
  for(auto i = 0; i < D; i++)
    {
      cout << s.Sd[i] << " ";
      cost += C[i][s.Sd[i]];
    }
  cout << endl;
  cout << "Ss:" << endl;

  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(i % 2 == 0)
        cout << "Depot " << i / 2 << ": "<< endl;

      cout << "\tRoute_" << i % 2 << " (" << i << "): ";
      if(s.Ss[i].size())
        {
          float temp = 0.0;
          for(unsigned j = 0; j < s.Ss[i].size() - 1; j++)
            {
              cout << s.Ss[i][j] << " ";
              temp += t[s.Ss[i][j]][s.Ss[i][j + 1]];
            }
          cout << s.Ss[i][s.Ss[i].size() - 1] << " ";
          cout << ": " << temp;
          cost += temp;
        }
      cout << endl;
    }

  cout << "noInSs = ";
  for(auto i = D; i < V; i++)
    if(!s.noInSs[i])
    cout << i << " ";
  cout << endl;

  cout << "cA = ";
  for(auto i = 0; i < N; i++)
    cout << s.coveredArea[i] << " ";
  cout << endl;
  cout << "Feasible = " << s.feasible << endl;
  cout << "cost = " << cost << endl;
}
