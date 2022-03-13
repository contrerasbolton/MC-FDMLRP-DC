/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "ils.hpp"

ILS::ILS(int seed, int N, int D, int We, int W, int S, int B, int **beta, int ND, float T, int V, float M,
         vector<vector<int> > &b, vector<vector<int> > &bi, float **t, int *mMax, int *mMin, float **C, float costUAV, float *parameters)
{
  this->seed = seed;
  this->N = N;
  this->D = D;
  this->We = We;
  this->W = W;
  this->WD = D + We + W;
  this->S = S;
  this->B = B;
  this->ND = ND;
  this->T = T;
  this->V = V;
  this->beta = b;
  this->bi = bi;
  this->t = t;
  this->mMax = mMax;
  this->mMin = mMin;
  this->C = C;
  this->M = M;
  this->b = beta;
  this->costUAV = costUAV;

  // Initialize
  generator.seed(seed);
  randValue = uniform_real_distribution<double>(0.0, 1.0);
  decipresi = 10;
  timeCplex = 10;
  // Parameters
  E  = parameters[0]; // size of elite solutions
  Kt = parameters[1]; // iterations number
  Km = parameters[2]; // number that every Km iterations calls to the MILP solver
  Kg = parameters[3]; // number that every Kg iterations calls to the global reset
  Kl = parameters[4]; // number that every Kl iterations calls to the local reset
  p1 = parameters[5]; // probability for perturbation 1
  p2 = parameters[6]; // probability for perturbation 2
  p3 = parameters[7]; // probability for perturbation 3
  p4 = parameters[8]; // probability for perturbation 4
  p5 = parameters[9]; // probability for perturbation 5
  wp  = parameters[10]; // probability to fix variables to zero in the MILP solver
  Tc = parameters[11]; // the time limit to solve by Cplex

  cout << "Parameters: " << endl;
  cout << "E  = " << E << endl;
  cout << "Kt = " << Kt << endl;
  cout << "Km = " << Km << endl;
  cout << "Kg = " << Kg << endl;
  cout << "Kl = " << Kl << endl;
  cout << "p1 = " << p1 << endl;
  cout << "p2 = " << p2 << endl;
  cout << "p3 = " << p3 << endl;
  cout << "p4 = " << p4 << endl;
  cout << "p5 = " << p5 << endl;
  cout << "w  = " << wp << endl;
  cout << "Tc = " << Tc << endl;
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

void ILS::printOutput(string nameInstance)
{
  map<float, Solution>::iterator it = elite.begin();

  Solution s = it->second;

  float cost = 0.0, bCost = 0.0, wteCost = 0.0, wtCost = 0.0;
  int ballons = 0, watchtowers = 0, sumL = 0, sumP = 0;
  int cover = 0;
  cout << "Sd:" << endl;
  for(auto i = 0; i < WD; i++)
    {
      cout << s.Sd[i] << " ";
      if(i < D && s.Sd[i])
        {
          bCost += C[i][1]; // ballon cost
          cost += C[i][1];
          ballons++;
        }
      else if(i >= D && i < D + We && s.Sd[i])
        {
          wteCost += C[i][2]; // wte cost
          cost += C[i][2];
        }
      else if(s.Sd[i])
        {
          wtCost = C[i][0]; // wt cost
          cost += C[i][0];
          watchtowers++;
        }
    }
  cout << endl;
  cout << "Ss:" << endl;

  vector<bool> coveredArea(N, 0);
  for(auto i = 0; i < WD; i++)
    {
      if(s.Sd[i])
        {
          for(unsigned k = 0; k < beta[i].size(); k++)
            {
              cover = beta[i][k];
              coveredArea[cover] = 1;
            }
        }
    }
  for(auto i = 0; i < N; i++)
    {
      if(coveredArea[i])
        sumL++;
    }

  int nDrones = 0;
  set<int> coverDrone;
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(i % 2 == 0)
        {
          cout << "Depot " << i / 2 << ": "<< endl;
        }

      cout << "\tRoute_" << i % 2 << " (" << i << "): ";
      if(s.Ss[i].size())
        {
          nDrones++;
          float temp = 0.0;
          for(unsigned j = 0; j < s.Ss[i].size() - 1; j++)
            {
              int node = s.Ss[i][j];
              cout << node << " ";
              temp += t[node][s.Ss[i][j + 1]];
              if( j > 0 && j < s.Ss[i].size() - 1)
                {
                  for(unsigned k = 0; k < beta[node].size(); k++)
                    {
                      cover = beta[node][k];
                      if(coveredArea[cover] == 0)
                        coverDrone.insert(cover);
                    }
                }
            }

          cout << s.Ss[i][s.Ss[i].size() - 1] << " ";
          cout << ": " << temp;
          cost += costUAV * temp + C[0][3];
        }
      cout << endl;
      //cout << i << endl;
    }
  sumP = coverDrone.size();
  cout << "noInSs = ";
  for(auto i = WD; i < V; i++)
    if(!s.noInSs[i])
      cout << i << " ";
  cout << endl;

  cout << "cA = ";
  for(auto i = 0; i < N; i++)
    cout << s.coveredArea[i] << " ";
  cout << endl;
  cout << "Feasible = " << s.feasible << endl;
  cout << "cost = " << cost << endl;


  float DroneCost = nDrones * C[0][3];
  int totalCover = sumL + sumP;
  cout << "Approach            = " << nameAlgorithm << endl;
  cout << "Instance            = " << nameInstance << endl;
  cout << "T                   = " << T << endl;
  cout << "Initial Cost        = " << initialCost << endl;
  cout << "Obtained Cost       = " << costFinal << endl;
  cout << "Total Cost w/o RC   = " << it->second.cost - it->second.Rcost << endl;
  cout << "Route Cost          = " << it->second.Rcost << endl;
  cout << "Watchtowers Cost    = " << wtCost << endl;
  cout << "Ballons Cost        = " << bCost << endl;
  cout << "Drones Cost         = " << DroneCost << endl;
  cout << "Drone number        = " << nDrones << endl;
  cout << "existing WT number  = " << We << endl;
  cout << "Watchtowers number  = " << watchtowers << endl;
  cout << "Ballons number      = " << ballons << endl;
  cout << "Covering            = " << totalCover << endl;
  cout << "Covering by Facili. = " << sumL << endl;
  cout << "Covering by Drones  = " << sumP << endl;
  cout << "Time                = " << timeF << endl;

}

bool ILS::initialSolution(Solution &s)
{
  // cout << "Start Initial Solution" << endl;
  // reset covered Area
  for(auto i = 0; i < N; i++)
    s.coveredArea[i] = 0;

  for(auto i = D; i < We; i++)
    s.Sd[i] = 1;

  // reset depot
  for(auto i = D + We; i < WD; i++)
    s.Sd[i] = rnd(0, 1);

  // at least one wt actived
  bool atLeast = true;
  int d = 0;
  for(auto i = D; i < WD; i++)
    {
      if(s.Sd[i])
        {
          atLeast = false;
          break;
        }
    }
  if(atLeast)
    {
      if(We > 0)
        d = D + We;
      else
        d = rnd(D, WD - 1);
      s.Sd[d] = 1;
    }

  // covering areas
  for(auto i = 0; i < WD; i++)
    if(s.Sd[i])
      coveringAreaByNode(s, i);

  // reset routes
  for(auto c = 0; c < WD * 2; c++)
    s.Ss[c].clear();

  // cout << "depots:" << endl;
  // for(auto i = 0; i < WD; i++)
  //   cout << s.Sd[i] << " ";
  // cout << endl;

  // cout << "only depot" << endl;

  // for(auto i = 0; i < N; i++)
  //   cout << s.coveredArea[i] << " ";
  // cout << endl;

  vector<int> nodes;
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
              if(node >= WD && !s.noInSs[node])
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

  // float currentCost = 0.0;
  int u = 0, v = 0;
  // cout << D << " " << WD << endl;

  for(auto c = D * 2; c < WD * 2; )
    {
      d = c / 2;
      vector<int> aux;
      float routeCost = 0.0;
      // cout << "Test in Sd: " << d << ": " << c << " = " << s.Sd[d] << endl;
      if(s.Sd[d] == 1 && nodes.size())
        {
          aux.push_back(d);
          u = nodes.back();
          routeCost += t[d][u];
          // cout << "("<< d << " " << u << ")* ";
          while(nodes.size())
            {
              u = nodes.back();
              if(nodes.size() > 1)
                v = nodes[nodes.size() - 2];
              else
                v = -1;
              //cout << u << " " << v << " " << endl;
              if(v != -1 && routeCost + t[u][v] <= T)
                {
                  routeCost += t[u][v];
                  coveringAreaByNode(s, u);
                  // cout << "("<< u << " " << v << ") ";
                  aux.push_back(u);
                  nodes.pop_back();

                  if(!nodes.size()) // there no nodes
                    {
                      // cout << "what" << endl;
                      routeCost += t[u][d];
                      aux.push_back(d);
                      s.Ss[c] = aux;
                      break;
                    }
                }
              else // returning to depot
                {
                  u = aux.back();
                  // cout << "var " << routeCost + t[u][d] << endl;
                  if(routeCost + t[u][d] <= T) // returning from last u
                    {
                      // cout << "("<< u << " " << d << ")_ " << endl;
                      aux.push_back(d);
                      // cout << "var1" << endl;
                      routeCost += t[u][d];
                      // cout << "var2 "<< c << " " << s.Ss.size() << endl;
                      s.Ss[c] = aux;
                      // cout << "var3" << endl;
                      if(nodes.size() == 1) // there no nodes
                        {
                          // cout << "fin" << endl;
                          //routeCost += t[u][d];
                          // aux.push_back(d);
                          // s.Ss[c] = aux;
                          nodes.pop_back();
                          break;
                        }
                    }
                  else //cannot returning from last u
                    {
                      if(aux.size() > 2) // at least a depot and two node, so remove the last node
                        {
                          // cout << "("<< u << " " << d << ")_ " << endl;
                          coveringAreaByNode(s, u);
                          aux.push_back(d);
                          routeCost += t[u][d];
                          s.Ss[c] = aux;
                        }
                      else // in this case the route is removed
                        {
                          // cout << "remove" << endl;
                          if(u >= WD)
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
  // cout << endl;

  // cout << "Ss:" << endl;
  // for(unsigned i = 0; i < s.Ss.size(); i++)
  //   {
  //     if(i % 2 == 0)
  //       cout << "Depot " << (i / 2) << ": "<< endl;

  //     cout << "\tRoute_" << i % 2 << ": ";
  //     if(s.Ss[i].size())
  //       {
  //         //float temp = 0.0;
  //         for(unsigned j = 0; j < s.Ss[i].size() - 1; j++)
  //           {
  //             cout << s.Ss[i][j] << " ";
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
  // cout << endl;

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
                  float minimum = Infinity;
                  int here = -1, route = 0;
                  w = nodes.back();
                  //cout << "arreglar " << w << endl;
                  for(auto c = D * 2; c < WD * 2; c++)
                    {
                      d = c / 2;
                      // cout << "r" << c << endl;
                      if(s.Sd[d] == 1)
                        {
                          if(s.Ss[c].size() > 0)
                            {
                              for(unsigned r = 0; r < s.Ss[c].size() - 1; r++)
                                {
                                  u = s.Ss[c][r];
                                  v = s.Ss[c][r + 1];

                                  //cout << u << " " << v << "-> ";
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

  // for(auto i = 0; i < N; i++)
  //   cout << s.coveredArea[i] << " ";
  // cout << endl;

  //cout << computeCost(s) << endl;
  computeCost(s);
  //print_sol(s);
  // cout << "End initial Solution" << endl;

  return true;
}

float ILS::coveringSolver(Solution &s, float add)
{
  IloEnv env;
  int iter = 1;
  // All locations are fixed to 1
  for(auto i = 0; i < WD; i++)
    s.Sd[i] = 1;

  // covering areas
  for(auto i = 0; i < WD; i++)
    if(s.Sd[i])
      coveringAreaByNode(s, i);

  vector<vector<int> > nodesByAreas(N);

  // Nodes subset according to uncovered areas
  for(auto i = 0; i < N; i++)
    {
      if(!s.coveredArea[i])
        {
          for(unsigned j = 0; j < bi[i].size(); j++)
            {
              int node = bi[i][j];
              if(node >= WD && !s.noInSs[node])
                {
                  // s.noInSs[node] = true;
                  nodesByAreas[i].push_back(node);
                }
            }
        }
    }
  int count = 0;
  for(auto i = 0; i < N; i++)
    if(nodesByAreas[i].size())
      count++;

  // for(auto i = 0; i < N; i++)
  //   {
  //     cout << "Area " << i << ": ";
  //     for(unsigned j = 0; j < nodesByAreas[i].size(); j++)
  //       {
  //         cout << nodesByAreas[i][j] << " ";
  //       }
  //     cout << endl;
  //   }
  try
    {
      IloModel model(env);
      stringstream name;

      IloNumVarArray y(env, V);
      for(auto i = 0; i < V; i++)
        {
          name << "y_" << i;
          y[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
          name.str("");
        }

      IloCplex cplex(model);
      IloExpr fo(env);


      for(auto i = 0; i < V; i++)
        fo += y[i];

      model.add(IloMinimize(env, fo));

      IloRangeArray constraint(env, count);
      for(auto i = 0, k = 0; i < N; i++)
        {
          if(nodesByAreas[i].size())
            {
              IloExpr exp(env);
              for(unsigned j = 0; j < nodesByAreas[i].size(); j++)
                exp += y[nodesByAreas[i][j]];
              name << "constraint3_" << i;
              constraint[k] = IloRange(env, 1, exp, IloInfinity, name.str().c_str());
              name.str("");
              k++;
            }
        }
      model.add(constraint);

      // cplex.exportModel("model2.lp");
      cplex.setParam(IloCplex::Param::Threads, 1);
      cplex.setOut(env.getNullStream());
      cplex.setParam(IloCplex::Param::MIP::Display, 0);

      if(!cplex.solve())
        {
        }
      int V2 = cplex.getObjValue();
      float routingTime = cplex.getTime();
      cout << "1) Solving Covering problem (minimum number of nodes that cover the zone):" << endl;
      cout << "UB = " << V2 << ", LB = " << cplex.getBestObjValue() << ", time = " << routingTime << endl;
      cout << "y: ";
      vector<int> nickname(WD, 0);
      for(auto i = 0; i < WD; i++)
        nickname[i] = i;
      for(auto i = 0; i < V; i++)
        {
          if(cplex.getValue(y[i]) > 0.9)
            {
              nickname.push_back(i);
              cout << i << " ";
            }
        }
      cout << endl;

      V2 = V2 + WD;
      // decision variables x
      IloArray<IloNumVarArray> x(env, V2);
      for(auto i = 0; i < V2; i++)
        {
          x[i] = IloNumVarArray(env, V2);
          for(auto j = 0; j < V2; j++)
            {
              name << "x_" << i << "_" << j;
              x[i][j] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
              name.str("");
            }
        }

      // decision variables w
      IloNumVarArray w(env, V2);
      for(auto i = 0; i < V2; i++)
        {
          name << "w_" << i;
          w[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Int, name.str().c_str());
          name.str("");
        }


      IloModel model2(env);
      IloCplex cplex2(model2);

      IloExpr fo2(env);
      // route cost
      for(auto i = D; i < V2; i++)
        for(auto j = D; j < V2; j++)
          if(i != j)
            fo2 +=  t[nickname[i]][nickname[j]] * x[i][j];

      model2.add(IloMinimize(env, fo2));

      // Drone routing constraints
      // Constraints 7) \Sum_{k in K} y_{dk} <= M_k, for all d in D
      IloRangeArray constraint7(env, WD - D);
      for(auto d = D; d < WD; d++)
      	{
      	  IloExpr exp(env);
          for(auto j = D; j < V2; j++)
            exp += x[d][j];

          if(d < D)
            exp -= 0;
          else
            exp -= ND;
      	  name << "constraint7_" << d;
      	  constraint7[d - D] = IloRange(env, -IloInfinity, exp, 0, name.str().c_str());
      	  name.str("");
      	}
      model2.add(constraint7);

      // Constraints 8) \Sum_{j in V} x_{dj} - \Sum_{j in V} x_{jd} = 0, for all d in D
      IloRangeArray constraint8(env, WD - D);
      for(auto d = D; d < WD; d++)
        {
          IloExpr sum1(env);
          for(auto j = D; j < V2; j++)
            sum1 += x[j][d];

          IloExpr sum2(env);
          for(auto j = D; j < V2; j++)
            sum2 += x[d][j];

          name << "constraint8_" << d;
      	  constraint8[d - D] = IloRange(env, 0, sum2 - sum1, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint8);

      // Constraints 9) \Sum_{j in V} x_{ij} = v_i, for all i in S
      for(auto i = WD; i < V2; i++)
        {
          IloExpr sum(env);
          for(auto j = D; j < V2; j++)
            sum += x[i][j];
          name << "constraint9_" << i;
          model2.add(IloRange(env, 1, sum, 1, name.str().c_str()));
          name.str("");
        }

      // Constraints 10) \Sum_{j in V} x_{ji} = v_i, for all i in S
      for(auto i = WD, c = 0; i < V2; i++, c++)
        {
          IloExpr sum(env);
          for(auto j = D; j < V2; j++)
            sum += x[j][i];
          name << "constraint10_" << i;
          model2.add(IloRange(env, 1, sum, 1, name.str().c_str()));
          name.str("");
        }

      // decision variables g
      IloArray<IloNumVarArray> g(env, V2);

      for(auto i = 0; i < V2; i++)
        {
          g[i] = IloNumVarArray(env, V2);
          for(auto j = 0; j < V2; j++)
            if(i != j)
              {
                name << "g_" << i << "_" << j;
                g[i][j] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, name.str().c_str());
                name.str("");
              }
        }

      // Constraints 11) \Sum_{j in V} g_{ij} - \Sum_{j in V} g_{ji} == \Sum_{j in V} t_{ij} * x_{ij}, for all i in S,x i != j
      IloRangeArray  constraint11(env, V2 - WD);
      for(auto i = WD, c = 0; i < V2; i++, c++)
        {
          IloExpr sum1(env);
          for(auto j = D; j < V2; j++)
            if(i != j)
              sum1 += g[i][j];

          IloExpr sum2(env);
          for(auto j = D; j < V2; j++)
            if(i != j)
              sum2 += g[j][i];

          IloExpr sum3(env);
          for(auto j = D; j < V2; j++)
            if(i != j)
              sum3 += t[nickname[i]][nickname[j]] * x[i][j];

          name << "constraint11_" << i;
          constraint11[c] = IloRange(env, 0, sum1 - sum2 - sum3, 0, name.str().c_str());
          name.str("");
        }
      model2.add(constraint11);

      IloArray<IloRangeArray> constraint12a(env, V2);
      for(auto i = WD; i < V2; i++)
        {
          constraint12a[i] = IloRangeArray(env, V2);
          for(auto j = WD; j < V2; j++)
            {
              if(i != j)
                {
                  name << "constraint12_UB" << i << "_" << j;
                  constraint12a[i][j] = IloRange(env, -IloInfinity, g[i][j] - (T + add) * x[i][j], 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint12a[i][j]);
                }
            }
        }

      IloArray<IloRangeArray> constraint12b(env, V2);
      for(auto i = WD; i < V2; i++)
        {
          constraint12b[i] = IloRangeArray(env, V2);
          for(auto j = WD; j < V2; j++)
            {
              if(i != j)
                {
                  name << "constraint12_LB_" << i << "_" << j;
                  constraint12b[i][j] = IloRange(env, 0, g[i][j] - x[i][j] * t[nickname[i]][nickname[j]], IloInfinity, name.str().c_str());
                  name.str("");
                  model2.add(constraint12b[i][j]);
                }
            }
        }

      IloArray<IloRangeArray> constraint13(env, WD);
      for(auto i = D; i < WD; i++)
        {
          constraint13[i] = IloRangeArray(env, V2);
          for(auto j = D; j < V2; j++)
            {
              if(i != j)
                {
                  name << "constraint13_LB_" << i << "_" << j;
                  constraint13[i][j] = IloRange(env, 0, g[i][j] - x[i][j] * t[nickname[i]][nickname[j]], 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint13[i][j]);
                }
            }
        }

      IloArray<IloRangeArray> constraint13b(env, V2);
      for(auto i = D; i < V2; i++)
        {
          constraint13b[i] = IloRangeArray(env, WD);
          for(auto j = D; j < WD; j++)
            {
              if(i != j)
                {
                  name << "constraint13_UB_" << i << "_" << j;
                  constraint13b[i][j] = IloRange(env, -IloInfinity, g[i][j] - x[i][j] * (T + add), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint13b[i][j]);
                }
            }
        }


      // // // Constraints 14) \Sum_{j in S} beta_{ji} * v_j >= P_i
      // // IloRangeArray constraint14(env, N);
      // // for(auto i = 0; i < N; i++)
      // //   {
      // //     IloExpr exp(env);
      // //     for(auto j = WD; j < V; j++)
      // //       exp += b[j][i] * v[j - WD];
      // //     exp -= 1;
      // //     name << "constraint14_" << i;
      // //     constraint14[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
      // //     name.str("");
      // //   }
      // // model2.add(constraint14);

      // Constraints 15) w_d = d, for all d in D
      IloRangeArray constraint15(env, WD);
      for(auto d = 0; d < WD; d++)
        {
          name << "constraint15_" << d;
      	  constraint15[d] = IloRange(env, d, w[d], d, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint15);

      // Constraints 16) w_i +- w_j < = (D - 1)  * (1 - x_{ji}), for all i in S, j in V, i != j
      IloArray<IloRangeArray> constraint16(env, V2);
      for(auto i = D; i < V2; i++)
        {
          constraint16[i] = IloRangeArray(env, V2);
          for(auto j = D; j < V2; j++)
            {
              if(i != j)
                {
                  name << "constraint16_" << i << " " << j;
                  constraint16[i][j] = IloRange(env, -IloInfinity, w[i] - w[j] - (WD - 1) * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint16[i][j]);
                }
            }
        }

      // Constraints 17) w_j +- w_i < = (D - 1)  * (1 - x_{ji}), for all i in S, j in V, i != j
      IloArray<IloRangeArray> constraint17(env, V2);
      for(auto i = D; i < V2; i++)
        {
          constraint17[i] = IloRangeArray(env, V2);
          for(auto j = D; j < V2; j++)
            {
              if(i != j)
                {
                  name << "constraint17_" << i << " " << j;
                  constraint17[i][j] = IloRange(env, -IloInfinity, w[j] - w[i] - (WD - 1) * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint17[i][j]);
                }
            }
        }

      // Constraints 18) x_{ij} = 0, for all i in D, j in D
      IloArray<IloRangeArray> constraint18(env, WD);
      for(auto i = D; i < WD; i++)
        {
          constraint18[i] = IloRangeArray(env, WD);
          for(auto j = D; j < WD; j++)
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
      // // // Constraints 19) L_i + P_i <= 1, for all i in N
      // // IloRangeArray constraint19(env, N);
      // // for(auto i = 0; i < N; i++)
      // //   {
      // //     name << "constraint19_" << i;
      // // 	  constraint19[i] = IloRange(env, -IloInfinity, L[i] + P[i], 1, name.str().c_str());
      // // 	  name.str("");
      // //   }
      // // model2.add(constraint19);

      for(auto i = D; i < V2; i++)
        model2.add(x[i][i] == 0);

      // cplex2.exportModel("model3.lp");
      cplex2.setParam(IloCplex::Param::Threads, 1);
      cplex2.setParam(IloCplex::Param::TimeLimit, timeCplex);
      cplex2.setOut(env.getNullStream());
      cplex2.setParam(IloCplex::Param::MIP::Display, 0);
      int routingCost = 0;
      cout << "2) Solving the multi-depot routing problem with minimum number of nodes with sigma = " << add << endl;
      if(!cplex2.solve())
        {
          do {
            iter += 1;
            add += 0.1;
            cout << "UB = " << routingCost << ", LB = " << cplex2.getBestObjValue() << ", time = " << routingTime << endl;
            cout << "can't solve, we tried with sigma " << add << endl;
            cout << "2) Solving the multi-depot routing problem with minimum number of nodes with sigma = " << add << endl;

            for(auto i = D; i < V2; i++)
              {
                //constraint13b[i] = IloRangeArray(env, WD);
                for(auto j = D; j < WD; j++)
                  {
                    if(i != j)
                      {
                        // name << "constraint13_UB_" << i << "_" << j;
                        // constraint13b[i][j] = IloRange(env, -IloInfinity, g[i][j] - x[i][j] * (T + add), 0, name.str().c_str());
                        // name.str("");
                        model2.remove(constraint13b[i][j]);
                      }
                  }
              }

            for(auto i = D; i < V2; i++)
              {
                constraint13b[i] = IloRangeArray(env, WD);
                for(auto j = D; j < WD; j++)
                  {
                    if(i != j)
                      {
                        name << "constraint13_UB_" << i << "_" << j;
                        constraint13b[i][j] = IloRange(env, -IloInfinity, g[i][j] - x[i][j] * (T + add), 0, name.str().c_str());
                        name.str("");
                        model2.add(constraint13b[i][j]);
                      }
                  }
              }
            if(add > 2.0)
              {
                env.end();
                return iter;
              }
          }
          while(!cplex2.solve());

          routingCost = cplex2.getObjValue();
          routingTime = cplex2.getTime();
          cout << "UB = " << routingCost << ", LB = " << cplex2.getBestObjValue() << ", time = " << routingTime << endl;
          env.end();
          return iter;
        }
      routingCost = cplex2.getObjValue();
      routingTime = cplex2.getTime();
      cout << "UB = " << routingCost << ", LB = " << cplex2.getBestObjValue() << ", time = " << routingTime << endl;

      for(auto i = D; i < V2; i++)
        {
          for(auto j = D; j < V2; j++)
            {
              if(cplex2.getValue(x[i][j]) > 0.9)
                {
                  cout << "x'[" << i << "][" << j << "] = " << "x[" << nickname[i] << "][" << nickname[j] << "] = "
                       << cplex2.getValue(x[i][j]) << " -> " << t[nickname[i]][nickname[j]] << endl;
                }
            }
        }

      vector<bool> check(V, true);
      // for(auto c = 0; c < WD * 2; c++)
      //   s.Ss[c].clear();
      int route = 0;
      for(auto i = D; i < WD; i++)
        {
          for(auto j = WD; j < V2; j++)
            {
              // cout << i << " " << j << endl;
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
                          if(v < WD)
                            {
                              // cout << "(" << u << ", " << v << ") " << endl;
                              // cout << "*" << endl;

                              s.Ss[route].push_back(v);
                              //i = route;

                              // route++;
                              break;
                            }
                          u = v;
                          v = D - 1;
                        }
                      v++;
                    }
                }
            }
        }
      for(auto d = D * 2; d < WD * 2; d++)
        {
          if(s.Ss[d].size())
            {
              for(unsigned i = 0; i < s.Ss[d].size() - 1; i++)
                {
                  s.Ss[d][i] = nickname[s.Ss[d][i]];
                }
            }
        }
      // computeCost(s);
      // print_sol(s);
      env.end();
      return iter;
    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
    exit(0);
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
  return iter;
}

float ILS::solver(vector<int> &initial, float currentUB, float timeLimit, Solution &s, bool flag, bool initSol, bool original, int tries, bool fixed)
{
  cout << "Solving the MILP model ..." << endl;
  cout << "Sd used:" << endl;

  for(auto i = 0; i < WD; i++)
    cout << initial[i] << " ";
  cout << endl;
  if(initSol)
    {
      for(auto i = 0; i < WD; i++)
        initial[i] = s.Sd[i];

      // at least one wt actived
      bool atLeast = true;
      int d = 0;
      for(auto i = D; i < WD; i++)
        {
          if(s.Sd[i])
            {
              atLeast = false;
              break;
            }
        }
      if(atLeast)
        {
          d = rnd(D, WD - 1);
          s.Sd[d] = 1;
        }

      cout << "new Sd used" << endl;
      for(auto i = 0; i < WD; i++)
        cout << initial[i] << " ";
      cout << endl;
    }

  IloEnv env;
  try
    {
      stringstream name;

      // decision variables L
      IloBoolVarArray L(env, N);
      for(auto i = 0; i < N; i++)
        {
          name << "L_" << i;
          L[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
          name.str("");
        }
      // decision variables y
      IloNumVarArray y(env, WD);
      if(original)
        {
          for(auto d = 0; d < WD; d++)
            {
              name << "y_" << d;
              y[d] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
              name.str("");
            }
        }
      else
        {
          for(auto d = 0; d < WD; d++)
            {
              name << "y_" << d;
              if(initial[d])
                y[d] = IloNumVar(env, 1, 1, IloNumVar::Bool, name.str().c_str());
              else
                y[d] = IloNumVar(env, 0, 0, IloNumVar::Bool, name.str().c_str());
              name.str("");
            }
          //temporal
        }

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

      IloModel model2(env);
      IloCplex cplex2(model2);

      IloExpr fo(env);

      // route cost
      for(auto i = D; i < V; i++)
        for(auto j = D; j < V; j++)
          if(i != j)
            fo += costUAV * t[i][j] * x[i][j];

      // drones cost
      for(auto d = D; d < WD; d++)
        for(auto j = WD; j < V; j++)
          fo += C[d][3] * x[d][j];

      // ballons cost
      for(auto d = 0; d < D; d++)
        // if(y[d])
        //   fo += C[d][1];
        fo += y[d] * C[d][1];

      // tower cost
      for(auto d = D + We; d < WD; d++)
        fo += y[d] * C[d][0];

      model2.add(IloMinimize(env, fo));

      // Original function objective is became to constraint
      // Constraints 2) \Sum_{i in N} L_i + L_i == N
      IloRangeArray constraint2(env, 1);
      IloExpr exp(env);
      for(auto i = 0; i < N; i++)
        exp += L[i] + P[i];
      name << "constraint2";
      constraint2[0] = IloRange(env, N, exp, N, name.str().c_str());
      name.str("");
      model2.add(constraint2);

      // Constraints 3) \Sum_{d in D} \Sum_{k in K} alpha_{dki} * y_{dk} >= L_i
      IloRangeArray constraint3(env, N);
      for(auto i = 0; i < N; i++)
        {
          IloExpr exp(env);
          for(auto d = 0; d < WD; d++)
            exp += b[d][i] * y[d];
          exp -= L[i];
          name << "constraint3_" << i;
          constraint3[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model2.add(constraint3);

      //Constraints 4) \Sum_{d in D} \Sum_{k in K} y_{dk} >= m_k, for all k in K
      if(We > 0)
        {
          IloRangeArray constraint4(env, We);
          IloExpr exp1(env);
          for(auto d = D, c = 0; d < D + We; d++, c++)
            {
              exp1 = y[d];
              name << "constraint4_0_" << d;
              constraint4[c] = IloRange(env, 1, exp1, 1, name.str().c_str());
              name.str("");
              exp1.clear();
            }
          model2.add(constraint4);
        }

      // Drone routing constraints
      // Constraints 7) \Sum_{k in K} y_{dk} <= M_k, for all d in D
      IloRangeArray constraint7(env, WD - D);
      for(auto d = D; d < WD; d++)
      	{
      	  IloExpr exp(env);
          for(auto j = D; j < V; j++)
            exp += x[d][j];

          if(d < D)
            exp -= 0;
          else
            exp -= ND * y[d];
      	  name << "constraint7_" << d;
      	  constraint7[d - D] = IloRange(env, -IloInfinity, exp, 0, name.str().c_str());
      	  name.str("");
      	}
      model2.add(constraint7);

      // Constraints 8) \Sum_{j in V} x_{dj} - \Sum_{j in V} x_{jd} = 0, for all d in D
      IloRangeArray constraint8(env, WD - D);
      for(auto d = D; d < WD; d++)
        {
          IloExpr sum1(env);
          for(auto j = D; j < V; j++)
            sum1 += x[j][d];

          IloExpr sum2(env);
          for(auto j = D; j < V; j++)
            sum2 += x[d][j];

          name << "constraint8_" << d;
      	  constraint8[d - D] = IloRange(env, 0, sum2 - sum1, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint8);

      // Constraints 9) \Sum_{j in V} x_{ij} = v_i, for all i in S
      IloRangeArray constraint9(env, S);
      for(auto i = WD, c = 0; i < V; i++, c++)
        {
          IloExpr sum(env);
          for(auto j = D; j < V; j++)
            sum += x[i][j];
          sum -= v[c];
          name << "constraint9_" << i;
      	  constraint9[c] = IloRange(env, 0, sum, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint9);

      // Constraints 10) \Sum_{j in V} x_{ji} = v_i, for all i in S
      IloRangeArray constraint10(env, S);
      for(auto i = WD, c = 0; i < V; i++, c++)
        {
          IloExpr sum(env);
          for(auto j = D; j < V; j++)
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
      for(auto i = WD, c = 0; i < V; i++, c++)
        {
          IloExpr sum1(env);
          for(auto j = D; j < V; j++)
            if(i != j)
              sum1 += g[i][j];

          IloExpr sum2(env);
          for(auto j = D; j < V; j++)
            if(i != j)
              sum2 += g[j][i];

          IloExpr sum3(env);
          for(auto j = D; j < V; j++)
            if(i != j)
              sum3 += t[i][j] * x[i][j];

          name << "constraint11_" << i;
          constraint11[c] = IloRange(env, 0, sum1 - sum2 - sum3, 0, name.str().c_str());
          name.str("");
        }
      model2.add(constraint11);

      IloArray<IloRangeArray> constraint12a(env, V);
      for(auto i = WD; i < V; i++)
        {
          constraint12a[i] = IloRangeArray(env, V);
          for(auto j = WD; j < V; j++)
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
      for(auto i = WD; i < V; i++)
        {
          constraint12b[i] = IloRangeArray(env, V);
          for(auto j = WD; j < V; j++)
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

      IloArray<IloRangeArray> constraint13(env, WD);
      for(auto i = D; i < WD; i++)
        {
          constraint13[i] = IloRangeArray(env, V);
          for(auto j = D; j < V; j++)
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
      for(auto i = D; i < V; i++)
        {
          constraint13b[i] = IloRangeArray(env, WD);
          for(auto j = D; j < WD; j++)
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
          for(auto j = WD; j < V; j++)
            exp += b[j][i] * v[j - WD];
          exp -= P[i];
          name << "constraint14_" << i;
          constraint14[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model2.add(constraint14);

      // Constraints 15) w_d = d, for all d in D
      IloRangeArray constraint15(env, WD);
      for(auto d = 0; d < WD; d++)
        {
          name << "constraint15_" << d;
      	  constraint15[d] = IloRange(env, d, w[d], d, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint15);

      // Constraints 16) w_i +- w_j < = (D - 1)  * (1 - x_{ji}), for all i in S, j in V, i != j
      IloArray<IloRangeArray> constraint16(env, V);
      for(auto i = D; i < V; i++)
        {
          constraint16[i] = IloRangeArray(env, V);
          for(auto j = D; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint16_" << i << " " << j;
                  constraint16[i][j] = IloRange(env, -IloInfinity, w[i] - w[j] - (WD - 1) * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint16[i][j]);
                }
            }
        }

      // Constraints 17) w_j +- w_i < = (D - 1)  * (1 - x_{ji}), for all i in S, j in V, i != j
      IloArray<IloRangeArray> constraint17(env, V);
      for(auto i = D; i < V; i++)
        {
          constraint17[i] = IloRangeArray(env, V);
          for(auto j = D; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint17_" << i << " " << j;
                  constraint17[i][j] = IloRange(env, -IloInfinity, w[j] - w[i] - (WD - 1) * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint17[i][j]);
                }
            }
        }

      // Constraints 18) x_{ij} = 0, for all i in D, j in D
      IloArray<IloRangeArray> constraint18(env, WD);
      for(auto i = D; i < WD; i++)
        {
          constraint18[i] = IloRangeArray(env, WD);
          for(auto j = D; j < WD; j++)
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
      	  constraint19[i] = IloRange(env, -IloInfinity, L[i] + P[i], 1, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint19);

      for(auto i = D; i < V; i++)
        model2.add(x[i][i] == 0);

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
          for(auto d = D * 2; d < WD * 2; d++)
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

          for(auto u = D; u < V; u++)
            {
              for(auto v = D; v < V; v++)
                {
                  startVar.add(x[u][v]);
                  startVal.add(xSol[u][v]);
                }
              xSol[u].end(); // free
            }
          xSol.end(); // free

          cplex2.addMIPStart(startVar, startVal, IloCplex::MIPStartAuto, "MIPStart");

          // exit(0);

        }
      // fixed variables
      if(fixed)
        {
          vector<float> values(V, -1);
          int nfixed = 0, ntotal = 0;
          for(unsigned i = D; i < s.noInSs.size(); i++)
            {
              if(!s.noInSs[i])
                {
                  values[i] = randValue(generator);
                  ntotal++;
                }
              //cout << i << " " << values[i] << endl;
            }
          for(auto i = D; i < V; i++)
            {
              if(values[i] > 0 && values[i] < wp)
                {
                  nfixed++;
                  for(auto j = D; j < V; j++)
                    model2.add(x[i][j] == 0);
                }
            }
          cout << "There are " << nfixed << " / " << ntotal << " variables fixed " << endl;
        }

      // Cplex Parameters
      //timeLimit = 10;
      // cplex2.exportModel("model2.lp");
      cplex2.setParam(IloCplex::Param::Threads, 1);
      cplex2.setParam(IloCplex::Param::TimeLimit, timeLimit);
      cplex2.setParam(IloCplex::Param::Emphasis::MIP, 1);
      cplex2.setParam(IloCplex::Param::MIP::Limits::Solutions, tries);
      // cplex2.setParam(IloCplex::Param::Benders::Strategy, 3);
      // cplex2.use(timeLimitCallback(env, cplex2, IloFalse, cplex2.getCplexTime(), timeLimit, currentUB));
      cplex2.setOut(env.getNullStream());
      cplex2.setParam(IloCplex::Param::MIP::Display, 0);
      // Solve
      if(!cplex2.solve())
        {
          cout << "UB = " << currentUB << endl;
          cout << "can't solve'" << endl;
          env.end();
          return Infinity;
        }
      float routingCost = cplex2.getObjValue();
      float routingTime = cplex2.getTime();
      cout << "\tUB = " << routingCost << endl;
      cout << "\tLB = " << cplex2.getBestObjValue() << endl;
      cout << "\ttime = " << routingTime << endl;

      if(flag)
        {
          vector<bool> check(V, true);
          for(auto c = 0; c < WD * 2; c++)
            s.Ss[c].clear();

          s.Sd.clear();

          int route = 0;
          for(auto i = 0; i < WD; i++)
            s.Sd.push_back(cplex2.getValue(y[i]));

          // cout << "y: ";
          // for(auto i = 0; i < WD; i++)
          //   cout << cplex2.getValue(y[i]) <<  " ";
          // cout << endl;
          // cout << "Sd: ";
          // for(auto i = 0; i < WD; i++)
          //   cout << s.Sd[i] << " ";
          // cout << endl;


          // cout << "x: " << endl;
          // for(auto i = D; i < V; i++)
          //   {
          //     for(auto j = D; j < V; j++)
          //       {
          //         if(cplex2.getValue(x[i][j]) > 0.9)
          //           {
          //             cout << "x[" << i << "][" << j << "] = "<< cplex2.getValue(x[i][j]) << " -> " << t[i][j] << endl;
          //           }
          //       }
          //   }

          for(auto i = D; i < WD; i++)
            {
              for(auto j = WD; j < V; j++)
                {
                  // cout << i << " " << j << endl;
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
                              if(v < WD)
                                {
                                  // cout << "(" << u << ", " << v << ") " << endl;
                                  // cout << "*" << endl;

                                  s.Ss[route].push_back(v);
                                  //i = route;

                                  // route++;
                                  break;
                                }
                              u = v;
                              v = D - 1;
                            }
                          v++;
                        }
                    }
                }
            }
        }

      // cout << "P: ";
      // for(auto i = 0; i < N; i++)
      //   cout << cplex2.getValue(P[i]) << " ";
      // cout << endl;

      // cout << "L: ";
      // for(auto i = 0; i < N; i++)
      //   cout << cplex2.getValue(L[i]) << " ";
      // cout << endl;

      // cout << "x: " << endl;
      // float routeCost = 0.0;
      // for(auto i = D; i < V; i++)
      //   {
      //     for(auto j = D; j < V; j++)
      //       {
      //         if(cplex2.getValue(x[i][j]) > 0.9)
      //           {
      //             routeCost += t[i][j];
      //             // cout << "x[" << i << "][" << j << "] = "<< cplex2.getValue(x[i][j]) << " -> " << t[i][j] << endl;
      //           }
      //       }
      //   }

      // cout << "v: ";
      // for(auto i = 0; i < S; i++)
      //   cout << cplex2.getValue(v[i]) << " ";
      // cout << endl;

      // float B2 = 0.0;
      float finalCost = cplex2.getObjValue();
      // print_sol(s);

      // cout << "Budget = " << B2 << endl;
      // cout << "Route cost  = " << routeCost << endl;
      // cout << "Total cost" << B2 + routeCost << " " << finalCost <<endl;

      // if(finalCost < M)
      //   {
      //     // cout << "M = " << M << endl;
      //     M = finalCost;
      //     // cout << "updated M = " << M << endl;
      //   }
      cout << "MILP model ends" << endl << endl;
      env.end();
      return finalCost;
      // cout << "Total time  = " << locationTime + routingTime << endl;
    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
    exit(0);
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
  return -1;
}

void ILS::coveringAreaByNode(Solution &s, int i)
{ // depot [0, WD]
  int cover = -1;
  for(unsigned k = 0; k < beta[i].size(); k++)
    {
      cover = beta[i][k];
      s.coveredArea[cover]++;
    }
}

// exchange two nodes randomly (one from Ss and another one from noInSs).
void ILS::perturbation2(Solution &s)
{
  vector<int> candidate;
  for(auto i = WD; i < V; i++)
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

      int i = rnd(0, candidate.size() - 1);
      int r = rnd(1, s.Ss[route].size() - 2);

      s.Ss[route].insert(s.Ss[route].begin() + r, candidate[i]);
      s.noInSs[candidate[i]] = true;
    }
}

// exchange two nodes randomly from Ss.
void ILS::perturbation3(Solution &s)
{
  int route = 0;
  do {
    route = rnd(0, s.Ss.size() - 1);
  } while(s.Ss[route].size() == 0);

  if(s.Ss[route].size() > 3)
    {
      int i = rnd(1, s.Ss[route].size() - 2);
      int j = rnd(1, s.Ss[route].size() - 2);
      int aux = s.Ss[route][i];

      s.Ss[route][i] = s.Ss[route][j];
      s.Ss[route][j] = aux;
    }
}

void ILS::perturbation4(Solution &s)
{
  int count = 0;
  vector<int> r;

  for(auto i = D * 2; i < WD * 2 && count < 2; i++)
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

      int i = rnd(1, s.Ss[route1].size() - 2);
      int j = rnd(1, s.Ss[route2].size() - 2);

      int aux = s.Ss[route1][i];
      s.Ss[route1][i] = s.Ss[route2][j];
      s.Ss[route2][j] = aux;
    }
}

void ILS::perturbation5(Solution &s)
{
  vector<int> index;
  int d = 0;
  int minK[2] = {0, 0};
  for(unsigned i = D * 2; i < s.Ss.size(); i += 2)
    {
      d = i / 2;
      if(s.Sd[d] == 1 && s.Ss[i].size())
        {
          minK[0]++;
          if(s.Ss[i].size() < 5 && s.Ss[i + 1].size() < 5)
            index.push_back(i);
        }
    }
  if(!index.size() || minK[0] < 2)
    return;

  int j = rnd(0, index.size() - 1);
  int route = index[j];
  d = route / 2;

  vector<int> nodes;
  for(unsigned i = 1; i < s.Ss[route].size() - 1; i++)
    nodes.push_back(s.Ss[route][i]);

  s.Ss[route].clear();
  if(s.Ss[route + 1].size())
    {
      for(unsigned i = 1; i < s.Ss[route + 1].size() - 1; i++)
        nodes.push_back(s.Ss[route + 1][i]);
      s.Ss[route + 1].clear();
    }

  for(unsigned i = 0; i < beta[d].size(); i++)
    {
      s.coveredArea[beta[d][i]]--;
    }

  shuffle(nodes.begin(), nodes.end(), default_random_engine(seed));
  s.Sd[d] = 0;
  int w = 0, u = 0, v = 0;
  for(unsigned i = 0; i < nodes.size(); i++)
    {
      w = nodes[i];
      float minimum = Infinity;
      int here = -1;
      for(auto c = D * 2; c < WD * 2; c++)
        {
          d = c / 2;
          // cout << "r" << c << endl;
          if(s.Sd[d] == 1)
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
          s.Ss[route].insert(s.Ss[route].begin() + here + 1, w);
          coveringAreaByNode(s, w);
        }
    }
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
              if(node >= WD && !mark[node])
                {
                  mark[node] = true;
                  missing.push_back(node);
                }
            }
        }
    }

  if(missing.size())
    {
      for(auto ii = 0; ii < N; ii++)
        {
          if(!s.coveredArea[ii])
            {
              for(auto i = missing.size() - 1; missing.size(); i--)
                {
                  int u, v, w;
                  float minimum = Infinity;
                  int here = -1, route = 0;
                  w = missing.back();
                  for(auto c = D * 2; c < WD * 2; c++)
                    {
                      d = c / 2;
                      // cout << "r" << c << endl;
                      if(s.Sd[d] == 1)
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
                      s.Ss[route].insert(s.Ss[route].begin() + here + 1, w);
                      coveringAreaByNode(s, w);
                      s.noInSs[w] = true;
                      missing.pop_back();
                      ii = 0;
                      break;
                    }
                }
            }
        }
    }
}
void ILS::computeRoute(Solution &s, vector<float> &c)
{
  c = vector<float>(WD * 2, 0);
  for(unsigned i = 0; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size())
        {
          for(unsigned j = 1; j < s.Ss[i].size(); j++)
            c[i] += t[s.Ss[i][j - 1]][s.Ss[i][j]];
        }
    }
}

bool ILS::LS_swap_outnodes(Solution &s)
{
  for(unsigned i = D * 2; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size())
        {
          for(unsigned j = 1; j < s.Ss[i].size() - 1; j++)
            {
              float currentCost = t[s.Ss[i][j - 1]][s.Ss[i][j]] + t[s.Ss[i][j]][s.Ss[i][j + 1]];
              for(int out = WD; out < V; out++)
                {
                  if(!s.noInSs[out])
                    {
                      float newCost = t[s.Ss[i][j - 1]][out] + t[out][s.Ss[i][j + 1]];
                      if(currentCost - newCost > 0)
                        {
                          //cout << i << " " << j << " " << s.Ss[i][j] << ", " << out << endl;
                          s.noInSs[out] = true;
                          s.noInSs[s.Ss[i][j]] = false;
                          s.Ss[i][j] = out;
                        }
                    }
                }
            }
        }
    }
  return true;
}

bool ILS::LS_swap(Solution &s)
{
  vector<float> Cost;
  computeRoute(s, Cost);
  bool flag = false;
  for(unsigned i = D * 2; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size())
        {
          for(unsigned j = 1; j < s.Ss[i].size() - 1; j++)
            {
              for(unsigned route = D * 2; route < s.Ss.size(); route++)
                {
                  if(s.Ss[route].size())
                    {
                      for(unsigned l = 1; l < s.Ss[route].size() - 1; l++)
                        {
                          if(i == route && j == l)
                            continue;
                          else
                            {
                              int opt = 0;
                              float a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = 0, h = 0;
                              if(s.Ss[route][l] != s.Ss[i][j + 1])
                                {
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
                                      if(i == route && (e - a +  h - d < 0))
                                        opt = 1;
                                      if(i != route && (e + f - a - b + Cost[i] < T) && (g + h - c - d + Cost[route] < T))
                                        opt = 1;
                                    }
                                  else
                                    {
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
                                  a = t[s.Ss[i][j - 1]][s.Ss[i][j]];
                                  b = t[s.Ss[route][l]][s.Ss[route][l + 1]];
                                  c = t[s.Ss[i][j - 1]][s.Ss[route][l]];
                                  d = t[s.Ss[i][j]][s.Ss[route][l + 1]];

                                  if(- a - b + c + d < 0)
                                    opt = 1;
                                }

                              if(opt)
                                {
                                  int aux = s.Ss[i][j];
                                  s.Ss[i][j] = s.Ss[route][l];
                                  s.Ss[route][l] = aux;

                                  computeCost(s);
                                  computeRoute(s, Cost);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
  return flag;
}

bool ILS::removeNode(Solution &s)
{
  if(!s.feasible)
    return false;
  int node = 0, area = 0;
  int route = 0;

  do {
    route = rnd(0, s.Ss.size() - 1);
  } while(s.Ss[route].size() == 0);
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
                }
              // cout << "que " << i << endl;
              s.Ss[route].erase(s.Ss[route].begin() + i); // remove node
              s.noInSs[node] = false; // add to list noInSs
              i--;
            }
        }
    }

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
  if(s.size() < 3)
    return false;
  float diff = 0;
  bool improve = false;
  // int depot = s[0];
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
  // if(depot != s[0])
  //   {
  //     cout << "must fix the depot at begin" << endl;
  //     exit(0);
  //   }

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
  for(auto i = 0; i < WD; i++)
    {
      if(i < D && s.Sd[i])
        {
          s.Lcost += C[i][1]; // ballon cost
          coveringAreaByNode(s, i);
        }
      else if(i >= D && i < D + We && s.Sd[i])
        {
          s.Lcost += C[i][2]; // wte cost
          coveringAreaByNode(s, i);
        }
      else if(s.Sd[i])
        {
          s.Lcost += C[i][0]; // wt cost
          coveringAreaByNode(s, i);
        }
    }

  //cout << "\nLcost: "<< s.Lcost << endl;
  for(auto i = WD; i < V; i++)
    s.noInSs[i] = false;

  // compute drones' routes
  int u = 0, v = 0;
  float routeCost = 0.0;
  float routeCostReal = 0.0;
  s.Rcost = 0.0;
  s.feasible = 1;

  for(unsigned i = D * 2; i < s.Ss.size(); i++)
    {
      if(s.Ss[i].size() > 0 && s.Sd[i / 2] != 1)
        {
          // cout << i << " " << s.Ss[i].size() << " " << s.Sd[i / 2]  << endl;
          // to do: memory leak with cplex
          // cout << "error this route is not watchtower" << endl;
          // print_sol(s);
          s.Sd[i / 2] = 1;
          s.feasible = 0;
          //exit(0);
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
              routeCostReal +=  costUAV * t[u][v];
              // cout << "covering node " << v << endl;
              coveringAreaByNode(s, v);
              s.noInSs[u] = 1;
              s.noInSs[v] = 1;
              // if(s.Ss[i][j + 1] < WD)
              //   {
              //     print_sol(s);
              //     cout << "error!!!!!!!!!!!!!!!!!!!!! it's a depot'" << endl;
              //     exit(0);
              //   }
            }
          u = s.Ss[i][s.Ss[i].size() - 2];
          v = s.Ss[i][s.Ss[i].size() - 1];
          s.noInSs[u] = 1;
          s.noInSs[v] = 1;

          // // Remove it after!
          // if(v != s.Ss[i][0])
          //   {
          //     print_sol(s);
          //     cout << "error!!!!!!!!!!!!!!!!!!!!!" << endl;
          //     exit(0);
          //   }

          // cout << u << " " << v << endl;
          routeCost += t[u][v];
          routeCostReal +=  costUAV * t[u][v];
          if(routeCost > T)
            s.feasible = 0;
          s.Rcost += routeCostReal;
        }
    }
  // final cost
  s.cost = s.Lcost + s.Rcost;
  // check if the solution is feasible
  for(auto i = 0; i < N; i++)
    if(!s.coveredArea[i])
      s.feasible = 0;

  // if(routeCost > 10000)
  //   {
  //     print_sol(s);
  //     cout << "error big number !!!!!!!!!!!!!!!!!!!!!" << endl;
  //     exit(0);
  //   }

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
  for(auto i = 0; i < WD; i++)
    {
      cout << s.Sd[i] << " ";
      if(i < D && s.Sd[i])
        cost += C[i][1]; // ballon cost
      else if(i >= D && i < D + We && s.Sd[i])
        cost += C[i][2]; // wte cost
      else if(s.Sd[i])
        cost += C[i][0]; // wt cost
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
          cost += costUAV * temp + C[0][3];
        }
      cout << endl;
      //cout << i << endl;
    }
  cout << "noInSs = ";
  for(auto i = WD; i < V; i++)
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

void ILS::updateElite(Solution &s1, float cost)
{
  if(s1.feasible)
    {
      if(elite.size() > 0)
        {
          map<float, Solution>::iterator rit = --elite.end();
          if(elite.size() < E || rit->first > cost)
            {
              if(elite.find(cost) == elite.end())
                {
                  elite.insert(make_pair(cost, s1));
                  //cout << "elite size" << elite.size() << endl;
                  while(elite.size() >= E)
                    elite.erase(--elite.end());
                }
            }
        }
      else
        elite.insert(make_pair(cost, s1));
    }
}

void ILS::runMH(int timeLimitAlgorithm)
{
  nameAlgorithm = "Simple-Matheuristic";
  cout << "Algorithm time limit = " << timeLimitAlgorithm << endl;
  cout << endl << nameAlgorithm << " starts" << endl;
  int iterMax = 10;
  if(timeLimitAlgorithm != -1)
    iterMax = 999999999;
  else
    timeLimitAlgorithm = 999999999;

  float timelimit = 10;

  Solution s(N, V, WD);
  Solution sBest = Solution();
  vector<int> initial(WD, 1);
  // build initial solution
  cout << "Initial solution starts" << endl;
  int nTimes = coveringSolver(s, 0);

  float currentCost = computeCost(s);
  cout << "Initial solution needed: " << nTimes << " iterations" << endl;
  cout << "Solution found is " << ((s.feasible == 1) ? "feasible" : "infeasible") << endl;
  if(s.feasible)
    cout << "Initial Cost: " << currentCost << endl;
  else
    {
      currentCost = Infinity;
      currentCost = solver(initial, currentCost, 1000, s, true, false, true, 1, false);
      currentCost = computeCost(s);
      cout << "Initial Cost: " << currentCost << endl;
    }

  initialCost = costFinal = currentCost;
  cout << "Initial solution ends" << endl;

  auto start = chrono::steady_clock::now();
  for(auto iter = 0; iter < iterMax; iter++)
    {
      auto end = chrono::steady_clock::now();
      auto diff = end - start;
      auto totalTime  = chrono::duration <double, std::ratio<1>> (diff).count();
      // cout << totalTime << endl;
      if(totalTime > timeLimitAlgorithm)
        break;

      vector<int> Sd(WD);
      for(auto i = 0; i < WD; i++)
        Sd[i] = s.Sd[i];

      currentCost = solver(s.Sd, currentCost, timelimit, s, true, true, true, 2, true);

      if(currentCost < costFinal)
        {
          costFinal = currentCost;
          updateElite(s, s.cost);
        }
    }

  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  auto totalTime  = chrono::duration <double, std::ratio<1>> (diff).count();
  cout << "Cost = " << costFinal << endl;
  cout << "Time = " << totalTime << endl;
  timeF = totalTime;
}

void ILS::run(int timeLimitAlgorithm, int variant)
{
  nameAlgorithm = "Matheuristic";
  cout << "Algorithm time limit = " << timeLimitAlgorithm << endl;
  cout << endl << nameAlgorithm << " starts" << endl;

  if(timeLimitAlgorithm != -1)
    Kt = 999999999;
  else
    timeLimitAlgorithm = 999999999;

  if(variant)
    Km = 999999999;

  auto start = chrono::steady_clock::now();

  Solution s(N, V, WD);
  Solution sBest = Solution();
  Solution sBestBest = Solution();
  // build initial solution
  cout << "Initial solution starts" << endl;
  int nTimes = coveringSolver(s, 0);

  float s_cost = computeCost(s);
  vector<int> initial(WD, 1);
  cout << "Initial solution needed: " << nTimes << " iterations" << endl;
  cout << "Solution found is " << ((s.feasible == 1) ? "feasible" : "infeasible") << endl;
  if(s.feasible)
    cout << "Initial Cost: " << s_cost << endl;
  else
    {
      s_cost = Infinity;
      s_cost = solver(initial, s_cost, 1000, s, true, false, true, 1, false);
      s_cost = computeCost(s);
      cout << "Initial Cost: " << s_cost << endl;
    }
  cout << "Initial solution ends" << endl;

  float s_bestCost = Infinity;
  float s_bestBestCost = Infinity;
  if(s.feasible)
    {
      s_bestBestCost = s_bestCost = s_cost;
      //cout << sBest.coveredArea.size();
      sBest = s;
      sBestBest = sBest;
    }

  cout << endl << "Matheuristic starts"<< endl;
  initialCost = s_bestBestCost;
  float increaseTime = Tc;
  float alpha = 0;

  for(unsigned iter = 1; iter < Kt; iter++)
    {
      auto end = chrono::steady_clock::now();
      auto diff = end - start;
      auto totalTime  = chrono::duration <double, std::ratio<1>> (diff).count();
      // cout << totalTime << endl;
      if(totalTime > timeLimitAlgorithm)
        break;
      if(iter % 20000 == 0)
        cout << "iter " << iter << endl;
      // Perturbations
      alpha = randValue(generator);
      if(alpha < p1 && elite.size()) // perturbation1
        {
          map<float, Solution>::iterator it = elite.begin();
          int position = rnd(0, elite.size() - 1);
          advance(it, position);
          s = it->second;
        }
      else if(alpha < p2)
        perturbation2(s);
      else if(alpha < p3)
        perturbation3(s);
      else if(alpha < p4)
        perturbation4(s);
      else // p5
        perturbation5(s);
      s_cost = computeCost(s);
      // Local searches
      for(auto c = D * 2; c < WD * 2; c++)
        {
          // 2-opt for each route
          if(s.Ss[c].size())
            twoOpt(s.Ss[c]);
        }
      s_cost = computeCost(s);
      removeNode(s);
      s_cost = computeCost(s);
      LS_swap(s);
      s_cost = computeCost(s);
      LS_swap_outnodes(s);
      s_cost = computeCost(s);
      // keep best solution
      if(s.feasible && (s_bestCost > s_cost))
        {
          s_bestCost = s_cost;
          sBest = s;
        }
      else if(iter % Kl == 0)  // Local Reset
        {
          initialSolution(s);
          s_cost = computeCost(s);
        }
      else if(!s.feasible && elite.size()) // Acceptance criterion: reject solution if is infeasible
        {
          map<float, Solution>::iterator it = elite.begin();
          int position = rnd(0, elite.size() - 1);
          // cout << "rnd " << position << endl;
          advance(it, position);
          s = it->second;
          s_cost = it->first;
          s_bestCost = s_cost; //mirar esto!!!!!
        }
      // Solving the MILP model
      if(iter % Km == 0)
        {
          cout << "BestBestCost " << s_bestBestCost << endl;
          cout << "BestCost " << s_bestCost << endl;
          cout << "MIP best" << endl;
          map<float, Solution>::iterator it;
          if(elite.size() > 0)
            {
              it = elite.begin();
              s = it->second;
              for(unsigned i = 0; i < initial.size(); i++)
                initial[i] = s.Sd[i];
              alpha = solver(initial, s_bestBestCost, increaseTime, s, true, true, true, 2, true);
            }
          if(alpha < s_bestCost)
            increaseTime = Tc;
          else
            increaseTime += 2;

          s_cost = computeCost(s);
        }
      else if(iter % Kg == 0) // Global reset
        {
          cout << "Intensification" << endl;
          // Intensification
          map<float, Solution> current_elite(elite);
          for(map<float, Solution>::iterator it = current_elite.begin(); it != current_elite.end(); ++it)
            {
              s = it->second;
              // 2-opt for each route
              for(auto c = D * 2; c < WD * 2; c++)
                {
                  if(s.Ss[c].size())
                    twoOpt(s.Ss[c]);
                }
              s_cost = computeCost(s);
              removeNode(s);
              s_cost = computeCost(s);
              LS_swap(s);

              if(s.feasible && (s_bestCost > s_cost))
                {
                  updateElite(s, s.cost);
                  s_bestCost = s_cost;
                  sBest = s;
                }
            }
          if(s_bestBestCost > s_bestCost)
            {
              s_bestBestCost = s_bestCost;
              sBestBest = sBest;
              cout << iter << " " << s_bestCost << endl;
            }
          if(elite.size())
            {
              map<float, Solution>::iterator it = elite.begin();
              int position = rnd(0, elite.size() - 1);
              // cout << "rnd " << position << endl;
              advance(it, position);
              s = it->second;
              s_cost = it->first;
              s_bestCost = s_cost;
            }
          else
            {
              initialSolution(s);
              s_cost = computeCost(s);
              s_bestCost = Infinity;
            }
          sBest = s;
        }
    }
  cout << "Best solution:" << endl;
  auto end = chrono::steady_clock::now();
  auto diff = end - start;
  auto totalTime  = chrono::duration <double, std::ratio<1>> (diff).count();

  s_bestBestCost = computeCost(sBestBest);
  cout << "Best cost " << setprecision(decipresi) << s_bestBestCost << endl;
  // print_sol(sBestBest);
  cout << "Best cost " << setprecision(decipresi) << s_bestBestCost << endl;

  if(elite.size())
    {
      cout << "Elite solution:" << endl;
      for(map<float, Solution>::iterator it = elite.begin(); it != elite.end(); ++it)
        cout << setprecision(decipresi) << it->first << endl;

      map<float, Solution>::iterator it = elite.begin();
      // print_sol(it->second);
      costFinal = it->first;
      timeF = totalTime;
    }
  else
    {
      cout << "Feasible solution is not found." << endl;

      costFinal = Infinity;
      timeF = totalTime;
    }
}
