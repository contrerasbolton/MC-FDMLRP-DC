/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "ils.hpp"

// Globals variables (parameters)
int N;
int D;
int K;
int S;
int B;
int ND;
float T;
int V;
int *mMax;
int *mMin;
int **C;
int **beta;
int ***alpha;
float **t;
float distSum;
string nameInstance;
vector<vector<vector<int> > > a; // simple alpha
vector<vector<int> > b; // simple beta
vector<vector<int> > bi; // inverse beta
void printError(string msg)
{
  cout << msg << endl;
  exit(0);
}

void readInstance(const char *instance)
{
  FILE *file;
  char temp[100];
  if((file = fopen(instance, "r")) == NULL)
    printError("Error in reading of the instance " + string(instance));

  /* reading n */
  if(!fscanf(file, "N=%d;\n", &N))
    printError("reading error in N");
  /* reading d */
  if(!fscanf(file, "D=%d;\n", &D))
    printError("reading error in D");
  /* reading k */
  if(!fscanf(file, "K=%d;\n", &K))
    printError("reading error in K");
  /* reading s */
  if(!fscanf(file, "S=%d;\n", &S))
    printError("reading error in s");
  /* reading B */
  if(!fscanf(file, "B=%d;\n", &B))
    printError("reading error in B");

  mMax = new int[K];
  mMin = new int[K];
  /* reading min m_k */
  if(!fscanf(file, "M=[%d %d];\n", &mMin[0], &mMin[1]))
    printError("reading error in m Min");
  /* reading min m_k */
  if(!fscanf(file, "MM=[%d %d];\n", &mMax[0], &mMax[1]))
    printError("reading error in m Max");
  /* reading ND */
  if(!fscanf(file, "ND=%d;\n", &ND))
    printError("reading error in ND");
  /* reading T */
  if(!fscanf(file, "T=%f;\n", &T))
    printError("reading error in T");
  /* reading C */
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in C");

  cout << "N = " << N << endl;
  cout << "D = " << D << endl;
  cout << "k = " << K << endl;
  cout << "S = " << S << endl;
  cout << "B = " << B << endl;
  cout << "m min = " << mMin[0] << " " << mMin[1] << endl;
  cout << "m max = " << mMax[0] << " " << mMax[1] << endl;
  cout << "ND = " << ND << endl;
  cout << "T = " << T << endl;

  cout << temp << endl;
  C = new int*[D];
  for(auto i = 0; i < D; i++)
    {
      C[i] = new int[2];
      if(!fscanf(file, "%d\t%d\n", &C[i][0], &C[i][1]))
        printError("reading error in C_" + to_string(i));
      cout << C[i][0] << " " <<  C[i][1] << endl;
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in C");
  cout << "It is loaded" << endl;
  cout << temp << endl;
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << temp << endl;

  beta = new int*[S];
  for(auto i = 0; i < D; i++)
    b.push_back(vector<int>());

  for(auto i = 0; i < S; i++)
    {
      beta[i] = new int[N];
      vector<int> aux;
      cout << i + D << " -> ";
      for(auto j = 0; j < N; j++)
        {
          if(!fscanf(file, "%d", &beta[i][j]))
            printError("reading error in beta " + to_string(i) + " " + to_string(j));
          if(beta[i][j])
            {
              aux.push_back(j);
              cout << j << " ";
            }
        }
      b.push_back(aux);
      cout << endl;
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << "It is loaded" << endl;
  cout << temp << endl;

  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << temp << endl;

  alpha = new int**[D];
  for(auto d = 0; d < D; d++)
    {
      vector<vector<int>> aux;
      alpha[d] = new int*[K];
      for(auto k = 0; k < K; k++)
        {
          vector<int> aux2;
          alpha[d][k] = new int[N];
          for(auto i = 0; i < N; i++)
            {
              if(!fscanf(file, "%d", &alpha[d][k][i]))
                printError("reading error in alpha " + to_string(d) + " " + to_string(k) + " " + to_string(i));
              if(alpha[d][k][i])
                {
                  aux2.push_back(i);
                  cout << d << " " << k << " " << i << endl;
                }
            }
          aux.push_back(aux2);
          //cout << endl;
        }
      a.push_back(aux);
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << "It is loaded" << endl;
  cout << temp << endl;

  if(!fscanf(file, "%s\n", temp))
    printError("reading error in t");
  cout << temp << endl;

  V = D + S;
  t = new float*[V];
  distSum = 0.0;
  for(auto i = 0; i < V; i++)
    {
      t[i] = new float[V];
      for(auto j = 0; j < V; j++)
        {
          if(i != j)
            {
              if(!fscanf(file, "%f", &t[i][j]))
                printError("reading error in t " + to_string(i) + " " + to_string(j));
              distSum += t[i][j];
            }
          else
            t[i][j] = 99999;
          cout << t[i][j] << " ";
        }
      cout << endl;
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in t");
  cout << "It is loaded" << endl;
  cout << temp << endl;
  fclose(file);
  // compute betaInv
  for(auto i = 0; i < N; i++)
    bi.push_back(vector<int>());

  for(unsigned i = 0; i < b.size(); i++)
    for(unsigned j = 0; j < b[i].size(); j++)
      bi[b[i][j]].push_back(i + D);

  for(unsigned i = 0; i < bi.size(); i++)
    {
      cout << i << " -> ";

      for(unsigned j = 0; j < bi[i].size(); j++)
        cout << bi[i][j] << " ";
      cout << endl;
    }
}

void solveMILP(int opt)
{
  IloEnv env;
  try
    {
      IloModel model(env);
      stringstream name;

      // Definition of the decision variables
      // decision variables L
      IloBoolVarArray L(env, N);
      for(auto i = 0; i < N; i++)
        {
          name << "L_" << i;
          L[i] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
          name.str("");
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

      // decision variables u
      IloNumVarArray u(env, V);
      // decision variables g
      IloArray<IloNumVarArray> g(env, V);

      if(opt % 2 == 0)
        {
          for(auto i = 0; i < V; i++)
            {
              name << "u_" << i;
              u[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Float, name.str().c_str());
              name.str("");
            }
        }
      else
        {
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
        }
      // decision variables w
      IloNumVarArray w(env, V);
      for(auto i = 0; i < V; i++)
        {
          name << "w_" << i;
          w[i] = IloNumVar(env, 0, IloInfinity, IloNumVar::Int, name.str().c_str());
          name.str("");
        }

      // decision variables y
      IloArray<IloNumVarArray> y(env, D);
      for(auto d = 0; d < D; d++)
        {
          y[d] = IloNumVarArray(env, N);
          for(auto k = 0; k < K; k++)
            {
              name << "y_" << d << "_" << k;
              y[d][k] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
              name.str("");
            }
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

      IloCplex cplex(model);
      if(opt == 0 || opt == 1) // Model 1
        {
          // Function objetive: 1) \Sum_{i in N} L_i + P_i
          IloExpr fo(env);
          for(auto i = 0; i < N; i++)
            fo += L[i] + P[i];
          model.add(IloMaximize(env, fo));
        }
      else if(opt == 2 || opt == 3) // Model 2
        {
          // Another function objetive: min route cost
          IloExpr fo(env);
          for(auto i = 0; i < V; i++)
            for(auto j = 0; j < V; j++)
              if(i != j)
                fo += t[i][j] * x[i][j];

          model.add(IloMinimize(env, fo));
        }
      else // Model 3
        {
          // Another function objective: min route cost
          IloExpr fo(env);
          for(auto i = 0; i < V; i++)
            for(auto j = 0; j < V; j++)
              if(i != j)
                fo += t[i][j] * x[i][j];

          for(auto d = 0; d < D; d++)
            for(auto k = 0; k < K; k++)
              fo += y[d][k] * C[d][k];
          model.add(IloMinimize(env, fo));
        }

      // Location-covering constraints
      // Constraints 2) \Sum_{d in D} \Sum_{k in K} alpha_{dki} * y_{dk} >= L_i
      IloRangeArray constraint2(env, N);
      for(auto i = 0; i < N; i++)
        {
          IloExpr exp(env);
          for(auto d = 0; d < D; d++)
            for(auto k = 0; k < K; k++)
              exp += alpha[d][k][i] * y[d][k];
          exp -= L[i];
          name << "constraint2_" << i;
          constraint2[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model.add(constraint2);

      if(opt == 0 || opt == 1)
        {
          // Constraints 3) \Sum_{d in D} \Sum_{k in K}  y_{dk} <= B, for all i in N
          IloRangeArray constraint3(env, 1);
          IloExpr exp(env);
          for(auto d = 0; d < D; d++)
            for(auto k = 0; k < K; k++)
              exp += y[d][k] * C[d][k];
          name << "constraint3";
          constraint3[0] = IloRange(env, -IloInfinity, exp, B, name.str().c_str());
          name.str("");
          model.add(constraint3);
        }
      else
        {
          // original function objective is became to constraint
          // Constraints 3) \Sum_{i in N} L_i + L_i == N
          IloRangeArray constraint3(env, 1);
          IloExpr exp(env);
          for(auto i = 0; i < N; i++)
            exp += L[i] + P[i];
          name << "constraint3";
          constraint3[0] = IloRange(env, N, exp, N, name.str().c_str());
          name.str("");
          model.add(constraint3);
        }

      // Constraints 4) \Sum_{d in D} \Sum_{k in K} y_{dk} >= m_k, for all k in K
      IloRangeArray constraint4(env, K);
      for(auto k = 0; k < K; k++)
        {
          IloExpr exp(env);
          for(auto d = 0; d < D; d++)
            exp += y[d][k];
          name << "constraint4_" << k;
          constraint4[k] = IloRange(env, mMin[k], exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model.add(constraint4);

      // Constraints 5) \Sum_{d in D} y_{dk} <= M_k, for all k in K
      IloRangeArray constraint5(env, K);
      for(auto k = 0; k < K; k++)
      	{
      	  IloExpr exp(env);
      	  for(auto d = 0; d < D; d++)
      	    exp += y[d][k];
      	  name << "constraint5_" << k;
      	  constraint5[k] = IloRange(env, -IloInfinity, exp, mMax[k], name.str().c_str());
      	  name.str("");
      	}
      model.add(constraint5);

      // Constraints 6) \Sum_{k in K} y_{dk} <= M_k, for all d in D
      IloRangeArray constraint6(env, D);
      for(auto d = 0; d < D; d++)
      	{
      	  IloExpr exp(env);
          for(auto k = 0; k < K; k++)
      	    exp += y[d][k];
      	  name << "constraint6_" << d;
      	  constraint6[d] = IloRange(env, -IloInfinity, exp, 1, name.str().c_str());
      	  name.str("");
      	}
      model.add(constraint6);

      // Drone routing constraints
      // Constraints 7) \Sum_{k in K} y_{dk} <= M_k, for all d in D
      IloRangeArray constraint7(env, D);
      for(auto d = 0; d < D; d++)
      	{
      	  IloExpr exp(env);
          for(auto j = 0; j < V; j++)
            exp += x[d][j];
          exp -= ND * y[d][0];
      	  name << "constraint7_" << d;
      	  constraint7[d] = IloRange(env, -IloInfinity, exp, 0, name.str().c_str());
      	  name.str("");
      	}
      model.add(constraint7);

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
      model.add(constraint8);

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
      model.add(constraint9);

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
      model.add(constraint10);

      // MTZ
      if(opt % 2 == 0)
        {
          // Constraints 11) u_i + t_{ij} < = u_j + Q * (1 - x_{ji}), for all i in S, j in V, i != j
          IloArray<IloRangeArray>  constraint11(env, S);
          for(auto i = D, c = 0; i < V; i++, c++)
            {
              constraint11[c] = IloRangeArray(env, V);
              for(auto j = 0; j < V; j++)
                {
                  if(i != j)
                    {
                      name << "constraint11_" << i << " " << j;
                      constraint11[c][j] = IloRange(env, -IloInfinity, u[i] - u[j] + t[i][j] - V * (1 - x[i][j]), 0, name.str().c_str());
                      name.str("");
                      model.add(constraint11[c][j]);
                    }
                }
            }

          // Constraints 12) u_j >= \Sum_{j in D} x_{ij} * t_{ij}, for all j in S
          IloRangeArray constraint12(env, S);
          for(auto j = 0; j < S; j++)
            {
              IloExpr sum(env);
              for(auto d = 0; d < D; d++)
                sum += x[d][j + D] * t[d][j + D];
              sum = -1 * sum + u[j + D];
              name << "constraint12_" << j;
              constraint12[j] = IloRange(env, 0, sum, IloInfinity, name.str().c_str());
              name.str("");
            }
          model.add(constraint12);

          // Constraints 13) u_i <= T, for all i in D
          IloRangeArray constraint13(env, D);
          for(auto d = 0; d < D; d++)
            {
              name << "constraint13_" << d;
              constraint13[d] = IloRange(env, -IloInfinity, u[d], T, name.str().c_str());
              name.str("");
            }
          model.add(constraint13);
        }
      else // GG
        {
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
          model.add(constraint11);

          IloArray<IloRangeArray> constraint12(env, V);
          for(auto i = 0; i < V; i++)
            {
              constraint12[i] = IloRangeArray(env, V);
              for(auto j = 0; j < V; j++)
                {
                  if(i != j)
                    {
                      name << "constraint12_" << i << "_" << j;
                      constraint12[i][j] = IloRange(env, -IloInfinity, g[i][j] - (distSum / 2.0) * x[i][j], 0, name.str().c_str());
                      name.str("");
                      model.add(constraint12[i][j]);
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
                    // model.add(g[i][j] >= x[i][j] * t[i][j]);
                    name << "constraint13_LB_" << i << "_" << j;
                    constraint13[i][j] = IloRange(env, 0, g[i][j] - x[i][j] * t[i][j], IloInfinity, name.str().c_str());
                    name.str("");
                    model.add(constraint13[i][j]);
                  }
                }
            }

          IloArray<IloRangeArray> constraint13b(env, D);
          for(auto i = 0; i < D; i++)
            {
              constraint13b[i] = IloRangeArray(env, V);
              for(auto j = 0; j < V; j++)
                {
                  if(i != j)
                  {
                    name << "constraint13_UB_" << i << "_" << j;
                    constraint13b[i][j] = IloRange(env, -IloInfinity, g[j][i], T, name.str().c_str());
                    name.str("");
                    model.add(constraint13b[i][j]);
                  }
                }
            }
        }

      // Constraints 14) \Sum_{j in S} beta_{ji} * v_j >= P_i
      IloRangeArray constraint14(env, N);
      for(auto i = 0; i < N; i++)
        {
          IloExpr exp(env);
          for(auto j = 0; j < S; j++)
            exp += beta[j][i] * v[j];
          exp -= P[i];
          name << "constraint14_" << i;
          constraint14[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model.add(constraint14);

      // Constraints 15) w_d = d, for all d in D
      IloRangeArray constraint15(env, D);
      for(auto d = 0; d < D; d++)
        {
          name << "constraint15_" << d;
      	  constraint15[d] = IloRange(env, d, w[d], d, name.str().c_str());
      	  name.str("");
        }
      model.add(constraint15);

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
                  model.add(constraint16[i][j]);
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
                  model.add(constraint17[i][j]);
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
                  model.add(constraint18[i][j]);
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
      model.add(constraint19);

      for(auto i = 0; i < V; i++)
        model.add(x[i][i] == 0);

      // Cplex Parameters
      int timeLimit = 3600;
      // cplex.exportModel("model.lp");
      cplex.setParam(IloCplex::Param::Threads, 1);
      cplex.setParam(IloCplex::Param::TimeLimit, timeLimit);
      cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, 4000);
      // Solve
      if(!cplex.solve())
        {
          cout << "can't solve'" << endl;
          throw(-1);
        }

      cout << "P: ";
      int sumP = 0;
      for(auto i = 0; i < N; i++)
        {
          sumP += cplex.getValue(P[i]);
          cout << cplex.getValue(P[i]) << " ";
        }
      cout << endl;

      cout << "L: ";
      int sumL = 0;
      for(auto i = 0; i < N; i++)
        {
          sumL += cplex.getValue(L[i]);
          cout << cplex.getValue(L[i]) << " ";
        }
      cout << endl;

      cout << "y: " << endl;
      int watchtowers = 0, ballons = 0;
      float B2 = 0.0;
      for(auto i = 0; i < D; i++)
        {
          for(auto j = 0; j < K; j++)
            {
              B2 += cplex.getValue(y[i][j]) * C[i][j];
              int y_ij = cplex.getValue(y[i][j]);
              cout << "y[" << i << "][" << j << "] = "<<  y_ij << endl;
              if(!j && y_ij)
                watchtowers++;
              else if(j && y_ij)
                ballons++;
            }
        }

      cout << "B = " << B2 << endl;
      float routeCost = 0.0;
      cout << "x: " << endl;
      for(auto i = 0; i < V; i++)
        {
          for(auto j = 0; j < V; j++)
            {
              if(cplex.getValue(x[i][j]) > 0.9)
                {
                  routeCost += t[i][j];
                  cout << "x[" << i << "][" << j << "] = "<< cplex.getValue(x[i][j]) << " -> " << t[i][j] << endl;
                }
            }
        }

      cout << "v: ";
      for(auto i = 0; i < S; i++)
        cout << cplex.getValue(v[i]) << " ";
      cout << endl;

      if(opt % 2 == 0)
        {
          cout << "u: ";
          for(auto i = 0; i < V; i++)
            cout << cplex.getValue(u[i]) << " ";
          cout << endl;
        }
      else
        {
          cout << "g: " << endl;
          for(auto i = 0; i < V; i++)
            {
              for(auto j = 0; j < V; j++)
                {
                  if(i != j && cplex.getValue(g[i][j]) > 0.9)
                    {
                      cout << "g[" << i << "][" << j << "] = "<< cplex.getValue(g[i][j]) << endl;
                    }
                }
            }
        }

      cout << "w: ";
      for(auto i = 0; i < V; i++)
        cout << cplex.getValue(w[i]) << " ";
      cout << endl;

      float totalTime = cplex.getTime();
      float totalCost = routeCost + B2;
      int totalCover = sumL + sumP;
      float gap = ((cplex.getStatus() == IloAlgorithm::Optimal) ? 0 : cplex.getMIPRelativeGap() * 100);
      IloAlgorithm::Status status = cplex.getStatus();

      cout << "Instance    = " << nameInstance << endl;
      cout << "Total Cost  = " << totalCost << endl;
      cout << "Route Cost  = " << routeCost << endl;
      cout << "Budget      = " << B2 << " / " << B << endl;
      cout << "Cover       = " << totalCover << endl;
      cout << "Cover L     = " << sumL << endl;
      cout << "Cover P     = " << sumP << endl;
      cout << "Watchtowers = " << watchtowers << endl;
      cout << "Ballons     = " << ballons << endl;
      cout << "Time        = " << totalTime << endl;
      cout << "Status      = " << status << endl;
      cout << "GAP (%)     = " << gap << endl;

      FILE *file;
      stringstream ss;
      string s;
      string output = "summary_" + to_string(opt) + ".txt";
      ss << status;
      ss >> s;
      if((file = fopen(output.c_str(), "a")) == NULL)
        {
          printf("Error in reading of the %s \n", output.c_str());
          exit(0);
        }
      fprintf(file, "%s\t%d\t%d\t%d\t%d\t%s\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",
              nameInstance.c_str(), N, S, D, B, s.c_str(), totalCost, routeCost, B2, totalCover, sumL, sumP, watchtowers, ballons, gap, totalTime);

    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
}

void callILS()
{
  cout << "ILS is running" << endl;
  int seed = 0;
  ILS *ils = new ILS(seed, N, D, K, S, B, ND, T, V, a, b, bi, t, mMax, mMin, C);
  ils->run();
  delete ils;
}

int main(int argc, char *argv[])
{
  cout << "An heuristic for MC-FDMLRP-DC" << endl;
  string folder = "instances/";
  nameInstance = argv[1];
  string name = folder + nameInstance;
  int opt = atoi(argv[2]);

  cout << name << endl;
  readInstance(name.c_str());

  if(opt == 6)
    callILS();
  else
    solveMILP(opt);

  return 0;
}
