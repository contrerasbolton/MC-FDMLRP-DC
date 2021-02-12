/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include <iostream>
#include <string>

#include <ilcplex/cplex.h>
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN

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
  cout << temp << endl;
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << temp << endl;

  beta = new int*[S];
  for(auto i = 0; i < S; i++)
    {
      beta[i] = new int[N];
      for(auto j = 0; j < N; j++)
        {
          if(!fscanf(file, "%d", &beta[i][j]))
            printError("reading error in beta " + to_string(i) + " " + to_string(j));
          cout << beta[i][j] << " ";
        }
      cout << endl;
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << temp << endl;

  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << temp << endl;

  alpha = new int**[D];
  for(auto d = 0; d < D; d++)
    {
      alpha[d] = new int*[K];
      for(auto k = 0; k < K; k++)
        {
          alpha[d][k] = new int[N];
          for(auto i = 0; i < N; i++)
            {
              if(!fscanf(file, "%d", &alpha[d][k][i]))
                printError("reading error in alpha " + to_string(d) + " " + to_string(k) + " " + to_string(i));
              cout << alpha[d][k][i] << " ";
            }
          cout << endl;
        }
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in beta");
  cout << temp << endl;

  if(!fscanf(file, "%s\n", temp))
    printError("reading error in t");
  cout << temp << endl;

  V = D + S;
  t = new float*[V];
  for(auto i = 0; i < V; i++)
    {
      t[i] = new float[V];
      for(auto j = 0; j < V; j++)
        {
          if(i != j)
            {
              if(!fscanf(file, "%f", &t[i][j]))
                printError("reading error in t " + to_string(i) + " " + to_string(j));
            }
          else
            t[i][j] = 99999;

          cout << t[i][j] << " ";
        }
      cout << endl;
    }
  if(!fscanf(file, "%s\n", temp))
    printError("reading error in t");
  cout << temp << endl;


  fclose(file);
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
      if(opt == 0)
        {
          // Function objetive: 1) \Sum_{i in N} L_i + P_i
          IloExpr fo(env);
          for(auto i = 0; i < N; i++)
            fo += L[i] + P[i];
          model.add(IloMaximize(env, fo));
        }
      else
        {
          // Another function objetive: min route cost
          IloExpr fo(env);
          for(auto i = 0; i < V; i++)
            {
              for(auto j = 0; j < V; j++)
                {
                  if(i != j)
                    {
                      fo += t[i][j] * x[i][j];
                    }
                }
            }
          model.add(IloMinimize(env, fo));

          // original function objetive is becamed to constraint
          IloExpr expsum(env);
          for(auto i = 0; i < N; i++)
            expsum += L[i] + P[i];
          model.add( expsum == N);
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

      // Constraints 3) \Sum_{d in D} \Sum_{k in K}  y_{dk} <= B, for all i in N
      IloRangeArray constraint3(env, 1);
      IloExpr exp2(env);
      for(auto d = 0; d < D; d++)
        for(auto k = 0; k < K; k++)
          exp2 += y[d][k] * C[d][k];
      name << "constraint3";
      constraint3[0] = IloRange(env, -IloInfinity, exp2, B, name.str().c_str());
      name.str("");
      model.add(constraint3);

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

      // Constraints 11) u_i + t_{ij} < = u_j + Q * (1 - x_{ji}), for all i in S, j in V, i != j
      int Q = V;
      IloArray<IloRangeArray>  constraint11(env, S);
      for(auto i = D, c = 0; i < V; i++, c++)
        {
          constraint11[c] = IloRangeArray(env, V);
          for(auto j = 0; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint11_" << i << " " << j;
                  constraint11[c][j] = IloRange(env, -IloInfinity, u[i] - u[j] + t[i][j] - Q * (1 - x[i][j]), 0, name.str().c_str());
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
            sum -= x[d][j + D] * t[d][j + D];
          sum += u[j + D];
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
      	  constraint13[d] = IloRange(env, -IloInfinity, u[d] - T, 0, name.str().c_str());
      	  name.str("");
        }
      model.add(constraint13);

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

      // Constraints 15) u_i <= T, for all i in D
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
        {
          model.add(x[i][i] == 0);
        }

      // Cplex Parameters
      cplex.exportModel("model.lp");

      // Solve
      if(!cplex.solve())
        {
          cout << "can't solve'" << endl;
          throw(-1);
        }


      cout << "P: ";
      for(auto i = 0; i < N; i++)
        cout << cplex.getValue(P[i]) << " ";
      cout << endl;

      cout << "L: ";
      for(auto i = 0; i < N; i++)
        cout << cplex.getValue(L[i]) << " ";
      cout << endl;

      cout << "y: " << endl;
      for(auto i = 0; i < D; i++)
        {
          for(auto j = 0; j < K; j++)
            cout << "y[" << i << "][" << j << "] = "<< cplex.getValue(y[i][j]) << endl;
        }

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

      cout << "u: ";
      for(auto i = 0; i < V; i++)
        cout << cplex.getValue(u[i]) << " ";
      cout << endl;

      cout << "w: ";
      for(auto i = 0; i < V; i++)
        cout << cplex.getValue(w[i]) << " ";
      cout << endl;

      int coverCost = cplex.getObjValue();
      float coverTime = cplex.getTime();
      cout << "Cover = " << coverCost << endl;
      cout << "Cost  = " << routeCost << endl;
      cout << "time  = " << coverTime << endl;


    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
}

void solveTwoStage(int opt)
{
  IloEnv env;
  try
    {
      IloModel model1(env);
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

      // Fist stage: Location
      IloCplex cplex(model1);
      // Function objetive: 1) \Sum_{i in N} L_i
      IloExpr fo(env);
      for(auto i = 0; i < N; i++)
        fo += L[i];
      model1.add(IloMaximize(env, fo));

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
      model1.add(constraint2);

      // Constraints 3) \Sum_{d in D} \Sum_{k in K}  y_{dk} <= B, for all i in N
      IloRangeArray constraint3(env, 1);
      IloExpr exp2(env);
      for(auto d = 0; d < D; d++)
        for(auto k = 0; k < K; k++)
          exp2 += y[d][k] * C[d][k];
      name << "constraint3";
      constraint3[0] = IloRange(env, -IloInfinity, exp2, B, name.str().c_str());
      name.str("");
      model1.add(constraint3);

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
      model1.add(constraint4);

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
      model1.add(constraint5);

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
      model1.add(constraint6);

      // Cplex Parameters
      cplex.exportModel("model1.lp");

      // Solve
      if(!cplex.solve())
        {
          cout << "can't solve'" << endl;
          throw(-1);
        }

      int locationCost = cplex.getObjValue();
      float locationTime = cplex.getTime();
      cout << "FO   = " << locationCost << endl;
      cout << "time = " << locationTime << endl;

      cout << "L: ";
      for(auto i = 0; i < N; i++)
        cout << cplex.getValue(L[i]) << " ";
      cout << endl;

      cout << "y: " << endl;
      for(auto i = 0; i < D; i++)
        {
          for(auto j = 0; j < K; j++)
            cout << "y[" << i << "][" << j << "] = "<< cplex.getValue(y[i][j]) << endl;
        }

      // Second stage
      IloModel model2(env);
      IloCplex cplex2(model2);
      // Function objetive: 1) \Sum_{i in N} P_i
      IloExpr fo2(env);
      for(auto i = 0; i < N; i++)
        fo2 += P[i];
      model2.add(IloMaximize(env, fo2));

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

      // Constraints 11) u_i + t_{ij} < = u_j + Q * (1 - x_{ji}), for all i in S, j in V, i != j
      int Q = V;
      IloArray<IloRangeArray>  constraint11(env, S);
      for(auto i = D, c = 0; i < V; i++, c++)
        {
          constraint11[c] = IloRangeArray(env, V);
          for(auto j = 0; j < V; j++)
            {
              if(i != j)
                {
                  name << "constraint11_" << i << " " << j;
                  constraint11[c][j] = IloRange(env, -IloInfinity, u[i] - u[j] + t[i][j] - Q * (1 - x[i][j]), 0, name.str().c_str());
                  name.str("");
                  model2.add(constraint11[c][j]);
                }
            }
        }

      // Constraints 12) u_j >= \Sum_{j in D} x_{ij} * t_{ij}, for all j in S
      IloRangeArray constraint12(env, S);
      for(auto j = 0; j < S; j++)
        {
          IloExpr sum(env);
          for(auto d = 0; d < D; d++)
            sum -= x[d][j + D] * t[d][j + D];
          sum += u[j + D];
          name << "constraint12_" << j;
      	  constraint12[j] = IloRange(env, 0, sum, IloInfinity, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint12);

      // Constraints 13) u_i <= T, for all i in D
      IloRangeArray constraint13(env, D);
      for(auto d = 0; d < D; d++)
        {
          name << "constraint13_" << d;
      	  constraint13[d] = IloRange(env, -IloInfinity, u[d] - T, 0, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint13);

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
      model2.add(constraint14);

      // Constraints 15) u_i <= T, for all i in D
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
          int Li = cplex.getValue(L[i]);
      	  constraint19[i] = IloRange(env, -IloInfinity, Li + P[i], 1, name.str().c_str());
      	  name.str("");
        }
      model2.add(constraint19);

      for(auto i = 0; i < V; i++)
        {
          model2.add(x[i][i] == 0);
        }

      // Cplex Parameters
      cplex2.exportModel("model2.lp");

      // Solve
      if(!cplex2.solve())
        {
          cout << "can't solve'" << endl;
          throw(-1);
        }

      int routingCost = cplex2.getObjValue();
      float routingTime = cplex2.getTime();
      cout << "FO   = " << routingCost << endl;
      cout << "time = " << routingTime << endl;

      cout << "P: ";
      for(auto i = 0; i < N; i++)
        cout << cplex2.getValue(P[i]) << " ";
      cout << endl;

      cout << "x: " << endl;
      float routeCost = 0.0;
      for(auto i = 0; i < V; i++)
        {
          for(auto j = 0; j < V; j++)
            {
              if(cplex2.getValue(x[i][j]) > 0.9)
                {
                  routeCost += t[i][j];
                  cout << "x[" << i << "][" << j << "] = "<< cplex2.getValue(x[i][j]) << " -> " << t[i][j] << endl;
                }
            }
        }

      cout << "v: ";
      for(auto i = 0; i < S; i++)
        cout << cplex2.getValue(v[i]) << " ";
      cout << endl;

      cout << "u: ";
      for(auto i = 0; i < V; i++)
        cout << cplex2.getValue(u[i]) << " ";
      cout << endl;

      cout << "w: ";
      for(auto i = 0; i < V; i++)
        cout << cplex2.getValue(w[i]) << " ";
      cout << endl;

      cout << "Total Cover = " << locationCost + routingCost << endl;
      cout << "Total Cost  = " << routeCost << endl;
      cout << "Total time  = " << locationTime + routingTime << endl;
    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
}


int main(int argc, char *argv[])
{
  cout << "An heuristic for MC-FDMLRP-DC" << endl;
  string folder = "instances/";
  string name = folder + argv[1];
  int opt = atoi(argv[2]);

  cout << name << endl;
  readInstance(name.c_str());

  if(opt % 2 == 0)
    solveMILP(opt);
  else
    solveTwoStage(opt);
  return 0;
}
