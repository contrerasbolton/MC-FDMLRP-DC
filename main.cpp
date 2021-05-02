/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2; -*- */
#include "ils.hpp"

// Globals variables (parameters)
int N;
int We;
int W;
int D;
int K; // eliminar
int WD;
int S;
int B;
int ND;
float T;
int V;
int *mMax;
int *mMin;
float *techCost;
int *WTe;
float **C;
int **beta;
int **alpha;
float **t;
float distSum;
string nameInstance;
vector<vector<vector<int> > > a; // simple alpha
vector<vector<int> > a2; // simple alpha
vector<vector<int> > b; // simple beta
vector<vector<int> > bi; // inverse beta

// Input parameters
int seed;
int timeLimit;
int memoryLimit;

void printError(string msg)
{
  cout << msg << endl;
  exit(0);
}

void readInstance(const char *instance)
{
  FILE *file;
  //char temp[100];
  int temp = 0;
  if((file = fopen(instance, "r")) == NULL)
    printError("Error in reading of the instance " + string(instance));

  /* reading N D We W S */
  if(!fscanf(file, "%d %d %d %d %d", &N, &D, &We, &W, &S))
    printError("reading error in N, W, D, S");

  WD = D + We + W;

  K = 1;
  mMax = new int[K];
  mMin = new int[K];
  if(!fscanf(file, "%d %f %d %d", &ND, &T, &mMin[0], &mMax[0]))
    printError("reading error in ND, T, Min and Max of watchtowers/ballons ");

  techCost = new float[4]; // Cost = {Watchtowers, Ballons, ?, drone}
  if(!fscanf(file, "%f %f %f %f", &techCost[0], &techCost[1], &techCost[2], &techCost[3]))
    printError("reading error in Costs");

  V = D + We + W + S;

  for(auto i = 0; i < V; i++)
    b.push_back(vector<int>());

  for(auto i = 0; i < V; i++)
    {
      int id = 0;
      int size = 0;

      if(!fscanf(file, "%d %d", &id, &size))
        printError("reading error in Covered areas: size and id");


      for(auto j = 0; j < size; j++)
        {
          if(!fscanf(file, "%d", &temp))
            printError("reading error in Covered areas: size and id");
          b[id].push_back(temp);
        }
    }

  for(auto i = 0; i < WD; i++)
    a2.push_back(vector<int>());

  t = new float*[V];
  for(auto i = 0; i < V; i++)
    {
      t[i] = new float[V];
      for(auto j = 0; j < V; j++)
        {
          if(!fscanf(file, "%f", &t[i][j]))
            printError("reading error in t " + to_string(i) + " " + to_string(j));
        }
    }

  cout << N << " " << D << " " << We << " " << W << " " << S << " " <<endl;
  cout << ND << " " << T << " " << mMin[0] << " " << mMax[0] << endl;
  cout << techCost[0] << " " << techCost[1] << " " << techCost[2] << " " << techCost[3] << endl;
  cout << "Node Type, id and areas that are covered:" << endl;
  for(auto i = 0; i < (int) b.size(); i++)
    {
      if(i < D)
        cout << "D  ";
      else if (i >= D && i < D + We)
        cout << "We ";
      else if (i >= D + We && i < D + We + W)
        cout << "W  ";
      else if (i >= WD)
        cout << "n  ";

      cout << i << ": ";
      for(unsigned j = 0; j < b[i].size(); j++)
        cout << b[i][j] << " ";
      cout << endl;
    }

  // create beta
  beta = new int*[V];
  for(auto i = 0; i < V; i++)
    {
      beta[i] = new int[N];
      for(auto j = 0; j < N; j++)
        beta[i][j] = 0;
      for(unsigned j = 0; j < b[i].size(); j++)
        beta[i][b[i][j]] = 1;
    }

  // compute betaInv
  for(auto i = 0; i < N; i++)
    bi.push_back(vector<int>());

  for(unsigned i = 0; i < b.size(); i++)
    for(unsigned j = 0; j < b[i].size(); j++)
      bi[b[i][j]].push_back(i);

  C = new float*[WD];
  for(auto i = 0; i < WD; i++)
    {
      C[i] = new float[4];
      C[i][0] = techCost[0];
      C[i][1] = techCost[1];
      C[i][2] = techCost[2];
      C[i][3] = techCost[3];
    }

  // for(auto i = 0; i < S; i++)
  //   {
  //     for(auto j = 0; j < N; j++)
  //       cout << beta[i][j] << " ";
  //     cout << endl;
  //   }

  cout << "Area id and its covered nodes:" << endl;
  for(unsigned i = 0; i < bi.size(); i++)
    {
      cout << i << ": ";
      for(unsigned j = 0; j < bi[i].size(); j++)
        cout << bi[i][j] << " ";
      cout << endl;
    }

  // for(auto i = 0; i < V; i++)
  //   {
  //     for(auto j = 0; j < V; j++)
  //       cout << t[i][j] << " ";
  //     cout << endl;
  //   }
  cout << "reading is ok" << endl;
}

void solveMILP(int opt)
{
  float costUAV = 0.02;
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
          for(auto i = 0; i < D; i++)
            {
              name << "u_" << i;
              u[i] = IloNumVar(env, 0, 0, IloNumVar::Float, name.str().c_str());
              name.str("");
            }
          for(auto i = D; i < V; i++)
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
      IloNumVarArray y(env, WD);
      for(auto d = 0; d < WD; d++)
        {
          name << "y_" << d;
          y[d] = IloNumVar(env, 0, 1, IloNumVar::Bool, name.str().c_str());
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
          // Another function objetive: min drones cost + ballons cost + tower cost
          IloExpr fo(env);

          // drones cost
          for(auto d = D; d < WD; d++)
            for(auto j = WD; j < V; j++)
              fo += techCost[3] * x[d][j];

          // ballons cost
          for(auto d = 0; d < D; d++)
            fo += y[d] * techCost[1];

          // tower cost
          for(auto d = D + We; d < WD; d++)
            fo += y[d] * techCost[0];

          model.add(IloMinimize(env, fo));
        }
      else // Model 3
        {
          // Another function objective: min (route cost + drones cost + ballons cost + tower cost)
          IloExpr fo(env);

          // route cost
          for(auto i = D; i < V; i++)
            for(auto j = D; j < V; j++)
              if(i != j)
                fo += costUAV * t[i][j] * x[i][j];

          // drones cost
          for(auto d = D; d < WD; d++)
            for(auto j = WD; j < V; j++)
              fo += techCost[3] * x[d][j];

          // ballons cost
          for(auto d = 0; d < D; d++)
            fo += y[d] * techCost[1];

          // tower cost
          for(auto d = D + We; d < WD; d++)
            fo += y[d] * techCost[0];

          model.add(IloMinimize(env, fo));
        }

      if(opt == 0 || opt == 1)
        {
          B = 200;
          // Constraints 2) \Sum_{d in D} \Sum_{k in K}  y_{dk} <= B, for all i in N
          IloRangeArray constraint2(env, 1);
          IloExpr exp(env);

          // drones cost
          for(auto d = D; d < WD; d++)
            for(auto j = WD; j < V; j++)
              exp += techCost[3] * x[d][j];

          // ballons cost
          for(auto d = 0; d < D; d++)
            exp += y[d] * techCost[1];

          // tower cost
          for(auto d = D + We; d < WD; d++)
            exp += y[d] * techCost[0];

          name << "constraint2";
          constraint2[0] = IloRange(env, -IloInfinity, exp, B, name.str().c_str());
          name.str("");
          model.add(constraint2);
        }
      else
        {
          // Original function objective is became to constraint
          // Constraints 2) \Sum_{i in N} L_i + L_i == N
          IloRangeArray constraint2(env, 1);
          IloExpr exp(env);
          for(auto i = 0; i < N; i++)
            exp += L[i] + P[i];
          name << "constraint2";
          constraint2[0] = IloRange(env, N, exp, N, name.str().c_str());
          name.str("");
          model.add(constraint2);
        }


      // Location-covering constraints
      // Constraints 3) \Sum_{d in D} \Sum_{k in K} alpha_{dki} * y_{dk} >= L_i
      IloRangeArray constraint3(env, N);
      for(auto i = 0; i < N; i++)
        {
          IloExpr exp(env);
          for(auto d = 0; d < WD; d++)
            exp += beta[d][i] * y[d];
          exp -= L[i];
          name << "constraint3_" << i;
          constraint3[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model.add(constraint3);

      // Constraints 4) \Sum_{d in D} \Sum_{k in K} y_{dk} >= m_k, for all k in K
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
          model.add(constraint4);
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
      model.add(constraint7);

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
      model.add(constraint8);

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
      model.add(constraint9);

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
      model.add(constraint10);

      // MTZ
      if(opt % 2 == 0)
        {
          // Constraints 11) u_i + t_{ij} < = u_j + Q * (1 - x_{ji}), for all i in S, j in V, i != j
          IloArray<IloRangeArray>  constraint11(env, V);
          for(auto i = WD, c = 0; i < V; i++, c++)
            {
              constraint11[c] = IloRangeArray(env, V);
              for(auto j = D; j < V; j++)
                {
                  if(i != j)
                    {
                      name << "constraint11_" << i << " " << j;
                      constraint11[c][j] = IloRange(env, -IloInfinity, u[i] - u[j] + t[i][j] - T * (1 - x[i][j]), 0, name.str().c_str());
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
              for(auto d = D; d < WD; d++)
                sum -= x[d][j + WD] * t[d][j + WD];
              sum += u[j + WD];
              name << "constraint12_" << j;
              constraint12[j] = IloRange(env, 0, sum, IloInfinity, name.str().c_str());
              name.str("");
            }
          model.add(constraint12);

          // Constraints 13) u_i <= T, for all i in D
          IloRangeArray constraint13(env, WD - D);
          for(auto d = D; d < WD; d++)
            {
              name << "constraint13_" << d;
              constraint13[d - D] = IloRange(env, -IloInfinity, u[d], T, name.str().c_str());
              name.str("");
            }
          model.add(constraint13);
        }
      else // GG
        {
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
          model.add(constraint11);

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
                      model.add(constraint12a[i][j]);
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
                    model.add(constraint12b[i][j]);
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
                      model.add(constraint13[i][j]);
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
          for(auto j = WD; j < V; j++)
            exp += beta[j][i] * v[j - WD];
          exp -= P[i];
          name << "constraint14_" << i;
          constraint14[i] = IloRange(env, 0, exp, IloInfinity, name.str().c_str());
          name.str("");
        }
      model.add(constraint14);

      // Constraints 15) w_d = d, for all d in D
      IloRangeArray constraint15(env, WD);
      for(auto d = 0; d < WD; d++)
        {
          name << "constraint15_" << d;
      	  constraint15[d] = IloRange(env, d, w[d], d, name.str().c_str());
      	  name.str("");
        }
      model.add(constraint15);

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
                  model.add(constraint16[i][j]);
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
                  model.add(constraint17[i][j]);
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

      for(auto i = D; i < V; i++)
        model.add(x[i][i] == 0);

      // Cplex Parameters
      //cplex.exportModel("model.lp");
      cplex.setParam(IloCplex::Param::Threads, 1);
      cplex.setParam(IloCplex::Param::TimeLimit, timeLimit);
      cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, memoryLimit);
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

      cout << "y: ";
      for(auto i = 0; i < WD; i++)
        cout << cplex.getValue(y[i]) << " ";
      cout << endl;

      int watchtowers = 0, ewatchtowers = 0, ballons = 0, nDrones = 0;
      float wtCost = 0.0, ballCost = 0.0, routeCost = 0.0;
      // count ballons
      for(auto i = 0; i < D; i++)
        {
          if(cplex.getValue(y[i]) > 0.9)
            {
              ballons++;
              ballCost += techCost[1];
            }
        }
      // count existing watchtowers
      for(auto i = D; i < D + We; i++)
        if(cplex.getValue(y[i]) > 0.9)
          ewatchtowers++;
      // count watchtowers
      for(auto i = D + We; i < WD; i++)
        {
          if(cplex.getValue(y[i]) > 0.9)
            {
              watchtowers++;
              wtCost += techCost[0];
            }
        }
      cout << "x: " << endl;
      for(auto i = D; i < V; i++)
        {
          for(auto j = D; j < V; j++)
            {
              if(cplex.getValue(x[i][j]) > 0.9)
                {
                  routeCost += t[i][j] * costUAV;
                  cout << "x[" << i << "][" << j << "] = "<< cplex.getValue(x[i][j]) << " -> " << t[i][j] << endl;
                  if(i < WD)
                    nDrones++;
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
          for(auto i = D; i < V; i++)
            cout << cplex.getValue(u[i]) << " ";
          cout << endl;
        }
      else
        {
          cout << "g: " << endl;
          for(auto i = D; i < V; i++)
            {
              for(auto j = D; j < V; j++)
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

      float LB = cplex.getBestObjValue();
      float totalTime = cplex.getTime();
      float DroneCost = nDrones * techCost[3];

      float totalCost = wtCost + DroneCost + ballCost;

      if(opt > 3)
        totalCost += routeCost;
      int totalCover = sumL + sumP;
      float gap = ((cplex.getStatus() == IloAlgorithm::Optimal) ? 0 : cplex.getMIPRelativeGap() * 100);
      IloAlgorithm::Status status = cplex.getStatus();

      cout << "Instance            = " << nameInstance << endl;
      cout << "LB                  = " << LB << endl;
      cout << "Total Cost          = " << totalCost << " " << cplex.getObjValue() << endl;
      if(opt > 3)
        cout << "Total Cost w/o RC   = " << totalCost - routeCost << endl;
      cout << "Route Cost          = " << routeCost << endl;
      cout << "Watchtowers Cost    = " << wtCost << endl;
      cout << "Ballons Cost        = " << ballCost << endl;
      cout << "Drones Cost         = " << DroneCost << endl;
      cout << "Drone number        = " << nDrones << endl;
      cout << "existing WT number  = " << ewatchtowers << endl;
      cout << "Watchtowers number  = " << watchtowers << endl;
      cout << "Ballons number      = " << ballons << endl;
      cout << "Covering            = " << totalCover << endl;
      cout << "Covering by Facili. = " << sumL << endl;
      cout << "Covering by Drones  = " << sumP << endl;
      cout << "Time                = " << totalTime << endl;
      cout << "Status              = " << status << endl;
      cout << "GAP (%)             = " << gap << endl;

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
      fprintf(file, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n",
              nameInstance.c_str(), N, S, D, We, W, B, s.c_str(), LB, totalCost, routeCost, wtCost, ballCost, DroneCost, totalCover, sumL, sumP, ballons, ewatchtowers,
              watchtowers, nDrones, gap, totalTime);
      fclose(file);

    } catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
  }
  catch (...) {
    cerr << "Error" << endl;
  }
  env.end();
}

void callILS(int opt)
{
  cout << "ILS is running" << endl;
  // initial cost, cost obtained, computing time
  float result[3] = {0.0, 0.0, 0.0};
  ILS *ils = new ILS(seed, N, D, K, S, B, beta, ND, T, V, distSum, a, b, bi, t, mMax, mMin, C);
  ils->run(&result[0]);
  delete ils;
  cout << result[0] << "\t" << result[1] << "\t" << result[2] << endl;
  FILE *file;
  stringstream ss;
  string s;
  string output = "summary_" + to_string(opt) + ".txt";
  if((file = fopen(output.c_str(), "a")) == NULL)
    {
      printf("Error in reading of the %s \n", output.c_str());
      exit(0);
    }
  fprintf(file, "%s\t%f\t%f\t%f\t\n", nameInstance.c_str(), result[0], result[1], result[2]);
  fclose(file);
}

void callMatheuristic()
{
  cout << "Matheuristic is running" << endl;
  ILS *ils = new ILS(seed, N, D, K, S, B, beta, ND, T, V, distSum, a, b, bi, t, mMax, mMin, C);
  ils->runMH();
  delete ils;
}


void solveTwoStage(int opt)
{

}

int main(int argc, char *argv[])
{
  cout << "An heuristic for MC-FDMLRP-DC" << endl;
  string folder = "instances/";
  nameInstance = argv[1];
  string name = folder + nameInstance;
  int opt = atoi(argv[2]);
  seed = atoi(argv[3]);
  timeLimit = atoi(argv[4]);
  if(argc > 5)
    memoryLimit = atoi(argv[5]);
  else
    memoryLimit = 5000;

  cout << name << endl;
  readInstance(name.c_str());

  if(opt == 7)
    callMatheuristic();
  else if(opt == 6)
    callILS(opt);
  else
    solveMILP(opt);

  return 0;
}
