#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<sstream>
#include<ctime>

#include"matplotlibcpp.h"

#include"coptcpp_pch.h"

#include "gurobi_c++.h" 

namespace plt = matplotlibcpp;

using namespace std;

//=======================
//data in use
string filename = "C101";
int N = 100;//number of customer
double capacity;
double* xcoord = nullptr;
double* ycoord = nullptr;
double* demand = nullptr;
double* ready_time = nullptr;
double* due_date = nullptr;
double* service_time = nullptr;

double** dis = nullptr;
double fakeinf;

//solution save
double** sol = nullptr;

//time consume
double copt_time;
double gurobi_time;
//=======================

void Init();

void Release();

void useGurobi();

void useCOPT();

void drawSolution();

string itos(int i) { stringstream s; s << i; return s.str(); }
string dtos(double d) { stringstream s; s << d; return s.str(); }

int main() {

	Init();

	useCOPT();

	useGurobi();

	drawSolution();

	Release();

	return 0;
}

void Init() {
	cout << "init data set ..." << endl;

	N++;

	xcoord = new double[N];
	ycoord = new double[N];
	demand = new double[N];
	ready_time = new double[N];
	due_date = new double[N];
	service_time = new double[N];

	ifstream fin(filename + ".txt", ios::in);
	string str;

	int i, j;

	int len1, len2;

	fakeinf = 0;

	if (filename.size() <= 5) { len1 = 1, len2 = 7; }
	else { len1 = 5, len2 = 12; }

	for (i = 0; i < len1; i++) fin >> str;

	fin >> capacity;

	for (i = 0; i < len2; i++) fin >> str;

	for (i = 0; i < N; i++) {
		fin >> str;
		fin >> xcoord[i] >> ycoord[i] >> demand[i] >> 
			ready_time[i] >> due_date[i] >> service_time[i];

		fakeinf += demand[i] + ready_time[i] + due_date[i] + service_time[i];

	}

	fin.close();

	sol = new double* [N];
	dis = new double* [N];
	for (i = 0; i < N; i++) {
		sol[i] = new double[N];
		dis[i] = new double[N];
		for (j = 0; j < N; j++) {
			sol[i][j] = 0;
			dis[i][j] = sqrt((xcoord[i] - xcoord[j]) * (xcoord[i] - xcoord[j]) 
				+ (ycoord[i] - ycoord[j]) * (ycoord[i] - ycoord[j]));

			fakeinf + dis[i][j];
		}
	}
}

void Release() {

	ofstream fout("record.csv", ios::app);

	fout << filename << ',' << N - 1 << ',' << copt_time << ',' << gurobi_time << endl;

	fout.close();

	delete[] xcoord;
	delete[] ycoord;
	delete[] demand;
	delete[] ready_time;
	delete[] due_date;
	delete[] service_time;

	for (int i = 0; i < N; i++) delete[] sol[i];
	delete[] sol;

	cout << "clear data set ..." << endl;
}

void useGurobi() {
	cout << "use Gurobi here ..." << endl;

	clock_t start_time = clock();

	try {
		GRBEnv* env = new GRBEnv();
		GRBVar** vars = new GRBVar * [N];
		GRBVar* tm = new GRBVar[N];
		GRBVar* cp = new GRBVar[N];

		for (int i = 0; i < N; i++) vars[i] = new GRBVar[N];

		GRBModel model = GRBModel(*env);

		//variables
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				vars[i][j] = model.addVar(0, 1, dis[i][j],
					GRB_BINARY, "x_" + itos(i) + "_" + itos(j));
			}
		}

		for (int i = 0; i < N; i++) {
			tm[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "t_" + itos(i));
			cp[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "c_" + itos(i));
		}
		//===========================================================================
		//constraints
		
		// Forbid edge from node back to itself
		tm[0].set(GRB_DoubleAttr_UB, 0);
		cp[0].set(GRB_DoubleAttr_UB, 0);
		for (int i = 0; i < N; i++) {
			vars[i][i].set(GRB_DoubleAttr_UB, 0);
		}

		//degree: one and balance
		for (int i = 0; i < N; i++) {
			GRBLinExpr expr1 = 0, expr2 = 0;
			for (int j = 0; j < N; j++) {
				expr1 += vars[i][j];
				expr2 += vars[j][i];
			}
			model.addConstr(expr1 == expr2, "balance_" + itos(i));
			if (i) model.addConstr(expr1 == 1, "degreeone_" + itos(i));
			else model.addConstr(expr1 >= 1, "degreeone_" + itos(i));
		}

		//connected graph
		for (int i = 0; i < N; i++) {
			for (int j = 1; j < N; j++) {
				model.addConstr(tm[j] >= tm[i] + service_time[i] + 
					dis[i][j] + fakeinf * (vars[i][j] - 1), "time_" + itos(i) + "_" + itos(j));
				model.addConstr(cp[j] >= cp[i] + demand[i] +
					fakeinf * (vars[i][j] - 1), "capacity_" + itos(i) + "_" + itos(j));
			}
		}

		//time window
		for (int i = 1; i < N; i++) {
			model.addConstr(tm[i] >= ready_time[i], "ready_time_" + itos(i));
			model.addConstr(tm[i] <= due_date[i], "due_date_" + itos(i));
		}

		//capacity
		for (int i = 1; i < N; i++) {
			model.addConstr(cp[i] + demand[i] <= capacity, "load_" + itos(i));
		}

		model.optimize();

		if (model.get(GRB_IntAttr_SolCount) > 0) {
			for (int i = 0; i < N; i++) {
				sol[i] = model.get(GRB_DoubleAttr_X, vars[i], N);
			}
		}

		delete[] cp;
		delete[] tm;
		for (int i = 0; i < N; i++) delete[] vars[i];
		delete[] vars;
		delete env;

	}
	catch (GRBException e) {
		cout << " Error number : " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << " Error during optimization " << endl;
	}

	gurobi_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;

	cout << endl << "Gurobi Time:\t" << gurobi_time << "s" << endl << endl;
}

void useCOPT() {
	cout << "use COPT here ..." << endl;

	clock_t start_time = clock();

	try {
		Envr env;
		Model model = env.CreateModel("VRPTW");

		//decision variables
		VarArray* vars = new VarArray[N];
		char* temp = nullptr;
		temp = new char[N];
		for (int i = 0; i < N; i++) temp[i] = COPT_BINARY;

		for (int i = 0; i < N; i++) {
			vars[i] = model.AddVars(N, nullptr, nullptr, dis[i], temp, "binaryvariables");
		}

		delete[] temp;

		VarArray tm = model.AddVars(N, COPT_CONTINUOUS, "time");
		VarArray cp = model.AddVars(N, COPT_CONTINUOUS, "capacity");

		//set objective
		{
			Expr expr(vars[0][0], 0);
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					expr += vars[i][j] * dis[i][j];
				}
			}
			//model.SetObjective(expr, COPT_MINIMIZE);
		}

		//===========================================================================
		//constraints
		// Forbid edge from node back to itself
		tm[0].Set(COPT_DBLINFO_UB, 0);
		cp[0].Set(COPT_DBLINFO_UB, 0);
		for (int i = 0; i < N; i++) {
			vars[i][i].Set(COPT_DBLINFO_UB, 0);
		}

		//degree: one and balance
		for (int i = 0; i < N; i++) {
			Expr expr1(vars[0][0],0), expr2(vars[0][0], 0);
			for (int j = 0; j < N; j++) {
				expr1 += vars[i][j];
				expr2 += vars[j][i];
			}
			model.AddConstr(expr1 == expr2, "balance");
			if (i) model.AddConstr(expr1 == 1, "degreeone");
			else model.AddConstr(expr1 >= 1, "degreeone");
		}

		//connected graph
		for (int i = 0; i < N; i++) {
			for (int j = 1; j < N; j++) {
				model.AddConstr(tm[j] >= tm[i] + service_time[i] +
					dis[i][j] + fakeinf * (vars[i][j] - 1), "time");
				model.AddConstr(cp[j] >= cp[i] + demand[i] +
					fakeinf * (vars[i][j] - 1), "capacity");
			}
		}

		//time window
		for (int i = 1; i < N; i++) {
			model.AddConstr(tm[i] >= ready_time[i], "ready_time");
			model.AddConstr(tm[i] <= due_date[i], "due_date");
		}

		//capacity
		for (int i = 1; i < N; i++) {
			model.AddConstr(cp[i] + demand[i] <= capacity, "load");
		}

		//===========================================================================

		model.Solve();

		if (model.GetIntAttr(COPT_INTATTR_HASLPSOL) > 0) {
			int k = 0;
			VarArray varsol = model.GetVars();
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					sol[i][j] = model.GetVars().GetVar(k++).Get(COPT_DBLINFO_VALUE);
				}
			}
		}
		
	}
	catch (CoptException e)
	{
		cout << "Error Code = " << e.GetCode() << endl;
		cout << e.what() << endl;
	}
	catch (...)
	{
		cout << "Unknown exception occurs!";
	}

	copt_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;

	cout << endl << "COPT Time:\t" << copt_time << "s" << endl << endl;
}

void drawSolution() {
	cout << "draw the solution ..." << endl;

	vector<double> xpr, ypr;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (sol[i][j] > 0.5) {
				xpr = { xcoord[i],xcoord[j] };
				ypr = { ycoord[i],ycoord[j] };
				plt::plot(xpr, ypr);
			}
		}
	}

	plt::axis("off");

	plt::title(filename + " N=" + itos(N - 1) + " COPT: " + dtos(copt_time) + "s Gurobi: " + dtos(gurobi_time) + "s");

	plt::show();
}

