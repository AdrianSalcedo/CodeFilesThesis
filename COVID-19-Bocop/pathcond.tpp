// Function for the path constraints of the problem
// a <= g(t,y,u,z,p) <= b

// The following are the input and output available variables 
// for the path constraints of your optimal control problem.

// Input :
// dim_path_constraints : number of path constraints
// time : current time (t)
// initial_time : time value on the first discretization point
// final_time : time value on the last discretization point
// dim_* is the dimension of next vector in the declaration
// state : vector of state variables
// control : vector of control variables
// algebraicvars : vector of algebraic variables
// optimvars : vector of optimization vector of optimization parameters
// constants : vector of constants

// Output :
// path_constraints : vector of path constraints expressions ("g" in the example above)

// The functions of your problem have to be written in C++ code
// Remember that the vectors numbering in C++ starts from 0
// (ex: the first component of the vector state is state[0])

// Tdouble variables correspond to values that can change during optimization:
// states, controls, algebraic variables and optimization parameters.
// Values that remain constant during optimization use standard types (double, int, ...).

#include "header_pathcond"
{
	// HERE : description of the path constraints
	// Please give a function or a value for each path constraint
	double kappa = constants[9];
	double B = constants[24];
	Tdouble L = state[0];
	Tdouble S = state[1];
	Tdouble E = state[2];
	Tdouble I_S = state[3];
	Tdouble I_A = state[4];
	Tdouble H = state[5];
	Tdouble R = state[6];
	Tdouble D = state[7];
	Tdouble V = state[8];
	Tdouble J = state[9];

	path_constraints[0] = L+S+E+I_S+I_A+H+R+V-1;
	path_constraints[1] = -L;
	path_constraints[2] = -S;
	path_constraints[3] = -E;
	path_constraints[4] = -I_S;
	path_constraints[5] = -I_A;
	path_constraints[1] = -H;
	path_constraints[7] = -R;
	path_constraints[8] = -D;
	path_constraints[9] = -V;
	path_constraints[10] = -J;
	path_constraints[11] = kappa*I_S - B;
}


