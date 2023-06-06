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
	//path_constraints[0] = ;
	Tdouble Sp = state[0];
	Tdouble Lp = state[1];
	Tdouble Ip = state[2];
	Tdouble Sv = state[3];
	Tdouble Iv = state[4];
	Tdouble J = state[5];
	
	path_constraints[0] = Sp+Lp+Ip-1;
//	path_constraints[1] = -Sp;
//	path_constraints[2] = -Lp;
//	path_constraints[3] = -Ip;
//	path_constraints[4] = -Sv;
//	path_constraints[5] = -Iv;
//	path_constraints[6] = -J;

}


