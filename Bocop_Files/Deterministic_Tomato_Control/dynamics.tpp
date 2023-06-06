// Function for the dynamics of the problem
// dy/dt = dynamics(y,u,z,p)

// The following are the input and output available variables 
// for the dynamics of your optimal control problem.

// Input :
// time : current time (t)
// normalized_time: t renormalized in [0,1]
// initial_time : time value on the first discretization point
// final_time : time value on the last discretization point
// dim_* is the dimension of next vector in the declaration
// state : vector of state variables
// control : vector of control variables
// algebraicvars : vector of algebraic variables
// optimvars : vector of optimization parameters
// constants : vector of constants

// Output :
// state_dynamics : vector giving the expression of the dynamic of each state variable.

// The functions of your problem have to be written in C++ code
// Remember that the vectors numbering in C++ starts from 0
// (ex: the first component of the vector state is state[0])

// Tdouble variables correspond to values that can change during optimization:
// states, controls, algebraic variables and optimization parameters.
// Values that remain constant during optimization use standard types (double, int, ...).

#include "header_dynamics"
{
	// HERE : description of the function for the dynamics
	// Please give a function or a value for the dynamics of each state variable
	double betap = constants[0];
	double r1 = constants[1];
	double r2 = constants[2];
	double b = constants[3];
	double betav = constants[4];
	double gamma = constants[5];
	double gammaf = constants[6];
	double theta = constants[7];
	double mu = constants[8];
	double A1 = constants[9];
	double A2 = constants[10];
	double A3 = constants[11];
	double c1 = constants[12];
	double c2 = constants[13];
	double c3 = constants[14];
	double Np = constants[15];
	double Nv = constants[16];

	Nv = mu /(gamma + gammaf);
	
	Tdouble Sp = state[0];
	Tdouble Lp = state[1];
	Tdouble Ip = state[2];
	Tdouble Sv = state[3];
	Tdouble Iv = state[4];
	Tdouble J = state[5];

	Tdouble u1 = control[0];
	Tdouble u2 = control[1];
	Tdouble u3 = control[2];

	state_dynamics[0] = -betap * Sp * Iv + (r1 + u1) * Lp + (r2 + u2) * Ip;
	state_dynamics[1] = betap * Sp * Iv - (b + r1 + u1) * Lp;
	state_dynamics[2] = b * Lp - (r2 + u2) * Ip;
	state_dynamics[3] = -betav * Sv * Ip -(gamma + gammaf + u3)*Sv + (1-theta)*mu;//+ ((1 - theta - Sv)* mu) / Nv;
	state_dynamics[4] = betav * Sv * Ip -(gamma + gammaf + u3)*Iv + theta*mu;//+ ((theta - Iv) * mu) / Nv;
	state_dynamics[5] = A1 * Lp + A2 * Ip + A3 * Iv + 0.5 * c1 * u1 * u1 + 0.5 * c2 * u2 * u2 + 0.5 * c3 * u3 * u3;
}


