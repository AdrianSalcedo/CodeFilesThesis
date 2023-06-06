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
	double theta = constants[0];
	double mu = constants[1];
	double N = constants[2];
	double epsilon = constants[3];
	double delta_L = constants[4];
	double delta_V = constants[5];
	double delta_R = constants[6];
	double lambda_V = constants[7];
	double varepsilon = constants[8];
	double kappa = constants[9];
	double p = constants[10];
	double gamma_S = constants[11];
	double mu_I_S = constants[12];
	double delta_H = constants[13];
	double gamma_A = constants[14];
	double gamma_H = constants[15];
	double mu_H = constants[16];
	double a_S = constants[17];
	double a_H = constants[18];
	double a_D = constants[19];
	double c_L = constants[20];
	double c_V = constants[21];
	double beta_S = constants[22];
	double beta_A = constants[23];
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
	
	Tdouble u_L = control[0];
	Tdouble u_V = control[1];
	
	state_dynamics[0] = theta*mu*N-epsilon*(beta_A*I_A + beta_S*I_S)/N*L-(u_L + delta_L)*L -mu*L;
	state_dynamics[1] = (1-theta)*mu*N+(u_L+delta_L)*L + delta_V*V+delta_R*R-((beta_A*I_A + beta_S*I_S)/N + (lambda_V +u_V)+mu)*S;
	state_dynamics[2] = (beta_A*I_A + beta_S*I_S)/N*(epsilon*L+(1-varepsilon)*V+S)-(kappa+mu)*E;
	state_dynamics[3] = p*kappa*E-(gamma_S+mu_I_S+delta_H+mu)*I_S;
	state_dynamics[4] = (1-p)*kappa*E-(gamma_A+mu)*I_A;
	state_dynamics[5] = delta_H*I_S-(gamma_H +mu_H +mu)*H;
	state_dynamics[6] = gamma_S*I_S+gamma_A*I_A+gamma_H*H-(delta_R+mu)*R;
	state_dynamics[7] = mu_I_S*I_S+mu_H*H;
	state_dynamics[8] = (lambda_V+u_V)*S-((1-varepsilon)*(beta_A*I_A + beta_S*I_S)/N+delta_V+mu)*V;
	state_dynamics[9] = a_S*p*kappa*E+a_H*delta_H*I_S+a_D*(mu_I_S*I_S+mu_H*H)+0.5*c_L*u_L*u_L+0.5*c_V*u_V*u_V;
	state_dynamics[10] = (u_V + lambda_V)*(L+S+E+I_A+R);
	state_dynamics[11] = p*kappa*E;

}


