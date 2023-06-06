// This code is published under the Eclipse Public License
// Authors: Daphne Giorgi, Benjamin Heymann, Jinyan Liu, Pierre Martinon, Olivier Tissot
// Inria Saclay and Cmap Ecole Polytechnique
// 2014-2017


// Function for the running cost
// Input :
// time : current time t
// initial_time : t0
// final_time : tf
// state : vector of state variables x
// control : vector of control variables u
// mode : mode of the system i
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// running_cost : running cost l(t,x,u,i)
#include "header_runningCost"
{
  double Sp = state[0];
  double Lp = state[1];
  double Ip = state[2];
  double Sv = state[3];
  double Iv = state[4];
  double u1 = control[0];
  double u2 = control[1];
  double u3 = control[2];
  running_cost = 0.4*Lp+0.3*Ip+0.2*Iv+0.5*(u1*u1 + u2*u2 + u3*u3);
}


// Function for the final cost
// Input :
// initial_time : t0
// final_time : tf
// final_state : vector of state variables x_f
// final_mode : final mode of the system i_f
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// final_cost : final cost g(t0,tf,x_f,i_f)
#include "header_finalCost"
{
  final_cost = 0e0;
}


// Function for the switching cost
// Input :
// current_mode : current mode
// next_mode : next mode
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// switching_cost : switching cost s(i_k,i_k+1)
#include "header_switchingCost"
{
  switching_cost = 0e0;
}

