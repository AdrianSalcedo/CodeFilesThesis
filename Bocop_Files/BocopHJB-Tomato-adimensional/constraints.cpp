// This code is published under the Eclipse Public License
// Authors: Daphne Giorgi, Benjamin Heymann, Jinyan Liu, Pierre Martinon, Olivier Tissot
// Inria Saclay and Cmap Ecole Polytechnique
// 2014-2017


// Function for the state admissibility
// Input :
// time : current time (t)
// state : vector of state variables (x)
// mode : current mode of the system (i)
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// true if the state is admissible
// false if it is not
#include "header_checkAdmissibleState"
{
  double Sp = state[0];
  double Lp = state[1];
  double Ip = state[2];
  double Sv = state[3];
  double Iv = state[4];
    if((0 <= Sp ) && (0 <= Lp )  && (0 <= Ip ) && (0 <= Sv)  && (0 <= Iv ))
  	  return true;
}


// Function for the (control,state) admissibility
// Input :
// time : current time (t)
// state : vector of state variables (x)
// control: vector of control variables (u)
// mode : current mode of the system (i)
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// true if the (control,state) pair is admissible
// false if it is not
#include "header_checkAdmissibleControlState"
{
  return true;
}


// Function for the final state admissibility
// Input :
// time : current time (t)
// final_state : vector of state variables (x)
// control: vector of control variables (u)
// mode : current mode of the system (i)
// constants : vector of constants
// dim_constant : dimension of the vector constants
// Output :
// true if the (control,state) pair is admissible
// false if it is not
#include "header_checkAdmissibleFinalState"
{
  return true;
}
