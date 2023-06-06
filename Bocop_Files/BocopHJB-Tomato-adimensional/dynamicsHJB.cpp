// This code is published under the Eclipse Public License
// Authors: Daphne Giorgi, Benjamin Heymann, Jinyan Liu, Pierre Martinon, Olivier Tissot
// Inria Saclay and Cmap Ecole Polytechnique
// 2014-2017

// General dynamics
// dy/dt = drift(t,y,u)dt + volatility(y,u)dWt where Wt is the standard Brownian motion

// Function for the drift (deterministic dynamics)
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
// state_dynamics : drift f(t,x,u,i) ie deterministic dynamics
#include "header_drift"
{

  //Example: double integrator, state = (x1 x2), control = (u) 
  //\dot x1 = x2
  //\dot x2 = u

  double Sp = state[0];
  double Lp = state[1];
  double Ip = state[2];
  double Sv = state[3];
  double Iv = state[4];
  //double Bp = state[5];
  //double Bv = state[6];
  double u1 = control[0];
  double u2 = control[1];
  double u3 = control[2];
  double betap = constants[0];
  double r1 = constants[1];
  double r2 = constants[2];
  double b = constants[3];
  double betav = constants[4];
  double gamma = constants[5];
  double gammaf = constants[6];
  double theta = constants[7];
  double mu = constants[8];
  double sigmaL = constants[9];
  double sigmaI = constants[10];
  double sigmaV = constants[11];
  double Nv = Sv+Iv;
  state_dynamics[0] = -betap*Sp*Iv +(r1+u1)*Lp +(r2+u2)*Ip;
  state_dynamics[1] = betap*Sp*Iv -(b+r1+u1)*Lp;
  state_dynamics[2] = b*Lp-(r2+u2)*Ip;
  state_dynamics[3] = -betav*Sv*Ip+(1-theta-Sv-u3)*mu/Nv;
  state_dynamics[4] = betav*Sv*Ip +(theta-Iv-u3)*mu/Nv;
  //state_dynamics[5] = 0;
  //state_dynamics[6] = 0;
}


// Function for the volatility (stochastic dynamics)
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
// volatility_dynamics : vector giving the volatility expression of the volatility
// Remember that this is a matrix of dimension dim_state x dim_brownian and you have to fill every coefficient.
#include "header_volatility"
{
  double sigmaL = constants[9];
  double sigmaI = constants[10];
  double sigmaV = constants[11];
  double Sp = state[0];
  double Lp = state[1];
  double Ip = state[2];
  double Sv = state[3];
  double Iv = state[4];
  double Bp = state[5];
  double Bv = state[6];
  volatility_dynamics[0][0] = Sp*(sigmaL*Lp+sigmaI*Ip);
  volatility_dynamics[0][1] = 0;
  volatility_dynamics[1][0] = -sigmaL*Sp*Lp;
  volatility_dynamics[1][1] = 0;
  volatility_dynamics[2][0] = -sigmaI*Sp*Ip;
  volatility_dynamics[2][1] = 0;
  volatility_dynamics[3][0] = 0;  
  volatility_dynamics[3][1] = 0;//-sigmaV*Sv;
  volatility_dynamics[4][0] = 0;
  volatility_dynamics[4][1] = 0;//-sigmaV*Iv;
  // volatility_dynamics[5][0] = 1;
  // volatility_dynamics[5][1] = 0;
  // volatility_dynamics[6][0] = 0;
  // volatility_dynamics[6][1] = 1;
}
