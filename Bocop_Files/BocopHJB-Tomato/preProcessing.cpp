// THIS FILE IS OPTIONAL !

// This code is published under the Eclipse Public License
// File: preProcessing.cpp
// Authors: Daphne Giorgi, Benjamin Heymann, Pierre Martinon, Olivier Tissot

// Function to make pre-processing, called before solving problem.
// Via this function, the user can access to the values of
// the definition and can change them.

// The following are the input and output available variables
// in post-processing function.
// Input :
// dim_constants : the number of constants
// constants     : the vector containing the constants

// Input/Output :
// starting_point: the vector for the starting point of the simulation
// starting_mode: the starting mode for the simulation

// Remember that the vectors numbering starts from 0
// (ex: the first component of the vector state is state[0])

#include "header_preProcessing"
{
    return 0;
}
