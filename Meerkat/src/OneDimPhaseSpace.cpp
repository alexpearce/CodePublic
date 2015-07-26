#include "OneDimPhaseSpace.hh"
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "TMath.h"

// Define some tolerance to be used to make sure that the lower and upper limits 
// fall into the phase space
#define TOLERANCE 1e-30

OneDimPhaseSpace::OneDimPhaseSpace(const char* phspName, 
                                   Double_t lower, 
                                   Double_t upper) : AbsPhaseSpace(phspName) {

  printf("%20.20s INFO: Creating 1D phase space in range (%f, %f)\n", m_name, lower, upper); 

  m_lowerLimit = lower;
  m_upperLimit = upper;
  
  if (m_upperLimit <= m_lowerLimit) {
    printf("%20.20s ERROR: Upper limit (%f) <= lowerLimit (%f) for flat phase space\n", 
           m_name, m_upperLimit, m_lowerLimit);
    abort(); 
  }

}

OneDimPhaseSpace::~OneDimPhaseSpace() {

}

Bool_t OneDimPhaseSpace::withinLimits(std::vector<Double_t> &x) {
  Double_t x1 = x[0];
//  if (x1 > m_upperLimit + TOLERANCE || x1 < m_lowerLimit - TOLERANCE) return 0;
  if (x1 > m_upperLimit || x1 < m_lowerLimit) return 0;
  return 1;
}

Double_t OneDimPhaseSpace::lowerLimit(UInt_t var) {
  if (var == 0) return m_lowerLimit;
  printf("%20.20s ERROR: var=%d for lower limit of 1D phase space\n", m_name, var);
  abort(); 
  return 0;
}

Double_t OneDimPhaseSpace::upperLimit(UInt_t var) {
  if (var == 0) return m_upperLimit;
  printf("%20.20s ERROR: var=%d for upper limit of 1D phase space\n", m_name, var);
  abort(); 
  return 0;

}

Bool_t OneDimPhaseSpace::limits(UInt_t var, __attribute__((unused)) std::vector<Double_t> &x, 
                              Double_t* lower, Double_t* upper) {

  if (var == 0) {
    *lower = m_lowerLimit;
    *upper = m_upperLimit;
    return 1;
  }
  printf("%20.20s ERROR: var=%d for limits of 1D phase space\n", m_name, var);
  abort(); 
  return 0;
}
