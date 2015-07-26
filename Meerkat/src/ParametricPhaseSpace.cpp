#include "ParametricPhaseSpace.hh"
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "TMath.h"
#include "TFormula.h"

ParametricPhaseSpace::ParametricPhaseSpace(const char* phspName, 
                                           AbsPhaseSpace* phsp, 
                                           const char* lowerFormula, 
                                           const char* upperFormula, 
                                           Double_t lower, 
                                           Double_t upper) : AbsPhaseSpace (phspName) {

  char upperName[32]; 
  char lowerName[32]; 

  printf("%20.20s INFO: Creating parametric phase space on top of base phase space \"%s\"\n", m_name, phsp->name()); 

  sprintf(upperName, "%20.20s_Upper", phspName); 
  sprintf(lowerName, "%20.20s_Lower", phspName); 

  m_lowerFormula = new TFormula(lowerName, lowerFormula);
  m_upperFormula = new TFormula(upperName, upperFormula);
  
  m_lowerLimit = lower; 
  m_upperLimit = upper; 
  
  m_phaseSpace = phsp; 
  
  if (m_phaseSpace->dimensionality() > 4) {
    printf("%20.20s ERROR: Dimensionality of the base phase space \"%s\" for parametric phase space exceeds 4\n", 
           m_name, m_phaseSpace->name() ); 
    abort(); 
  }
  m_dim = m_phaseSpace->dimensionality() + 1; 
  
  if (m_upperLimit <= m_lowerLimit) {
    printf("%20.20s ERROR: Upper limit (%f) <= lowerLimit (%f) for parametric phase space\n", 
           m_name, m_upperLimit, m_lowerLimit);
    abort(); 
  }

}

ParametricPhaseSpace::~ParametricPhaseSpace() {

}

Bool_t ParametricPhaseSpace::withinLimits(std::vector<Double_t> &x) {
  Double_t x1 = x[m_dim-1];
  if (x1 > m_upperLimit || x1 < m_lowerLimit) return 0;
  if (m_phaseSpace->withinLimits(x) == 0) return 0; 
  
  Double_t arg[4]; 
  UInt_t i; 
  for (i=0; i<4; i++) arg[i] = 0; 
  for (i=0; i<m_dim-1; i++) {
    arg[i] = x[i]; 
  }

  Double_t upper = m_upperFormula->Eval(arg[0], arg[1], arg[2], arg[3]); 
  Double_t lower = m_lowerFormula->Eval(arg[0], arg[1], arg[2], arg[3]); 
  if (upper < lower) {
    printf("%20.20s WARNING: Upper limit (%f) < lower limit (%f)\n for parametric phase space\n", m_name, upper, lower); 
  }

  if (x1 > upper || x1 < lower) return 0;
  return 1;
}

Double_t ParametricPhaseSpace::lowerLimit(UInt_t var) {
  if (var == m_dim-1) return m_lowerLimit;
  else if (var < m_dim-1) return m_phaseSpace->lowerLimit(var); 
  printf("%20.20s ERROR: var=%d for lower limit of parametric phase space\n", m_name, var);
  abort(); 
  return 0;
}

Double_t ParametricPhaseSpace::upperLimit(UInt_t var) {
  if (var == m_dim-1) return m_upperLimit;
  else if (var < m_dim-1) return m_phaseSpace->upperLimit(var); 
  printf("%20.20s ERROR: var=%d for upper limit of parametric phase space\n", m_name, var);
  abort(); 
  return 0;

}

Bool_t ParametricPhaseSpace::limits(UInt_t var, std::vector<Double_t> &x, 
                              Double_t* lower, Double_t* upper) {
  if (var == m_dim-1) {
    if (m_phaseSpace->withinLimits(x)) {
      Double_t arg[4]; 
      UInt_t i; 
      for (i=0; i<4; i++) arg[i] = 0; 
      for (i=0; i<m_dim-1; i++) {
        arg[i] = x[i]; 
      }
      *upper = m_upperFormula->Eval(arg[0], arg[1], arg[2], arg[3]); 
      *lower = m_lowerFormula->Eval(arg[0], arg[1], arg[2], arg[3]); 
      return 1; 
    }
    return 0; 
  }
  else if (var < m_dim-1) {
    return limits(var, x, lower, upper);
  }
  printf("%20.20s ERROR: var=%d for limits of parametric phase space\n", m_name, var);
  abort(); 
  return 0;
}

