#include <stdio.h>
#include <vector>

#include "TMath.h"

#include "CombinedPhaseSpace.hh"

CombinedPhaseSpace::CombinedPhaseSpace(const char* phspName, 
                                       std::vector<AbsPhaseSpace*> &phspList) : AbsPhaseSpace(phspName) {
  init(phspList); 
}

CombinedPhaseSpace::CombinedPhaseSpace(const char* phspName, 
                                       AbsPhaseSpace *phaseSpace1, 
                                       AbsPhaseSpace *phaseSpace2, 
                                       AbsPhaseSpace *phaseSpace3, 
                                       AbsPhaseSpace *phaseSpace4) : AbsPhaseSpace(phspName) {

  std::vector<AbsPhaseSpace*> phspVector(2); 
  phspVector[0] = phaseSpace1; 
  phspVector[1] = phaseSpace2; 
  if (phaseSpace3) {
    phspVector.push_back(phaseSpace3); 
  }
  if (phaseSpace4) {
    phspVector.push_back(phaseSpace4); 
  }

  init(phspVector); 

}

CombinedPhaseSpace::~CombinedPhaseSpace() {

}

void CombinedPhaseSpace::init(std::vector<AbsPhaseSpace*> &phspList) {

  m_phspList = phspList; 
  m_dim = 0;
  
  std::vector<AbsPhaseSpace*>::iterator i;
  
  printf("%20.20s INFO: Creating combined phase space of %d components\n", m_name, (UInt_t)m_phspList.size() ); 

  // Calculate dimensionality as a sum of component phase spaces  
  for (i=m_phspList.begin(); i != m_phspList.end(); i++) {
    UInt_t dim = (*i)->dimensionality(); 
    printf("%20.20s INFO: Adding component phase space \"%s\", dim=%d\n", m_name, (*i)->name(), dim); 
    m_dim += dim;
  }
  printf("%20.20s INFO: Resulting dimensionality is %d\n", m_name, m_dim); 

  m_lowerLimit.resize(m_dim);
  m_upperLimit.resize(m_dim);

  // Cache lower and upper limits for variables from component phase spaces    
  UInt_t n = 0;
  for (i=m_phspList.begin(); i != m_phspList.end(); i++) {
    UInt_t dim = (*i)->dimensionality(); 
    UInt_t j;
    for (j=0; j<dim; j++) {
      m_lowerLimit[n] = (*i)->lowerLimit(j);
      m_upperLimit[n] = (*i)->upperLimit(j);
      n++;
    }
  }
}


UInt_t CombinedPhaseSpace::dimensionality() {
  return m_dim;
}

Bool_t CombinedPhaseSpace::withinLimits(std::vector<Double_t> &x) {
  if (x.size() != m_dim) {
    printf("%20.20s WARNING: Dimensionality of vector (%d) does not correspond to phase space definition (%d)\n", 
       m_name, (UInt_t)x.size(), m_dim);
    return 0;
  }
  
  UInt_t n=0;
  std::vector<AbsPhaseSpace*>::iterator i; 
  for (i=m_phspList.begin(); i != m_phspList.end(); i++) {
    UInt_t dim = (*i)->dimensionality(); 
    UInt_t j; 
    std::vector<Double_t> subset(dim); 
    for (j=0; j<dim; j++) {
      subset[j] = x[n];
      n++;
    }
    if ( !(*i)->withinLimits(subset) ) return 0; 
  }
  
  return 1;
}

Double_t CombinedPhaseSpace::lowerLimit(UInt_t var) {
  if (var < m_lowerLimit.size()) {
    return m_lowerLimit[var];
  }

  printf("%20.20s ERROR: var=%d for lower limit of combined phase space of dim=%d\n", m_name, var, m_dim);
  abort(); 
  
  return 0;
}

Double_t CombinedPhaseSpace::upperLimit(UInt_t var) {
  if (var < m_upperLimit.size()) {
    return m_upperLimit[var];
  }

  printf("%20.20s ERROR: var=%d for upper limit of combined phase space of dim=%d\n", m_name, var, m_dim);
  abort(); 

  return 0;
}

Bool_t CombinedPhaseSpace::limits(UInt_t var, std::vector<Double_t> &x, 
                                  Double_t* lower, Double_t* upper) {

  UInt_t n = 0;
  UInt_t m = 0;
  std::vector<AbsPhaseSpace*>::iterator i; 

  for (i=m_phspList.begin(); i != m_phspList.end(); i++) {
    UInt_t dim = (*i)->dimensionality(); 
    std::vector<Double_t> subset(dim); 
    UInt_t j;

    for (j=0; j<dim; j++) {
      subset[j] = x[m];
      m++;
    }

    for (j=0; j<dim; j++) {
      if (x[n] < m_lowerLimit[n] || x[n] > m_upperLimit[n]) 
        return 0;
      if (n == var) {
        Bool_t res = (*i)->limits(j, subset, lower, upper);
        if (!res) return 0;
      }
      n++;
    }
  }
  return 1;

}
