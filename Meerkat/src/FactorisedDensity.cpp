#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TMath.h"

#include "AbsPhaseSpace.hh"
#include "AbsDensity.hh"
#include "FactorisedDensity.hh"

FactorisedDensity::FactorisedDensity(const char* pdfName, 
                                     AbsPhaseSpace* thePhaseSpace, 
                                     std::vector<AbsDensity*> &densityComponents) : AbsDensity(pdfName) {

  init(thePhaseSpace, densityComponents); 
}

FactorisedDensity::FactorisedDensity(const char* pdfName, 
                         AbsPhaseSpace* thePhaseSpace,
                         AbsDensity* d1, AbsDensity* d2, 
                         AbsDensity* d3, AbsDensity* d4) : AbsDensity(pdfName) {

  std::vector<AbsDensity*> densityVector(2); 
  densityVector[0] = d1; 
  densityVector[1] = d2; 
  if (d3) {
    densityVector.push_back(d3); 
  }
  if (d4) {
    densityVector.push_back(d4); 
  }
  init(thePhaseSpace, densityVector); 
}

FactorisedDensity::~FactorisedDensity() {

}

void FactorisedDensity::init(AbsPhaseSpace* thePhaseSpace, 
                std::vector<AbsDensity*> &densityComponents) {

  m_phaseSpace = thePhaseSpace;
  m_densityComponents = densityComponents; 
  
  printf("%20.20s INFO: Creating factorised density of %d components\n", m_name, (UInt_t)densityComponents.size() ); 

  UInt_t sumDim = 0; 
  std::vector<AbsDensity*>::iterator i; 
  for (i=densityComponents.begin(); i != densityComponents.end(); i++) {
    sumDim += (*i)->phaseSpace()->dimensionality(); 
  }

  m_dim = m_phaseSpace->dimensionality(); 
  if (sumDim != m_dim) {
    printf("%20.20s ERROR: Dimensionalities of the factorised phase space and sum of its components differ\n", m_name);
    abort(); 
  }

}

Double_t FactorisedDensity::density(std::vector<Double_t> &x) {  

  if (x.size() != m_dim) {
    printf("%20.20s ERROR: Dimensionality of vector (%d) does not correspond to phase space definition (%d)\n", 
       m_name, (UInt_t)x.size(), m_dim);
    abort(); 
  }
  
  Double_t prod = 1.;
  UInt_t n=0;
  std::vector<AbsDensity*>::iterator i; 
  for (i=m_densityComponents.begin(); i != m_densityComponents.end(); i++) {
    UInt_t dim = (*i)->phaseSpace()->dimensionality(); 
    std::vector<Double_t> subset(dim); 
    UInt_t j; 
    for (j=0; j<dim; j++) {
      subset[j] = x[n];
      n++;
    }
    Double_t d = (*i)->density(subset); 
    
    prod *= d;
    
//    printf("Component %d, density=%f, prod=%f\n", n, d, prod);
  }
  return prod; 

}
