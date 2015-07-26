#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TMath.h"

#include "DalitzPhaseSpace.hh"

DalitzPhaseSpace::DalitzPhaseSpace(const char* phspName, 
                                   Double_t md, 
                                   Double_t ma, 
                                   Double_t mb, 
                                   Double_t mc) : AbsPhaseSpace(phspName) {
                                   
  m_a = ma;
  m_b = mb;
  m_c = mc;
  m_d = md;
  
  printf("%20.20s INFO: Creating Dalitz phase space for mother of mass %f and daughters of masses %f, %f, %f\n", 
         m_name, m_d, m_a, m_b, m_c); 
                                   
  m_a2 = m_a*m_a;
  m_b2 = m_b*m_b;
  m_c2 = m_c*m_c;
  m_d2 = m_d*m_d;
  m_SqSum = m_d2 + m_a2 + m_b2 + m_c2;
  
  m_MinAB = TMath::Power(m_a + m_b, 2);
  m_MaxAB = TMath::Power(m_d - m_c, 2);
  m_MinBC = TMath::Power(m_b + m_c, 2);
  m_MaxBC = TMath::Power(m_d - m_a, 2);
}

DalitzPhaseSpace::~DalitzPhaseSpace() {

}

Bool_t DalitzPhaseSpace::withinLimits(std::vector<Double_t> &x) {

  Double_t m2ab = x[0];
  Double_t m2bc = x[1];

//  printf("%f %f %f %f %f %f\n", m2ab, m2bc, m_MinAB, m_MinBC, m_MaxAB, m_MaxBC);

  if (m2ab < m_MinAB || m2ab > m_MaxAB || m2bc < m_MinBC || m2bc > m_MaxBC) return 0;

//  printf("1\n");

  Double_t eb = (m2ab - m_a2 + m_b2)/2./TMath::Sqrt(m2ab);
  Double_t ec = (m_d2 - m2ab - m_c2)/2./TMath::Sqrt(m2ab);

  Double_t p2b = TMath::Power(eb, 2) - m_b2;
  Double_t p2c = TMath::Power(ec, 2) - m_c2;

  if (p2c<0 || p2b<0) return 0;

//  printf("2\n");

  Double_t m2bc_max = TMath::Power(eb+ec, 2) - TMath::Power(sqrt(p2b) - sqrt(p2c), 2);
  Double_t m2bc_min = TMath::Power(eb+ec, 2) - TMath::Power(sqrt(p2b) + sqrt(p2c), 2);

  if (m2bc < m2bc_min || m2bc > m2bc_max) return 0;

//  printf("3\n");

  return 1;
}

Double_t DalitzPhaseSpace::lowerLimit(UInt_t var) {
  switch(var) {
    case 0: return m_MinAB; 
    case 1: return m_MinBC;
    default: 
      printf("%20.20s ERROR: var=%d for lower limit of 2D phase space\n", m_name, var);
      abort(); 
  }
}

Double_t DalitzPhaseSpace::upperLimit(UInt_t var) {
  switch(var) {
    case 0: return m_MaxAB; 
    case 1: return m_MaxBC;
    default: 
      printf("%20.20s ERROR: var=%d for upper limit of 2D phase space\n", m_name, var);
      abort(); 
  }
}

Bool_t DalitzPhaseSpace::limits(UInt_t var, std::vector<Double_t> &x, 
                                Double_t* lower, Double_t* upper) {

  if (var == 0) { // AB 
    Double_t m2bc = x[1];
    if (m2bc < m_MinBC || m2bc > m_MaxBC) return 0;

    Double_t eb = (m2bc + m_b2 - m_c2)/2./TMath::Sqrt(m2bc);
    Double_t ea = (m_d2 - m2bc - m_a2)/2./TMath::Sqrt(m2bc);

    Double_t p2b = TMath::Power(eb, 2) - m_b2;
    Double_t p2a = TMath::Power(ea, 2) - m_a2;

    if (p2a<0 || p2b<0) return 0;

//  printf("2\n");

    *upper = TMath::Power(ea+eb, 2) - TMath::Power(sqrt(p2a) - sqrt(p2b), 2);
    *lower = TMath::Power(ea+eb, 2) - TMath::Power(sqrt(p2a) + sqrt(p2b), 2);

  }
  else if (var == 1) { // BC

    Double_t m2ab = x[0];
    if (m2ab < m_MinAB || m2ab > m_MaxAB) return 0;

    Double_t eb = (m2ab - m_a2 + m_b2)/2./TMath::Sqrt(m2ab);
    Double_t ec = (m_d2 - m2ab - m_c2)/2./TMath::Sqrt(m2ab);

    Double_t p2b = TMath::Power(eb, 2) - m_b2;
    Double_t p2c = TMath::Power(ec, 2) - m_c2;

    if (p2c<0 || p2b<0) return 0;

//  printf("2\n");

    *upper = TMath::Power(eb+ec, 2) - TMath::Power(sqrt(p2b) - sqrt(p2c), 2);
    *lower = TMath::Power(eb+ec, 2) - TMath::Power(sqrt(p2b) + sqrt(p2c), 2);

  } else {
    printf("%20.20s ERROR: var=%d for limits of 2D phase space\n", m_name, var);
    abort();
  }
  
  return 1; 

}
