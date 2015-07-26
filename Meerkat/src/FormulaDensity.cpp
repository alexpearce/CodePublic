#include <stdio.h>
#include <vector>

#include "TFormula.h"

#include "AbsPhaseSpace.hh"
#include "AbsDensity.hh"
#include "FormulaDensity.hh"

FormulaDensity::FormulaDensity(const char* pdfName, AbsPhaseSpace* thePhaseSpace, const char* formula) : AbsDensity(pdfName) {

  m_phaseSpace = thePhaseSpace;
  m_dim = m_phaseSpace->dimensionality(); 

  if (m_dim > 4) {
    printf("%20.20s ERROR: Dimensionality of the FormulaDensity phase space is %d (must not exceed 4)\n", 
           m_name, m_dim);
    abort(); 
  }

  m_formula = new TFormula(name(), formula); 
}

FormulaDensity::FormulaDensity(const char* pdfName, AbsPhaseSpace* thePhaseSpace, TFormula* formula) : AbsDensity(pdfName) {

  m_phaseSpace = thePhaseSpace;
  m_dim = m_phaseSpace->dimensionality(); 

  if (m_dim > 4) {
    printf("%20.20s ERROR: Dimensionality of the FormulaDensity phase space is %d (must not exceed 4)\n", 
           m_name, m_dim);
    abort(); 
  }

  m_formula = formula; 
}

FormulaDensity::~FormulaDensity() {

}

Double_t FormulaDensity::density(std::vector<Double_t> &x) {  

  UInt_t i;
  Double_t arg[4]; 
  for (i=0; i<4; i++) arg[i] = 0; 

  for (i=0; i<m_dim; i++) {
    arg[i] = x[i]; 
  }

  return m_formula->Eval(arg[0], arg[1], arg[2], arg[3]); 
}

