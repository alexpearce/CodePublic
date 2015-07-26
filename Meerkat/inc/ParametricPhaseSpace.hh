#ifndef PARAMETRIC_PHASE_SPACE
#define PARAMETRIC_PHASE_SPACE

#include "AbsPhaseSpace.hh"

#include "TMath.h"
#include "TFormula.h"

#include <vector>

/// Class which describes the parametric phase space
/// It adds one dimension on top of another phase space
/// with the upper and lower limits determined by TFormula 
/// expressions. 

class ParametricPhaseSpace : public AbsPhaseSpace {

  public: 
  
    //! Constructor
    /*! 
         \param [in] phspName phase space name
         \param [in] basePhaseSpace base phase space
         \param [in] lowerFormula formula expression for lower limit of additional variable
         \param [in] upperFormula formula expression for upper limit of additional variable
         \param [in] lowerLimit global lower limit of the additional variable
         \param [in] upperLimit global upper limit of the additional variable
    */ 
    ParametricPhaseSpace(const char* phspName, 
                         AbsPhaseSpace* basePhaseSpace, 
                         const char* lowerFormula, 
                         const char* upperFormula, 
                         Double_t lowerLimit, 
                         Double_t upperLimit); 

    //! Destructor
    virtual ~ParametricPhaseSpace(); 
    
    //! Get dimensionality of the phase space. 
    UInt_t dimensionality() { return m_dim; }

    //! Check if the point is within the phase space limits
    Bool_t withinLimits(std::vector<Double_t> &x); 

    //! Get lower limit
    Double_t lowerLimit(UInt_t var); 
    
    //! Get upper limit
    Double_t upperLimit(UInt_t var);

    //! Return limits (lower and upper) for the variable at the certain point of the phase space
    /*! 
      \param [in] var number of the variable
      \param [in] x point in the phase space
      \param [out] lower lower limit
      \param [out] upper upper limit
    */ 
    Bool_t limits(UInt_t var, std::vector<Double_t> &x, Double_t* lower, Double_t* upper); 

  private: 
  
    /// Lower limit in additional variable
    Double_t m_lowerLimit; 

    /// Upper limit in additional variable
    Double_t m_upperLimit; 

    /// Lower limit formula
    TFormula* m_lowerFormula; 
    
    /// Upper limit formula
    TFormula* m_upperFormula; 
    
    /// Cached dimensionality
    UInt_t m_dim; 

    /// Base phase space
    AbsPhaseSpace* m_phaseSpace; 
}; 

#endif

