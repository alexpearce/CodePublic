#ifndef ONEDIM_PHASE_SPACE
#define ONEDIM_PHASE_SPACE

#include "AbsPhaseSpace.hh"

#include "TMath.h"

#include <vector>

/// Class which describes one-dimensional flat phase space
/// (i.e. the finite range in one variable)

class OneDimPhaseSpace : public AbsPhaseSpace {

  public: 
  
    //! Constructor
    /*!  
         \param [in] phaseSpaceName name of the phase space
         \param [in] lowerLimit lower limit of the phase space
         \param [in] upperLimit upper limit of the phase space
    */ 
    OneDimPhaseSpace(const char* phaseSpaceName, 
                     Double_t lowerLimit, 
                     Double_t upperLimit); 

    //! Destructor
    virtual ~OneDimPhaseSpace(); 
    
    //! Get dimensionality of the phase space. Always equals to one. 
    UInt_t dimensionality() { return 1; }

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
  
    /// Lower limit
    Double_t m_lowerLimit; 

    /// Upper limit    
    Double_t m_upperLimit; 

}; 

#endif
