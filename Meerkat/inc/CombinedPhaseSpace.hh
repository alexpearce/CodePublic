#ifndef COMBINED_PHASE_SPACE 
#define COMBINED_PHASE_SPACE

#include "AbsPhaseSpace.hh"

#include "TMath.h"

#include <vector>

/// Class that defines the phase space which is a combination (direct product) of several 
/// phase spaces. The dimensionality of the resulting phase space is a sum of dimensionalities 
/// of the component phase spaces. 

class CombinedPhaseSpace : public AbsPhaseSpace {

  public: 

    //! Constructor from the vector of component phase spaces
    /*! 
        \param [in] phaseSpaceName name of the phase space
        \param [in] phspList vector of component phase spaces. 
    */ 
    CombinedPhaseSpace(const char* phaseSpaceName, 
                       std::vector<AbsPhaseSpace*> &phspList); 

    //! Constructor for combined phase space of up to four components
    /*! 
        \param [in] phaseSpaceName name of the phase space
        \param [in] phaseSpace1 1st phase space
        \param [in] phaseSpace2 2nd phase space
        \param [in] phaseSpace3 3rd phase space
        \param [in] phaseSpace4 4th phase space
    */ 
    CombinedPhaseSpace(const char* phaseSpaceName, 
                       AbsPhaseSpace *phaseSpace1, 
                       AbsPhaseSpace *phaseSpace2, 
                       AbsPhaseSpace *phaseSpace3=0, 
                       AbsPhaseSpace *phaseSpace4=0); 

    //! Destructor
    virtual ~CombinedPhaseSpace(); 
    
    //! Get dimensionality of phase space
    /*! 
      \return dimensionality of the phase space
    */ 
    UInt_t dimensionality(); 

    //! Check if the point is within phase space limits
    /*! 
      \param [in] x point vector
      \return true if the point is within phase space limits, false otherwise
    */ 
    Bool_t withinLimits(std::vector<Double_t> &x); 

    //! Return lower allowed limit of the variable
    /*! 
      \param [in] var number of the variable 
      \return lower allowed limit
    */ 
    Double_t lowerLimit(UInt_t var); 

    //! Return upper allowed limit of the variable
    /*! 
      \param [in] var number of the variable 
      \return upper allowed limit
    */ 
    Double_t upperLimit(UInt_t var); 

    //! Return limits (lower and upper) for the variable at the certain point of the phase space
    /*! 
      \param [in] var number of the variable
      \param [in] x point in the phase space
      \param [out] lower lower limit
      \param [out] upper upper limit
    */ 
    Bool_t limits(UInt_t var, 
                  std::vector<Double_t> &x, 
                  Double_t* lower, 
                  Double_t* upper); 

  private: 

    //! Common initialise function for all constructors
    void init(std::vector<AbsPhaseSpace*> &phspList); 

    //! Cached dimensionality of the combined phase space
    UInt_t m_dim;
  
    //! vector of component phase spaces
    std::vector<AbsPhaseSpace*> m_phspList; 

    //! Vector of cached lower limits in each variable
    std::vector<Double_t> m_lowerLimit;
    
    //! Vector of cached upper limits in each variable
    std::vector<Double_t> m_upperLimit;


}; 

#endif
