#ifndef ABS_PHASE_SPACE 
#define ABS_PHASE_SPACE

#include "TMath.h"

#include <vector>

/// Abstract class which defines phase space interface. 

class AbsPhaseSpace {

  public: 
  
    //! Constructor
    /*! 
      \param [in] phaseSpaceName name of the phase space
    */ 
    AbsPhaseSpace(const char* phaseSpaceName); 

    //! Destructor
    virtual ~AbsPhaseSpace(); 

    //! Get dimensionality of phase space
    /*! 
      \return dimensionality of the phase space
    */ 
    virtual UInt_t dimensionality() = 0; 

    //! Check if the point is within phase space limits
    /*! 
      \param [in] x point vector
      \return true if the point is within phase space limits, false otherwise
    */ 
    virtual Bool_t withinLimits(std::vector<Double_t> &x) = 0; 

    //! Return lower allowed limit of the variable
    /*! 
      \param [in] var number of the variable 
      \return lower allowed limit
    */ 
    virtual Double_t lowerLimit(UInt_t var) = 0; 

    //! Return upper allowed limit of the variable
    /*! 
      \param [in] var number of the variable 
      \return upper allowed limit
    */ 
    virtual Double_t upperLimit(UInt_t var) = 0;

    //! Return limits (lower and upper) for the variable at a certain point of the phase space
    /*! 
      \param [in] var number of the variable
      \param [in] x point in the phase space
      \param [out] lower lower limit
      \param [out] upper upper limit
    */ 
    virtual Bool_t limits(UInt_t var, std::vector<Double_t> &x, Double_t* lower, Double_t* upper) = 0; 

    //! Return the name of the phase space
    /*! 
          \return Name of the phase space
    */
    const char* name(void) { return m_name; }

  protected: 
  
    //! Phase space name
    char m_name[256]; 

}; 

#endif
