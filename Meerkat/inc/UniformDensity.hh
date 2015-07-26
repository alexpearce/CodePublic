#ifndef UNIFORM_DENSITY
#define UNIFORM_DENSITY

#include "AbsDensity.hh"
#include "AbsPhaseSpace.hh"

#include <vector>

/// Class that describes the uniform density over any phase space

class UniformDensity : public AbsDensity {

  public: 
  
    //! Constructor 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space definition
    */ 
    UniformDensity(const char* pdfName, AbsPhaseSpace* thePhaseSpace);

    //! Destructor
    virtual ~UniformDensity(); 

    //! Calculate PDF value at the given point
    /*! 
        \param [in] x the point at which to calculate the PDF
        \return PDF value (equals 1 inside phase space, 0 outside it)
    */
    Double_t density(std::vector<Double_t> &x);

    //! Return phase space definition for this PDF
    /*!    
       \return PDF phase space
    */
    AbsPhaseSpace* phaseSpace() { return m_phaseSpace; }
    
  private: 
  
    //! Reference to phase space
    AbsPhaseSpace* m_phaseSpace; 
    
}; 

#endif
