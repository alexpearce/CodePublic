#ifndef FACTORISED_DENSITY
#define FACTORISED_DENSITY

#include "AbsDensity.hh"
#include "AbsPhaseSpace.hh"

#include "TMath.h"

#include <vector>

/// Class that describes the factorised density: the product of densities in two or more component phase spaces.

class FactorisedDensity : public AbsDensity {

  public: 
  
    //! Constructor of factorised density of an arbitrary number of density components
    /*! 
      \param [in] pdfName PDF name
      \param [in] thePhaseSpace phase space. Dimensionality of the phase space should be equal to the sum of dimensionalities of all density components. 
      \param [in] densityComponents vector of density components. 
    */ 
    FactorisedDensity(const char* pdfName, 
                      AbsPhaseSpace* thePhaseSpace, 
                      std::vector<AbsDensity*> &densityComponents);

    //! Constructor of factorised density of up to four density components
    /*! 
      \param [in] pdfName PDF name
      \param [in] thePhaseSpace phase space. Dimensionality of the phase space should be equal to the sum of dimensionalities of all density components. 
      \param [in] d1 1st density component. 
      \param [in] d2 2nd density component. 
      \param [in] d3 3rd density component. 
      \param [in] d4 4th density component. 
    */ 
    FactorisedDensity(const char* pdfName, 
                      AbsPhaseSpace* thePhaseSpace, 
                      AbsDensity* d1, 
                      AbsDensity* d2, 
                      AbsDensity* d3 = 0, 
                      AbsDensity* d4 = 0);

    //! Destructor
    virtual ~FactorisedDensity(); 

    //! Calculate PDF value at the given point
    /*! 
        \param [in] x the point at which to calculate the PDF
        \return PDF value
    */
    Double_t density(std::vector<Double_t> &x);

    //! Return phase space definition for this PDF
    /*!
       \return PDF phase space
    */
    AbsPhaseSpace* phaseSpace() { return m_phaseSpace; }
    
  private: 
    
    //! Common initialise function used by all constructors. 
    /*! 
      \param [in] thePhaseSpace phase space
      \param [in] densityComponents vector of density components. 
    */ 
    void init(AbsPhaseSpace* thePhaseSpace, 
              std::vector<AbsDensity*> &densityComponents); 
  
    //! Reference to phase space
    AbsPhaseSpace* m_phaseSpace; 
    
    //! Vector of density components
    std::vector<AbsDensity*> m_densityComponents; 
    
    //! Cached dimensionality of the phase space
    UInt_t m_dim; 

}; 

#endif
