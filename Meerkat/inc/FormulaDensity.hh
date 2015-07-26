#ifndef FORMULA_DENSITY
#define FORMULA_DENSITY

#include "AbsDensity.hh"
#include "AbsPhaseSpace.hh"

#include "TFormula.h"

#include <vector>

/// Class that describes the density defined by ROOT formula
/// The phase space dimensionality must not exceed four

class FormulaDensity : public AbsDensity {

  public:

    //! Constructor 
    /*!
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space definition
        \param [in] formula ROOT formula string
    */
    FormulaDensity(const char* pdfName, 
                   AbsPhaseSpace* thePhaseSpace, 
                   const char* formula);

    //! Constructor
    /*!
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space definition
        \param [in] formula ROOT TFormula object
    */
    FormulaDensity(const char* pdfName, 
                   AbsPhaseSpace* thePhaseSpace, 
                   TFormula* formula);

    virtual ~FormulaDensity(); 

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

    //! ROOT formula
    TFormula* m_formula; 
    
    //! Cached dimensionality of the phase space
    UInt_t m_dim; 
}; 

#endif

