#ifndef POLYNOMIAL_DENSITY
#define POLYNOMIAL_DENSITY

#include "AbsDensity.hh"
#include "AbsPhaseSpace.hh"
#include "OneDimPhaseSpace.hh"

#include <vector>

/// Class that describes the polynomial density
/// which can be fitted to the data sample from NTuple
/// The phase space dimensionality must be either 1 or 2

class PolynomialDensity : public AbsDensity {

  public: 
  
    //! Constructor for 1D polynomial density. Integration over the 1D phase space is performed 
    //! analytically
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space definition
        \param [in] maxPower maximum power of the polynomial
        \param [in] tree ROOT tree containing data
        \param [in] var name of the tree variable
        \param [in] maxEvents maximum number of events to read from tree
    */ 
    PolynomialDensity(const char* pdfName, 
                      OneDimPhaseSpace* thePhaseSpace, 
                      UInt_t maxPower, 
                      TTree* tree, 
                      const char* var, 
                      UInt_t maxEvents = 0);

    //! Constructor for 2D polynomial density. Integration over the 2D phase space is 
    //! performed numerically (using MC integration). 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space definition
        \param [in] maxPower maximum power of the polynomial
        \param [in] tree ROOT tree containing data
        \param [in] var1 name of the 1st tree variable
        \param [in] var2 name of the 2nd tree variable
        \param [in] integEvents number of events for MC integration
        \param [in] maxEvents maximum number of events to read from tree
    */ 
    PolynomialDensity(const char* pdfName, 
                      AbsPhaseSpace* thePhaseSpace, 
                      UInt_t maxPower, 
                      TTree* tree, 
                      const char* var1, 
                      const char* var2, 
                      UInt_t integEvents, 
                      UInt_t maxEvents = 0);

    //! Destructor
    virtual ~PolynomialDensity(); 

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

    //! Read Ntuple contaning data points to be used for the fit. 
    void readTuple(TTree* tree, std::vector<TString> &vars, UInt_t maxEvents = 0);
  
    //! Reference to phase space
    AbsPhaseSpace* m_phaseSpace; 

    //! Cached dimensionality of the phase space
    UInt_t m_dim; 
    
    //! Vector of fitted parameters
    Double_t m_par[100]; 
    
    //! Power of the polynomial
    UInt_t m_power; 
}; 

#endif

