#ifndef ABS_DENSITY
#define ABS_DENSITY

#include "AbsPhaseSpace.hh"

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TRandom3.h"

#include <vector>

/// Abstract class which defines probability density interface

class AbsDensity {

  public: 

    //! Constructor. 
    /*! 
        \param [in] pdfName the PDF name
    */
    AbsDensity(const char* pdfName);
    
    //! Destructor
    virtual ~AbsDensity(); 
    
    //! Calculate PDF value at the given point
    /*! 
        \param [in] x the point at which to calculate the PDF
        \return PDF value
    */
    virtual Double_t density(std::vector<Double_t> &x) = 0;    

    //! Return phase space definition for this PDF
    /*!    
       \return PDF phase space
    */
    virtual AbsPhaseSpace* phaseSpace() = 0;

    //! Calculate 1D slice of the PDF
    /*!
       \param [in]  x    the point at which to calculate the slice
       \param [in]  num  number of the phase space variable to scan (value of x for this variable will be ignored)
       \param [out] hist output 1D ROOT histogram
       \param [in]  printout Enable printout
    */
    void slice(std::vector<Double_t> &x, UInt_t num, TH1F* hist, Bool_t printout = true);

    //! Calculate 2D slice of the PDF
    /*! 
       \param [in]  x    the point at which to calculate the slice
       \param [in]  numx number of the 1st phase space variable to scan (value of x for this variable will be ignored)
       \param [in]  numy number of the 2nd phase space variable to scan (value of x for this variable will be ignored)
       \param [out] hist output 2D ROOT histogram
       \param [in]  printout Enable printout
       \param [in]  inPhaseSpace If False, estimates the PDF somewhat outside of the phase space. Otherwise the bins outside of phase space are assigned zero.
    */
    void slice(std::vector<Double_t> &x, UInt_t numx, UInt_t numy, TH2F* hist, Bool_t printout = true, Bool_t inPhaseSpace = true);
    
    //! Calculate projection of the 1D PDF
    /*!
       \param [out] hist output 1D ROOT histogram
       \param [in]  printout Enable printout
    */
    void project(TH1F* hist, Bool_t printout = true);

    //! Calculate projection of the 2D PDF
    /*! 
       \param [out] hist output 2D ROOT histogram
       \param [in]  printout Enable printout
       \param [in]  inPhaseSpace If False, estimates the PDF somewhat outside of the phase space. Otherwise the bins outside of phase space are assigned zero.
    */
    void project(TH2F* hist, Bool_t printout = true, Bool_t inPhaseSpace = true);
    
    //! Set majorant for accept-reject method
    /*! 
       \param [in] majorant the majorant
    */ 
    void setMajorant(Double_t majorant) { m_majorant = majorant; }
    
    //! Set maximum number of tries for accept-reject method
    /*! 
       \param [in] maxTries maximum number of tries
    */ 
    void setMaxTries(UInt_t maxTries) { m_maxTries = maxTries; }
    
    //! Generate a single point within the phase space according to the PDF using accept-reject method. 
    /*! 
        \param [out] x the generated point
        \return the PDF value at this point. 
    */ 
    Double_t generate(std::vector<Double_t> &x); 
    
    //! Generate a sample of points within the phase space according to the PDF using accept-reject method. 
    /*! 
        Store them in the ROOT NTuple. 
          \param [out] tree ROOT NTuple to fill
          \param [in]  numEvents number of events to generate
    */ 
    void generate(TNtuple* tree, UInt_t numEvents); 
    
    //! Return the name of the PDF 
    /*! 
          \return PDF name
    */
    const char* name(void) { return m_name; }
    
  protected:

    //! PDF name
    char m_name[256];

    //! PDF majorant (maximum PDF value needed for accept-reject).
    Double_t m_majorant; 

    //! Maximum number of tries for accept-reject method
    UInt_t m_maxTries; 

    //! Random number generator
    TRandom3 m_rnd;
};

#endif
