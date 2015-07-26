#ifndef BINNED_DENSITY
#define BINNED_DENSITY

#include "AbsDensity.hh"

#include "TMath.h"

#include <vector>

/// A class that describes the PDF which is a multilinear interpolation 
/// on the binned grid. 

class BinnedDensity : public AbsDensity {

  public: 

    //! Constructor that creates the binned density from any AbsDensity of the dimensionality up to five. 
    /*!  
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] d input density
        \param [in] bins1 number of bins in 1st variable
        \param [in] bins2 number of bins in 2nd variable
        \param [in] bins3 number of bins in 3rd variable
        \param [in] bins4 number of bins in 4th variable
        \param [in] bins5 number of bins in 5th variable
    */ 
    BinnedDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  AbsDensity* d, 
                  UInt_t bins1, 
                  UInt_t bins2 = 0, 
                  UInt_t bins3 = 0, 
                  UInt_t bins4 = 0, 
                  UInt_t bins5 = 0);

    //! Constructor that creates the binned density from any AbsDensity of arbitrary dimensionality. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] binning vector of bin numbers for each variable. Vector size should match the dimensionality of the phase space. 
        \param [in] d input density
    */ 
    BinnedDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  std::vector<UInt_t> &binning, 
                  AbsDensity* d);

    //! Constructor that reads the binned density from a file. The dimensionality of the density stored in the file should match the dimensionality of the phase space. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] fileName input file name
    */ 
    BinnedDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  const char* fileName); 

    //! Destructor
    virtual ~BinnedDensity(); 

    //! Read the binned density from a file, ROOT or text depending on the extension. The dimensionality of the density stored in the file should match the dimensionality of the phase space. 
    /*! 
        \param [in] fileName input file name
    */ 
    void readFromFile(const char* fileName); 

    //! Read the binned density from a text file. The dimensionality of the density stored in the file should match the dimensionality of the phase space. 
    /*! 
        \param [in] fileName input file name
    */ 
    void readFromTextFile(const char* fileName); 

    //! Read the binned density from a ROOT file. The dimensionality of the density stored in the file should match the dimensionality of the phase space. 
    /*! 
        \param [in] fileName input file name
    */ 
    void readFromRootFile(const char* fileName); 

    //! Write the binned density into a file (ROOT or text format depending on the extension)
    /*!
        \param [in] fileName output file name
    */ 
    void writeToFile(const char* fileName); 

    //! Write the binned density into a file
    /*!
        \param [in] fileName output file name
    */ 
    void writeToTextFile(const char* fileName); 

    //! Write the binned density into a ROOT file
    /*!
        \param [in] fileName output ROOT file name
    */ 
    void writeToRootFile(const char* fileName); 

    //! Return the value of the PDF in a point
    /*! 
        \param [in] x point
        \return PDF value
    */ 
    Double_t density(std::vector<Double_t> &x);

    //! Return the phase space
    /*! 
        \return phase space
    */ 
    AbsPhaseSpace* phaseSpace() { return m_phaseSpace; }

  private: 

    //! Common initalisation function used by both constructors from an AbsDensity. 
    /*! 
        \param [in] thePhaseSpace phase space
        \param [in] binning vector of bin numbers for each variable. Vector size should match the dimensionality of the phase space. 
        \param [in] d input density
    */ 
    void init(AbsPhaseSpace* thePhaseSpace, 
              std::vector<UInt_t> &binning, 
              AbsDensity* d);

    //! Map of PDF values in bins
    std::vector<Double_t> m_map;

    //! Vector of bin numbers for each variable
    std::vector<UInt_t> m_binning;

    //! Reference to the phase space definition
    AbsPhaseSpace* m_phaseSpace; 

    //! Reference to the input density
    AbsDensity* m_density; 

};

#endif
