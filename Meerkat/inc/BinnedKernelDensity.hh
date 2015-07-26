#ifndef BINNED_KERNEL_DENSITY
#define BINNED_KERNEL_DENSITY

#include "AbsDensity.hh"

#include "TMath.h"

#include <vector>

class BinnedKernelDensity : public AbsDensity {

  public: 

    //! Constructor for 1-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 variable name
        \param [in] bins1 number of bins for the binned interpolation
        \param [in] width1 kernel width
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, 
                  UInt_t bins1, 
                  Double_t width1, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0, 
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 1-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple with weight. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 variable name
        \param [in] weight weight variable name
        \param [in] bins1 number of bins for the binned interpolation
        \param [in] width1 kernel width
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, 
                  const char *weight, 
                  UInt_t bins1, 
                  Double_t width1, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0, 
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 2-dimensional kernel PDF with binned interpolation from the sample of points in an NTUple. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, 
                  UInt_t bins1, UInt_t bins2, 
                  Double_t width1, Double_t width2, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 2-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple with weight. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] weight weight variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, 
                  const char* weight, 
                  UInt_t bins1, UInt_t bins2, 
                  Double_t width1, Double_t width2, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 3-dimensional kernel PDF with binned interpolation from the sample of points in an NTUple. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] vars3 3rd variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] bins3 number of bins for the 3rd variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] width3 kernel width for the 3rd variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 

    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, const char *vars3, 
                  UInt_t bins1, UInt_t bins2, UInt_t bins3, 
                  Double_t width1, Double_t width2, Double_t width3, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 3-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple with weight. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] vars3 3rd variable name
        \param [in] weight weight variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] bins3 number of bins for the 3rd variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] width3 kernel width for the 3rd variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 

    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, const char *vars3, 
                  const char *weight, 
                  UInt_t bins1, UInt_t bins2, UInt_t bins3, 
                  Double_t width1, Double_t width2, Double_t width3, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 4-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] vars3 3rd variable name
        \param [in] vars4 4th variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] bins3 number of bins for the 3rd variable for the binned interpolation
        \param [in] bins4 number of bins for the 4th variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] width3 kernel width for the 3rd variable
        \param [in] width4 kernel width for the 4th variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, const char *vars3, const char *vars4, 
                  UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, 
                  Double_t width1, Double_t width2, Double_t width3, Double_t width4, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 4-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple with weight. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] vars3 3rd variable name
        \param [in] vars4 4th variable name
        \param [in] weight weight variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] bins3 number of bins for the 3rd variable for the binned interpolation
        \param [in] bins4 number of bins for the 4th variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] width3 kernel width for the 3rd variable
        \param [in] width4 kernel width for the 4th variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, const char *vars3, const char *vars4, 
                  const char *weight, 
                  UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, 
                  Double_t width1, Double_t width2, Double_t width3, Double_t width4, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 5-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] vars3 3rd variable name
        \param [in] vars4 4th variable name
        \param [in] vars5 5th variable name
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] bins3 number of bins for the 3rd variable for the binned interpolation
        \param [in] bins4 number of bins for the 4th variable for the binned interpolation
        \param [in] bins5 number of bins for the 5th variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] width3 kernel width for the 3rd variable
        \param [in] width4 kernel width for the 4th variable
        \param [in] width5 kernel width for the 5th variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, const char *vars3, const char *vars4, const char *vars5, 
                  UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, UInt_t bins5, 
                  Double_t width1, Double_t width2, Double_t width3, Double_t width4, Double_t width5, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for 5-dimensional kernel PDF with binned interpolation from the sample of points in an NTuple with weight. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars1 1st variable name
        \param [in] vars2 2nd variable name
        \param [in] vars3 3rd variable name
        \param [in] vars4 4th variable name
        \param [in] vars5 5th variable name
        \param [in] weight weight
        \param [in] bins1 number of bins for the 1st variable for the binned interpolation
        \param [in] bins2 number of bins for the 2nd variable for the binned interpolation
        \param [in] bins3 number of bins for the 3rd variable for the binned interpolation
        \param [in] bins4 number of bins for the 4th variable for the binned interpolation
        \param [in] bins5 number of bins for the 5th variable for the binned interpolation
        \param [in] width1 kernel width for the 1st variable
        \param [in] width2 kernel width for the 2nd variable
        \param [in] width3 kernel width for the 3rd variable
        \param [in] width4 kernel width for the 4th variable
        \param [in] width5 kernel width for the 5th variable
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  const char *vars1, const char *vars2, const char *vars3, const char *vars4, const char *vars5, 
                  const char *weight, 
                  UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, UInt_t bins5, 
                  Double_t width1, Double_t width2, Double_t width3, Double_t width4, Double_t width5, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0,
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Constructor for kernel PDF with binned interpolation of arbitrary dimensionality from the sample of points in an NTuple. 
    /*! 
        \param [in] pdfName PDF name
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars vector of variable names. The size of vector should match the dimensionality of phase space or be larger by one.
        \param [in] binning vector of numbers of bins for the binned interpolation. The size of vector should match the dimensionality of phase space.
        \param [in] width vector of kernel widths. The size of vector should match the dimensionality of phase space.
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    BinnedKernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  std::vector<TString> &vars, 
                  std::vector<UInt_t> &binning, 
                  std::vector<Double_t> &width, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0, 
                  UInt_t maxEvents = 0,
                  UInt_t skipEvents = 0
                  );

    //! Destructor
    virtual ~BinnedKernelDensity();

    //! Write the binned PDF to a file. Can be used later by BinnedDensity. 
    //! Format of the file (text or ROOT) is decided based on the extension of the file name. 
    /*! 
        \param [in] fileName file name
    */ 
    void writeToFile(const char* fileName); 

    //! Write the binned PDF to a text file.
    /*! 
        \param [in] fileName file name
    */ 
    void writeToTextFile(const char* fileName); 

    //! Write the binned density into a ROOT file
    /*!
        \param [in] fileName output ROOT file name
    */ 
    void writeToRootFile(const char* fileName); 

    //! Calculate PDF density at the point. Multilinear binned interpolation is used which is calculated in the constructor, so this call is fast. 
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

    //! Set fractional mode of operation. 
    /*! 
      \param mode If mode==true, the density returned will be the ratio of the kernel and approximation density
    */ 
    void setFractionalMode(Bool_t mode = true) { m_fractionalMode = mode; }
    
    //! Normalise the PDF such that the average PDF value over the allowed phase space equals to 1
    void normalise(void); 

  private: 

    //! Common initialise method used by all constructors. 
    /*! 
        \param [in] thePhaseSpace phase space
        \param [in] tree ROOT NTuple 
        \param [in] vars vector of variable names. The size of vector should match the dimensionality of phase space.
        \param [in] binning vector of numbers of bins for the binned interpolation. The size of vector should match the dimensionality of phase space.
        \param [in] width vector of kernel widths. The size of vector should match the dimensionality of phase space.
        \param [in] approx approximation PDF. Uniform approximation is used for approx=0
        \param [in] toyEvents number of toy events for MC convolution of the approximation PDF. Use binned convolution if toyEvents=0
        \param [in] maxEvents maximum number of events to read from NTuple. Read all events if maxEvents=0
        \param [in] skipEvents number of NTuple events to skip from the beginning
    */ 
    void init(AbsPhaseSpace* thePhaseSpace, 
                  TTree* tree, 
                  std::vector<TString> &vars, 
                  std::vector<UInt_t> &binning, 
                  std::vector<Double_t> &width, 
                  AbsDensity* approx = 0, 
                  UInt_t toyEvents = 0, 
                  UInt_t maxEvents = 0, 
                  UInt_t skipEvents = 0 
                  );

    //! Convert an N-dimensional iterator vector into a linear bin index in the bin map
    /*! 
        \param [in] iter iterator vector
        \return bin index in the map
    */ 
    UInt_t iterToIndex( std::vector<UInt_t> &iter ); 

    //! Add a kernel density of a single data point to the binned map. 
    /*! 
        \param [in] map pointer to the binned map
        \param [in] point data point coordinates
        \param [in] weight point weight
    */ 
    void addToMap(std::vector<Double_t> &map, std::vector<Double_t> &point, Double_t weight = 1.); 

    void fillMapFromTree( TTree* tree, std::vector<TString> &vars, UInt_t maxEvents = 0, UInt_t skipEvents = 0);

    void fillMapFromDensity(AbsDensity* density, UInt_t toyEvents = 0); 

    //! Calculate the raw density using the binned map (estimated or approximation) at a given point
    /*! 
        \param [in] map reference to the bin map
        \param [in] x point
    */ 
    Double_t mapDensity(std::vector<Double_t> &map, std::vector<Double_t> &x);

    /// Bin map of estimated PDF
    std::vector<Double_t> m_map;

    /// Bin map of approximation PDF convolved with the kernel
    std::vector<Double_t> m_approxMap; 

    /// Vector of bin numbers in each variable
    std::vector<UInt_t> m_binning;

    /// Reference to phase space
    AbsPhaseSpace* m_phaseSpace; 

    /// Reference to approximation PDF
    AbsDensity* m_approxDensity; 

    /// Vector of kernel widths
    std::vector<Double_t> m_width;

    /// Cached dimensionality of the phase space
    UInt_t m_dim;

    /// Fractional mode flag
    Bool_t m_fractionalMode; 

};

#endif
