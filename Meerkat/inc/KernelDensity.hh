#ifndef KERNEL_DENSITY
#define KERNEL_DENSITY

#include "AbsDensity.hh"
#include "AbsPhaseSpace.hh"

#include "TMath.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"

#include <vector>

typedef std::vector<Double_t> TPhspVector; 
typedef std::vector<TPhspVector> TCell; 

/// Class that describes the unbinned kernel density. 

class KernelDensity : public AbsDensity {

  public: 

    KernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace,  
                  std::vector<Double_t> &width, 
                  UInt_t approxSize, 
                  AbsDensity* approxDensity = 0);

    KernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace,  
                  UInt_t approxSize, 
                  Double_t width1, 
                  Double_t width2 = 0, 
                  Double_t width3 = 0, 
                  Double_t width4 = 0, 
                  Double_t width5 = 0);

    KernelDensity(const char* pdfName, 
                  AbsPhaseSpace* thePhaseSpace,  
                  UInt_t approxSize, 
                  AbsDensity* approxDensity, 
                  Double_t width1, 
                  Double_t width2 = 0, 
                  Double_t width3 = 0, 
                  Double_t width4 = 0, 
                  Double_t width5 = 0);

    virtual ~KernelDensity(); 

    void setWidth(std::vector<Double_t> &width);

    Bool_t generateApproximation(UInt_t approxSize);

    Bool_t readTuple(TTree* tree, std::vector<TString> &vars, UInt_t maxEvents = 0);

    Bool_t readTuple(TTree* tree, const char* var1, UInt_t maxEvents = 0);

    Bool_t readTuple(TTree* tree, const char* var1, const char* var2, UInt_t maxEvents = 0);

    Bool_t readTuple(TTree* tree, const char* var1, const char* var2, 
                                  const char* var3, UInt_t maxEvents = 0);

    Bool_t readTuple(TTree* tree, const char* var1, const char* var2, 
                                  const char* var3, const char* var4, UInt_t maxEvents = 0);

    Bool_t readTuple(TTree* tree, const char* var1, const char* var2, 
                                  const char* var3, const char* var4, 
                                  const char* var5, UInt_t maxEvents = 0);

    Bool_t readTuple(TTree* tree, const char* var1, const char* var2, 
                                  const char* var3, const char* var4, 
                                  const char* var5, const char* var6, UInt_t maxEvents = 0);

    Double_t density(std::vector<Double_t> &x);

    AbsPhaseSpace* phaseSpace() { return m_phaseSpace; }

  private: 

    UInt_t numCells(void); 

    Int_t cellIndex(std::vector<Double_t> &x); 

    Double_t rawDensity(std::vector<Double_t> &x, std::vector<TCell> &vector); 

    AbsPhaseSpace* m_phaseSpace; 

    AbsDensity* m_approxDensity; 

    TPhspVector m_width;

    std::vector<TCell> m_apprVector; 

    std::vector<TCell> m_dataVector; 

}; 

#endif
