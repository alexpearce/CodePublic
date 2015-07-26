#ifndef Roo1DBinnedDensity_h
#define Roo1DBinnedDensity_h

#include "TMath.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"

#include <vector>

class RooRealVar;

/// A class that describes a previously created KDE (of any type), which
/// has been stored in the a root file. It implements a RooFit pdf
/// interface to access the created pdf in fits and plotting. Once created
/// the object owns all necessary information to be saved to and loaded from
/// root files and RooWorkspaces.
///
/// Only one-dimensional densities are supported!
///
/// author: Dominik Muller dominik.muller@cern.ch

// class Roo1DBinnedDensity : public AbsDensity {
class Roo1DBinnedDensity : public RooAbsPdf {
  public:
  Roo1DBinnedDensity(){}
  //! Constructor that reads the binned density from a file. The dimensionality of the density stored in the file should match the dimensionality of
  //the phase space.
  /*!
      \param [in] pdfName PDF name
      \param [in] thePhaseSpace phase space
      \param [in] fileName input file name
  */
  Roo1DBinnedDensity(const char* pdfName, const char* title, RooRealVar& _x, Double_t down, Double_t up, const char* fileName);

  //! Destructor
  virtual ~Roo1DBinnedDensity();
  Roo1DBinnedDensity(const Roo1DBinnedDensity& other, const char* name = 0);
  virtual Roo1DBinnedDensity* clone(const char* newname) const { return new Roo1DBinnedDensity(*this, newname); }

  //! Read the binned density from a ROOT file. The dimensionality of the density stored in the file should match the dimensionality of the phase
  //space.
  /*!
      \param [in] fileName input file name
  */
  void readFromRootFile(const char* fileName);

  Double_t density(Double_t x) const;

  protected:
  virtual Double_t evaluate() const;

  private:
  RooRealProxy m_var;

  //! Map of PDF values in bins
  std::vector<Double_t> m_map;

  //! Vector of bin numbers for each variable
  std::vector<UInt_t> m_binning;

  //! PDF name
  char m_name[256];

  //! Get rid of AbsPhaseSpace dependence in all
  // but the constructor
  UInt_t m_dim = 1;
  Double_t m_var_down = 0.;
  Double_t m_var_up = 1.;

  ClassDef(Roo1DBinnedDensity, 1);
};

#endif
