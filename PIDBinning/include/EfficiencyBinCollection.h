#ifndef EfficiencyBinCollection_H
#define EfficiencyBinCollection_H

#include "EfficiencyBin.h"
#include <deque>
#include "TMath.h"

class TH2Poly;

class EfficiencyBinCollection {
  public:
  EfficiencyBinCollection() {};
  virtual ~EfficiencyBinCollection() {};

  // Adds bin to the list of bins
  virtual void AddBin(EfficiencyBin*);

  // Prints all bins, contents and efficiencies to the stdout
  virtual void PrintBins();

  // Fills the collection, similar to normal histograms
  // Returns a pointer to the filled bin (no idea if that's needed at some point)
  virtual EfficiencyBin* Fill(Double_t x, Double_t y, bool _passed, Double_t _weight = 1.0);

  // Removes a bin from the collection
  // Returns bool if successful
  virtual bool RemoveBin(EfficiencyBin*);

  // Updates all efficiencies
  virtual void UpdateEfficiencies(bool final = false);
  virtual void SmearEfficiencies();

  // Makes something plotable
  virtual TH2Poly* MakePH2Poly(bool _binmap = false);

  // Function to fill the neighbour informations
  virtual void BuiltNeighbourhood();
  virtual void MergeBins(Double_t _max_kappa);

  private:
  /* data */
  std::deque<EfficiencyBin*> m_bins;
  Bool_t IsIntersecting(EfficiencyBin* bin1, EfficiencyBin* bin2);
};


#endif
