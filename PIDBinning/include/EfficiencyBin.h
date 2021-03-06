#ifndef EfficiencyBin_H
#define EfficiencyBin_H

#include "TMultiGraph.h"
#include <iostream>
#include <set>
#include "TF1.h"

class EfficiencyBinCollection;

class EfficiencyBin : public TMultiGraph {
  // class EfficiencyBin{
  public:
  virtual ~EfficiencyBin() {};

  // Function to add a bin to the list of bins in this object.
  virtual void AddBin(TGraph* _graph);
  virtual void AddBin(TMultiGraph* _graph);

  // Additional function in case of EfficiencyBins being added.
  // This flattens the list of graphs in the main object and
  // updates the number of events and Sumw2()
  virtual void AddBin(EfficiencyBin* _effbin);

  // Checks if the object is inside this bin.
  // virtual Int_t IsInside(Double_t x, Double_t y) const override;
  // virtual Int_t IsInside(Double_t x, Double_t y) const;

  // Some getter functions
  Double_t GetPassed() const { return m_sum_passed; }
  Double_t GetTotal() const { return m_sum_total; }
  Double_t GetPassedSumw2() const { return m_sum_Sq_total; }
  Double_t GetTotalSumw2() const { return m_sum_Sq_total; }
  Double_t GetEfficiency() const { return m_efficiency; }
  Double_t GetEfficiencyError() const { return m_efficiency_error; }
  Double_t GetEfficiencyErrorUp() const { return m_efficiency_error_up; }
  Double_t GetEfficiencyErrorDown() const { return m_efficiency_error_down; }

  Int_t GetBinXID(){return m_iX;}
  Int_t GetBinYID(){return m_iY;}
  void SetBinXID(Int_t iX){m_iX = iX;}
  void SetBinYID(Int_t iY){m_iY = iY;}

  // Fill this bin, flag indicates if the event passed the selection
  virtual void Fill(bool _passed = false, Double_t _weight = 1.0);
  // Updates the efficiencies. Uncertainty is -1000 if the bin is empty
  // and -1 if there is something weird due to weights. This ensures a
  // priority when merging bins.
  virtual void UpdateEfficiency();
  virtual void UpdateNeighbours();

  // Some static factory functions for typical shapes
  static EfficiencyBin* RectangularBin(Double_t x1, Double_t y1, Double_t x2, Double_t y2);
  static EfficiencyBin* SquareBin(Double_t x1, Double_t y1, Double_t edge);

  // All the functions we have to replace as our object is not actually a TGraph but a list of TGraphs...
  virtual Double_t Area();

  friend EfficiencyBinCollection;

  private:
  /* data */
  // Storage for all the boundary defining objects
  Double_t m_sum_passed = 0;
  Double_t m_sum_total = 0;
  Double_t m_sum_Sq_passed = 0;
  Double_t m_sum_Sq_total = 0;

  Double_t m_efficiency = 0;
  Double_t m_efficiency_error = 0;
  Double_t m_efficiency_error_up = 0;
  Double_t m_efficiency_error_down = 0;

  Int_t m_iX = 0;
  Int_t m_iY = 0;

  // For now, only some shapes are supported. So forbid manual construction.
  EfficiencyBin() {};
  EfficiencyBin(TGraph* _graph) : TMultiGraph() {
    if (_graph == nullptr) {
      std::cout << "EfficiencyBin::EfficiencyBin(TGraph*) called with NULL pointer!" << std::endl;
      exit(-1);
    }
    Add(_graph);
    //All bins are also a collection of subbins, always including themself
    m_sub_bins.insert(this);
  };

  // A reference to the mother bin if the bin got added to another somewhere.
  EfficiencyBin* mother_bin = nullptr;
  std::set<EfficiencyBin*> m_neighbours;
  std::set<EfficiencyBin*> m_sub_bins;
  EfficiencyBin* GetTopLevel() {
    if (mother_bin == nullptr) {
      return this;
    } else {
      return mother_bin->GetTopLevel();
    }
  }

  //Store the initial neighours for later shape merging
  EfficiencyBin *left = nullptr;
  EfficiencyBin *right = nullptr;
  EfficiencyBin *top = nullptr;
  EfficiencyBin *bottom = nullptr;
};

namespace AdaptiveBinning {
Double_t Kappa(EfficiencyBin* bin1, EfficiencyBin* bin2);
}

#endif
