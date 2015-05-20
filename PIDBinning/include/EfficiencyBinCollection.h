#ifndef EfficiencyBinCollection_H
#define EfficiencyBinCollection_H

#include "EfficiencyBin.h"
#include <deque>
#include "TMath.h"
#include "functional"
#include <unordered_map>
#include <string>

typedef std::unordered_map<int, EfficiencyBin*> t_binmap;
typedef std::pair<int, EfficiencyBin*> t_ibin_pair;

class TH2Poly;

class EfficiencyBinCollection : public TNamed {
  public:
  EfficiencyBinCollection(std::string name, std::string title, Int_t nbins_x = 1, Double_t x1 = 0, Double_t x2 = 1, Int_t nbins_y = 1,
                          Double_t y1 = 0, Double_t y2 = 1)
      : TNamed(name.data(), title.data()), m_x1(x1), m_x2(x2), m_y1(y1), m_y2(y2), m_nbins_x(nbins_x), m_nbins_y(nbins_y) {
    EfficiencyBin* tmp = nullptr;
    for (int iX = 0; iX < m_nbins_x; iX++) {
      for (int iY = 0; iY < m_nbins_y; iY++) {
        Double_t eta_low = m_y1 + (m_y2 - m_y1) / m_nbins_y * iY;
        Double_t eta_high = m_y1 + (m_y2 - m_y1) / m_nbins_y * (iY + 1);
        Double_t p_low = m_x1 + (m_x2 - m_x1) / m_nbins_x * iX;
        Double_t p_high = m_x1 + (m_x2 - m_x1) / m_nbins_x * (iX + 1);
        AddBin(tmp = EfficiencyBin::RectangularBin(p_low, eta_low, p_high, eta_high));
        tmp->SetBinXID(iX);
        tmp->SetBinYID(iY);
        bin_map.insert(t_ibin_pair(iX + m_nbins_x* iY, tmp));
      }
    }
  };
  virtual ~EfficiencyBinCollection() {};

  // Prints all bins, contents and efficiencies to the stdout
  virtual void PrintBins();

  // Fills the collection, similar to normal histograms
  // Returns a pointer to the filled bin (no idea if that's needed at some point)
  virtual EfficiencyBin* Fill(Double_t x, Double_t y, bool _passed, Double_t _weight = 1.0);

  // Removes a bin from the collection
  // Returns bool if successful
  virtual bool RemoveBin(EfficiencyBin*);

  // Updates all efficiencies
  virtual void UpdateEfficiencies();

  // Makes something plotable
  virtual TH2Poly* MakePH2Poly(bool _binmap = false, bool funny_colors = false);

  // Function to fill the neighbour informations
  virtual void BuiltNeighbourhood();
  virtual TH1* MergeBins(Double_t _max_kappa);

  virtual Int_t GetNBins(){return m_bins.size();}

  virtual void CleanShapes();
  virtual void CleanNeighbours();

  private:
  /* data */

  // Adds bin to the list of bins
  // Class only supports rectangular bins so user is not allowed to
  // add bins manually.
  virtual void AddBin(EfficiencyBin*);
  virtual EfficiencyBin* FindBin(Double_t x, Double_t y);
  virtual EfficiencyBin* GetBin(Int_t iX, Int_t iY);

  std::deque<EfficiencyBin*> m_bins;
  //Bin to collect anything which might fall out of the actual range.
  EfficiencyBin *m_overflow = new EfficiencyBin();
  Bool_t IsIntersecting(EfficiencyBin* bin1, EfficiencyBin* bin2);
  Double_t m_x1 = 0.;
  Double_t m_x2 = 0.;
  Double_t m_y1 = 0.;
  Double_t m_y2 = 0.;
  Int_t m_nbins_x = 0.;
  Int_t m_nbins_y = 0.;
  t_binmap bin_map;
};

#endif
