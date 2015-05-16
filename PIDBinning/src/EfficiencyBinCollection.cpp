#include "EfficiencyBinCollection.h"
#include <iostream>
#include "TH2Poly.h"
#include "TGraph.h"

#include <algorithm>

/*virtual*/ void EfficiencyBinCollection::AddBin(EfficiencyBin* _effbin) {
  if (_effbin == nullptr) {
    std::cout << "EfficiencyBinCollection::AddBin(EfficiencyBin*) called with NULL pointer!" << std::endl;
    exit(-1);
  }
  m_bins.push_back(_effbin);
}

/*virtual*/ void EfficiencyBinCollection::PrintBins() {
  std::cout << "Printing bins in EfficiencyBinCollection:" << std::endl;
  int iBin = 0;
  for (auto& bin : m_bins) {
    std::cout << "Bin " << iBin << ": Total = " << bin->GetTotal() << ", Passed = " << bin->GetPassed() << ", Eff = " << bin->GetEfficiency()
              << " +- " << bin->GetEfficiencyError() << std::endl;
    iBin++;
  }
}

/*virtual*/ EfficiencyBin* EfficiencyBinCollection::Fill(Double_t x, Double_t y, bool _passed, Double_t _weight) {
  EfficiencyBin* selected_bin = nullptr;
  for (auto& bin : m_bins) {
    if (bin->IsInside(x, y)) {
      selected_bin = bin;
      selected_bin->Fill(_passed, _weight);
      break;
    }
  }

  return selected_bin;
}

/*virtual*/ bool EfficiencyBinCollection::RemoveBin(EfficiencyBin* _effbin) {
  auto frst = begin(m_bins);
  auto lst = end(m_bins);
  while (frst != lst) {
    if (*frst == _effbin) {
      m_bins.erase(frst);
    }
    frst++;
  }

  return frst == lst;
}

/*virtual*/ void EfficiencyBinCollection::UpdateEfficiencies() {
  for (auto& bin : m_bins) {
    bin->UpdateEfficiency();
  }
}

// Makes something plotable
/*virtual*/ TH2Poly* EfficiencyBinCollection::MakePH2Poly() {
  auto poly = new TH2Poly();
  // sort(begin(m_bins), end(m_bins), [](EfficiencyBin* left, EfficiencyBin* right){return left->Area() > right->Area();});
  for (auto& bin : m_bins) {
    auto iBin = poly->AddBin(bin);
    poly->SetBinContent(iBin, bin->GetEfficiency());
    poly->SetBinError(iBin, bin->GetEfficiencyError());
  }
  return poly;
}

//Bool_t IsIntersecting(Int_t bn, Double_t* x, Double_t* y, Double_t xclipl, Double_t xclipr, Double_t yclipb, Double_t yclipt) {
Bool_t EfficiencyBinCollection::IsIntersecting(EfficiencyBin* bin1, EfficiencyBin* bin2) {
  auto x = ((TGraph*) bin1->GetListOfGraphs()->At(0))->GetX();
  auto y = ((TGraph*) bin1->GetListOfGraphs()->At(0))->GetY();
  auto bn = ((TGraph*) bin1->GetListOfGraphs()->At(0))->GetN();

  auto xclipl = ((TGraph*) bin2->GetListOfGraphs()->At(0))->GetX()[0];
  auto yclipb = ((TGraph*) bin2->GetListOfGraphs()->At(0))->GetY()[0];
  auto xclipr = ((TGraph*) bin2->GetListOfGraphs()->At(0))->GetX()[0];
  auto yclipt = ((TGraph*) bin2->GetListOfGraphs()->At(0))->GetY()[0];

   for (int counter = 0; counter < (bn-1); counter++) {
     if(x[counter] < xclipl || x[counter]> xclipr){continue;}
     if(y[counter] < yclipb || y[counter]> yclipt){continue;}
     return true;
   }
   return false;
}
