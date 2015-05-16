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
  Double_t total = 0.;
  for (auto& bin : m_bins) {
    total += bin->GetTotal();
    std::cout << "Bin " << iBin << ": Total = " << bin->GetTotal() << ", Passed = " << bin->GetPassed() << ", Eff = " << bin->GetEfficiency()
              << " +- " << bin->GetEfficiencyError() << " neighbours = " << bin->m_neighbours.size() << std::endl;
    iBin++;
  }
  std::cout << "Total number of entries: " << total << std::endl;
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

/*virtual*/ void EfficiencyBinCollection::SmearEfficiencies() {
  for (auto& bin : m_bins) {
    bin->SmearEfficiency();
  }
}

// Makes something plotable
/*virtual*/ TH2Poly* EfficiencyBinCollection::MakePH2Poly(bool _binmap) {
  auto poly = new TH2Poly();
  // sort(begin(m_bins), end(m_bins), [](EfficiencyBin* left, EfficiencyBin* right){return left->Area() > right->Area();});
  Int_t iColor = 0;
  poly->SetFloat();
  for (auto& bin : m_bins) {
    auto iBin = poly->AddBin(bin);
    if (_binmap) {
      poly->SetBinContent(iBin, iColor);
      poly->SetBinError(iBin, 0.);
    } else {
      poly->SetBinContent(iBin, bin->GetEfficiency());
      poly->SetBinError(iBin, bin->GetEfficiencyError());
    }
    iColor++;
  }
  return poly;
}

// Bool_t IsIntersecting(Int_t bn, Double_t* x, Double_t* y, Double_t xclipl, Double_t xclipr, Double_t yclipb, Double_t yclipt) {
Bool_t EfficiencyBinCollection::IsIntersecting(EfficiencyBin* bin1, EfficiencyBin* bin2) {
  // auto x = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetX();
  // auto y = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetY();
  // auto bn = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetN();

  auto xrefl = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetX()[0];
  auto xrefr = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetX()[2];
  auto yrefb = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetY()[0];
  auto yreft = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetY()[2];

  auto xclipl = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetX()[0] - 0.00000001;
  auto yclipb = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetY()[0] - 0.00000001;
  auto xclipr = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetX()[2] + 0.00000001;
  auto yclipt = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetY()[2] + 0.00000001;

  // for (int counter = 0; counter < (bn - 1); counter++) {
  // if (x[counter] < xclipl || x[counter] > xclipr) {
  // continue;
  //}
  // if (y[counter] < yclipb || y[counter] > yclipt) {
  // continue;
  //}
  // return true;
  //}
  // return false;

  // We want to count the following setup as neighbouring:
  //  #
  // #C#
  //  #
  // Assume the clip bin is in the center
  Double_t xrefc = (xrefl + xrefr) / 2.;
  Double_t yrefc = (yrefb + yreft) / 2.;

  bool in_x = xrefc > xclipl && xrefc < xclipr;
  bool in_y = yrefc > yclipb && yrefc < yclipt;

  // Bin ontop of C
  if (in_x && yrefc > yclipt) {
    // Check if lower boundary is in C
    if (yrefb <= yclipt) {
      return true;
    }
  }

  // Bin below of C
  if (in_x && yrefc < yclipb) {
    // Check if upper boundary is in C
    if (yreft >= yclipb) {
      return true;
    }
  }

  // Bin right of C
  if (in_y && xrefc > xclipr) {
    // Check if left boundary is in C
    if (xrefl <= xclipr) {
      return true;
    }
  }

  // Bin left of C
  if (in_y && xrefc < xclipl) {
    // Check if right boundary is in C
    if (xrefr >= xclipl) {
      return true;
    }
  }

  return false;
}

/*virtual*/ void EfficiencyBinCollection::BuiltNeighbourhood() {
  for (size_t iBin = 0; iBin < m_bins.size(); iBin++) {
    for (size_t yBin = 0; yBin < m_bins.size(); yBin++) {
      if (iBin == yBin) {
        continue;
      }
      if (IsIntersecting(m_bins.at(iBin), m_bins.at(yBin))) {
        m_bins.at(iBin)->m_neighbours.insert(m_bins.at(yBin));
      }
    }
  }
}

/*virtual*/ void EfficiencyBinCollection::MergeBins(Double_t _max_kappa) {
  EfficiencyBin* top = nullptr;
  EfficiencyBin* bottom = nullptr;
  Double_t smallest_kappa = 99999999.;
  Double_t temp_kappa = 0.;
  for (int ismear = 0; ismear < 20; ismear++) {
    for (auto& bin : m_bins) {
      for (auto& neigh : bin->m_neighbours) {
        auto neigh_toplevel = neigh->GetTopLevel();
        temp_kappa = AdaptiveBinning::Kappa(bin, neigh_toplevel);
        if (temp_kappa < smallest_kappa) {
          top = bin;
          bottom = neigh_toplevel;
          smallest_kappa = temp_kappa;
        }
      }
    }
    SmearEfficiencies();
  }

  if (smallest_kappa < _max_kappa) {
    top->AddBin(bottom);
    RemoveBin(bottom);
    MergeBins(_max_kappa);
  }
}
