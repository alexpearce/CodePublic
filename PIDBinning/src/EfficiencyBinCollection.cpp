#include "EfficiencyBinCollection.h"
#include <iostream>
#include "TH2Poly.h"
#include "TGraph.h"
#include "TH1D.h"

#include <algorithm>
#include <iomanip>

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
  std::cout << "Bin" << std::setw(18) << "Total" << std::setw(18) << "Passed" << std::setw(18) << "Efficiency" << std::setw(18) << "Error"
            << std::setw(18) << "Up" << std::setw(18) << "Down" << std::setw(18) << "Neighbours"<< std::setw(18) << "Subbins" << std::endl;
  for (auto& bin : m_bins) {
    total += bin->GetTotal();
    std::cout << iBin << std::setw(18) << bin->GetTotal() << std::setw(18) << bin->GetPassed() << std::setw(18) << bin->GetEfficiency()
              << std::setw(18) << bin->GetEfficiencyError() << std::setw(18) << bin->GetEfficiencyErrorUp() << std::setw(18)
              << bin->GetEfficiencyErrorDown() << std::setw(18) << bin->m_neighbours.size()<< std::setw(18) << bin->m_sub_bins.size() << std::endl;
    iBin++;
  }
  std::cout << "Total number of entries: " << total << std::endl;
}

/*virtual*/ EfficiencyBin* EfficiencyBinCollection::Fill(Double_t x, Double_t y, bool _passed, Double_t _weight) {
  EfficiencyBin* selected_bin = FindBin(x, y);
  selected_bin->Fill(_passed, _weight);

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
/*virtual*/ TH2Poly* EfficiencyBinCollection::MakePH2Poly(bool _binmap, bool funny_colors) {
  auto poly = new TH2Poly(GetName(), GetTitle(), m_x1, m_x2, m_y1, m_y2);
  Int_t iColor = 0;
  for (auto& bin : m_bins) {
    auto iBin = poly->AddBin(bin);
    if (_binmap) {
      poly->SetBinContent(iBin, (funny_colors ? ((iColor % 2 == 0) ? iColor : GetNBins() - iColor) : iColor));
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

  auto xrefl = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetX()[0];
  auto xrefr = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetX()[2];
  auto yrefb = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetY()[0];
  auto yreft = ((TGraph*)bin1->GetListOfGraphs()->At(0))->GetY()[2];

  auto xclipl = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetX()[0] - 0.00000001;
  auto yclipb = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetY()[0] - 0.00000001;
  auto xclipr = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetX()[2] + 0.00000001;
  auto yclipt = ((TGraph*)bin2->GetListOfGraphs()->At(0))->GetY()[2] + 0.00000001;

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
  // for (size_t iBin = 0; iBin < m_bins.size(); iBin++) {
  // for (size_t yBin = 0; yBin < m_bins.size(); yBin++) {
  // if (iBin == yBin) {
  // continue;
  //}
  // if (IsIntersecting(m_bins.at(iBin), m_bins.at(yBin))) {
  // m_bins.at(iBin)->m_neighbours.insert(m_bins.at(yBin));
  //}
  //}
  //}

  for (auto& bin : m_bins) {
    Int_t iX = bin->GetBinXID();
    Int_t iY = bin->GetBinYID();
    if (iX > 0) {
      bin->left = GetBin(iX - 1, iY);
    }
    if (iX < m_nbins_x - 1) {
      bin->right = GetBin(iX + 1, iY);
    }
    if (iY > 0) {
      bin->bottom = GetBin(iX, iY - 1);
    }
    if (iY < m_nbins_y - 1) {
      bin->top = GetBin(iX, iY + 1);
    }
    for (auto& neigh : {bin->left, bin->right, bin->top, bin->bottom}) {
      if (neigh == nullptr) {
        continue;
      }
      bin->m_neighbours.insert(neigh);
    }
  }
}

/*virtual*/ TH1* EfficiencyBinCollection::MergeBins(Double_t _max_kappa) {
  EfficiencyBin* top = nullptr;
  EfficiencyBin* bottom = nullptr;
  Double_t smallest_kappa = 99999999.;
  Double_t temp_kappa = 0.;
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

  if (smallest_kappa < _max_kappa) {
    top->AddBin(bottom);
    RemoveBin(bottom);
    //Update the neighbouring relations of the merged bin
    for(auto & neigh : bottom->m_neighbours){
      neigh->UpdateNeighbours();
    }
    return MergeBins(_max_kappa);
  }
  TH1* hist = new TH1D("", "", 40, 0, 20);
  for (auto& bin : m_bins) {
    for (auto& neigh : bin->m_neighbours) {
      auto neigh_toplevel = neigh->GetTopLevel();
      temp_kappa = AdaptiveBinning::Kappa(bin, neigh_toplevel);
      hist->Fill(temp_kappa);
    }
  }

  return hist;
}

/*virtual*/ EfficiencyBin* EfficiencyBinCollection::GetBin(Int_t iX, Int_t iY) { return bin_map.find(iX + m_nbins_x * iY)->second; }

/*virtual*/ EfficiencyBin* EfficiencyBinCollection::FindBin(Double_t x, Double_t y) {
  // Do boundary checks:
  if (x < m_x1 || x >= m_x2 || y < m_y1 || y >= m_y2) {
    return m_overflow;
  }

  Int_t iX = (x - m_x1) / (m_x2 - m_x1) * m_nbins_x;
  Int_t iY = (y - m_y1) / (m_y2 - m_y1) * m_nbins_y;
  auto bin = GetBin(iX, iY);

  // Consistency check
  // This might fail due to floating point issues if the event is right on the corner.
  // However, it should still fill the correct bin.
  if (bin->IsInside(x, y) == false) {
    std::cout << "Warning! Due to floating point issues, x = " << x << ", y = " << y << " is reported to be outside the area." << std::endl;
    std::cout << "Please Check. Corners of the bin are:" << std::endl;
    auto graph = (TGraph*)bin->GetListOfGraphs()->At(0);
    for (int blub = 0; blub < graph->GetN(); blub++) {
      std::cout << "P" << blub << " x = " << graph->GetX()[blub] << ", y = " << graph->GetY()[blub] << std::endl;
    }
  };

  return bin->GetTopLevel();
}

/*virtual*/ void EfficiencyBinCollection::CleanShapes(){


}

/*virtual*/ void EfficiencyBinCollection::CleanNeighbours(){
  for(auto & bin : m_bins){
    bin->UpdateNeighbours();
  }
}
