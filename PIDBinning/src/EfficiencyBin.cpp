#include "EfficiencyBin.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TList.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/DistFunc.h"

/*virtual*/ void EfficiencyBin::AddBin(TGraph* _graph) {
  if (_graph == nullptr) {
    std::cout << "EfficiencyBin::EfficiencyBin(TGraph*) called with NULL pointer!" << std::endl;
    exit(-1);
  }
  Add(_graph);
}

/*virtual*/ void EfficiencyBin::AddBin(TMultiGraph* _graph) {
  if (_graph == nullptr) {
    std::cout << "EfficiencyBin::EfficiencyBin(TGraph*) called with NULL pointer!" << std::endl;
    exit(-1);
  }
  Add(_graph);
}

/*virtual*/ void EfficiencyBin::AddBin(EfficiencyBin* _effbin) {
  Add(_effbin);
  this->m_sum_passed += _effbin->m_sum_passed;
  this->m_sum_total += _effbin->m_sum_total;
  this->m_sum_Sq_passed += _effbin->m_sum_Sq_passed;
  this->m_sum_Sq_total += _effbin->m_sum_Sq_total;
  _effbin->mother_bin = this;
  // Flatten list of subbins
  for (auto& subbin : _effbin->m_sub_bins) {
    m_sub_bins.insert(subbin);
  }

  // Update the relations
  std::set<EfficiencyBin*> updated;
  for (auto& neigh : m_neighbours) {
    auto bin = neigh->GetTopLevel();
    if (bin == this) {
      continue;
    }
    updated.insert(bin);
  }
  for (auto& neigh : _effbin->m_neighbours) {
    auto bin = neigh->GetTopLevel();
    if (bin == this) {
      continue;
    }
    updated.insert(bin);
  }

  m_neighbours = updated;

  UpdateEfficiency();
}

//[>virtual<] Int_t EfficiencyBin::IsInside(Double_t x, Double_t y) const {
// for (auto& obj : m_bounds) {
// if (obj->IsInside(x, y)) {
// return true;
//};
//}
// return false;
//}

/*virtual*/ void EfficiencyBin::Fill(bool _passed, Double_t _weight) {
  m_sum_total += _weight;
  m_sum_Sq_total += _weight * _weight;
  if (_passed) {
    m_sum_passed += _weight;
    m_sum_Sq_passed += _weight * _weight;
  }
}

/*virtual*/ Double_t EfficiencyBin::Area() {
  Double_t area = 0.0;
  TList* gl = GetListOfGraphs();
  TGraph* g;
  TIter next(gl);
  while ((g = (TGraph*)next())) {
    area += g->Integral();
  }
  return area;
}

/*virtual*/ void EfficiencyBin::UpdateEfficiency() {
  if (m_sum_total <= 0) {
    m_efficiency = -0.000001;
    m_efficiency_error = -1.;
    m_efficiency_error_up = -1.;
    m_efficiency_error_down = -1.;
    return;
  }
  m_efficiency = m_sum_passed / m_sum_total;
  if (m_efficiency < 0. || m_efficiency > 1) {
    m_efficiency = -0.000001;
    m_efficiency_error = -1.;
    m_efficiency_error_up = -1.;
    m_efficiency_error_down= -1.;
    return;
  }

  m_efficiency_error = TMath::Sqrt(
      TMath::Abs(((1 - 2 * m_efficiency) * m_sum_Sq_passed + m_efficiency * m_efficiency * m_sum_Sq_total) / (m_sum_total * m_sum_total)));
  if (m_efficiency != 0 and m_efficiency != 1) {
    m_efficiency_error_down =
        m_efficiency - (m_sum_passed * ROOT::Math::fdistribution_quantile(0.16, 2 * m_sum_passed, 2 * (m_sum_total - m_sum_passed + 1))) /
                           ((m_sum_total - m_sum_passed) + 1 +
                            (m_sum_passed) * ROOT::Math::fdistribution_quantile(0.16, 2 * m_sum_passed, 2 * (m_sum_total - m_sum_passed + 1)));
    m_efficiency_error_up =
        ((m_sum_passed + 1) * ROOT::Math::fdistribution_quantile(1 - 0.16, 2 * (m_sum_passed + 1), 2 * (m_sum_total - m_sum_passed))) /
            ((m_sum_total - m_sum_passed) +
             (m_sum_passed + 1) * ROOT::Math::fdistribution_quantile(1 - 0.16, 2 * (m_sum_passed + 1), 2 * (m_sum_total - m_sum_passed))) -
        m_efficiency;
  } else if (m_efficiency == 0) {
    m_efficiency_error_down = 0;
    m_efficiency_error_up =
        ((m_sum_passed + 1) * ROOT::Math::fdistribution_quantile(1 - 0.16, 2 * (m_sum_passed + 1), 2 * (m_sum_total - m_sum_passed))) /
            ((m_sum_total - m_sum_passed) +
             (m_sum_passed + 1) * ROOT::Math::fdistribution_quantile(1 - 0.16, 2 * (m_sum_passed + 1), 2 * (m_sum_total - m_sum_passed))) -
        m_efficiency;

  } else {
    m_efficiency_error_down =
        m_efficiency - (m_sum_passed * ROOT::Math::fdistribution_quantile(0.16, 2 * m_sum_passed, 2 * (m_sum_total - m_sum_passed + 1))) /
                           ((m_sum_total - m_sum_passed) + 1 +
                            (m_sum_passed) * ROOT::Math::fdistribution_quantile(0.16, 2 * m_sum_passed, 2 * (m_sum_total - m_sum_passed + 1)));
    m_efficiency_error_up = 0;
  }
  //Just checking if the confidence intervals make sense
  if(m_efficiency_error_up < 0 ||  m_efficiency_error_down < 0){
    std::cout << "Uncertainty has the wrong sign!" << std::endl;
  }

}

/*static*/ EfficiencyBin* EfficiencyBin::RectangularBin(Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
  Double_t x[] = {x1, x1, x2, x2, x1};
  Double_t y[] = {y1, y2, y2, y1, y1};
  TGraph* g = new TGraph(5, x, y);
  return new EfficiencyBin(g);
}

/*static*/ EfficiencyBin* EfficiencyBin::SquareBin(Double_t x1, Double_t y1, Double_t edge) { return RectangularBin(x1, y1, x1 + edge, y1 + edge); }

Double_t AdaptiveBinning::Kappa(EfficiencyBin* bin1, EfficiencyBin* bin2) {
  Double_t eff1 = bin1->GetEfficiency();
  Double_t eff2 = bin2->GetEfficiency();
  Double_t err1 = 0.;
  Double_t err2 = 0.;
  if (eff1 < eff2) {
    err1 = bin1->GetEfficiencyErrorUp();
    err2 = bin2->GetEfficiencyErrorDown();
  } else {
    err1 = bin1->GetEfficiencyErrorDown();
    err2 = bin2->GetEfficiencyErrorUp();
  }

  return std::fabs(eff1 - eff2) / std::sqrt(err1 * err1 + err2 * err2);
}
/*virtual*/ void EfficiencyBin::UpdateNeighbours(){
    std::set<EfficiencyBin*> updated;
    for(auto & neigh : m_neighbours){
      updated.insert(neigh->GetTopLevel());
    }
    m_neighbours = updated;
}
