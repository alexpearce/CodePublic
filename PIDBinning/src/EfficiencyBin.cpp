#include "EfficiencyBin.h"
#include "TEfficiency.h"
#include "TGraph.h"
#include "TList.h"

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
}

//[>virtual<] Int_t EfficiencyBin::IsInside(Double_t x, Double_t y) const {
  //for (auto& obj : m_bounds) {
    //if (obj->IsInside(x, y)) {
      //return true;
    //};
  //}
  //return false;
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
    TList *gl = GetListOfGraphs();
    TGraph *g;
    TIter next(gl);
    while ((g = (TGraph*) next())) {
      area += g->Integral();
    }
    return area;
}

/*virtual*/ void EfficiencyBin::UpdateEfficiency() {
  if (int(m_sum_total) == 0) {
    m_efficiency = 0.5;
    m_efficiency_error = -1000.;
    return;
  }
  m_efficiency = m_sum_passed / m_sum_total;
  if (m_efficiency < 0. || m_efficiency > 1) {
    m_efficiency = 0.5;
    m_efficiency_error = -1.;
    return;
  }

  Double_t upper = TEfficiency::AgrestiCoull(int(m_sum_total), int(m_sum_passed), 0.68, true);
  Double_t lower = TEfficiency::AgrestiCoull(int(m_sum_total), int(m_sum_passed), 0.68, false);
  m_efficiency_error = std::max(upper - m_efficiency, m_efficiency - lower);
}

/*static*/ EfficiencyBin* EfficiencyBin::RectangularBin(Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
  Double_t x[] = {x1, x1, x2, x2, x1};
  Double_t y[] = {y1, y2, y2, y1, y1};
  TGraph* g = new TGraph(5, x, y);
  return new EfficiencyBin(g);
}

/*static*/ EfficiencyBin* EfficiencyBin::SquareBin(Double_t x1, Double_t y1, Double_t edge) { return RectangularBin(x1, y1, x1 + edge, y1 + edge); }
