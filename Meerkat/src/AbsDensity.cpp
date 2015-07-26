#ifndef Roo1DBinnedSensity_h
#define Roo1DBinnedSensity_h

#include <stdio.h>
#include <vector>

#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

#include "AbsPhaseSpace.hh"
#include "AbsDensity.hh"

#include "Timer.hh"

AbsDensity::AbsDensity(const char* pdfName) {
  m_maxTries = 100000;
  m_majorant = 0. ;
  strncpy(m_name, pdfName, 255);
  m_name[255] = 0;  
}

AbsDensity::~AbsDensity() {

}

void AbsDensity::slice(std::vector<Double_t> &x, UInt_t num, TH1F* hist, Bool_t printout) {

  std::vector<Double_t> point = x; 

  UInt_t bins = hist->GetNbinsX(); 

  if ( printout )
    printf("%20.20s INFO: filling 1D slice in variable %d, hist \"%s\" (%d bins)\n", 
         m_name, num, hist->GetName(), bins ); 
  
  UInt_t i;
  set_timer();
  for (i=1; i<=bins; i++) {
    if (printout && timer(1)) printf("%20.20s INFO: filling bin %d\n", m_name, i);
    Double_t xi = hist->GetBinCenter(i); 
    point[num] = xi; 
    Double_t d = 0;
    if (phaseSpace()->withinLimits(point)) 
      d = density(point); 
    if (d<0) d = 0;
    hist->SetBinContent(i, d); 
  }
}

void AbsDensity::slice(std::vector<Double_t> &x, UInt_t numx, UInt_t numy, TH2F* hist, Bool_t printout, Bool_t inPhaseSpace) {

  std::vector<Double_t> point = x; 

  UInt_t xbins = hist->GetNbinsX(); 
  UInt_t ybins = hist->GetNbinsY(); 

  if (printout) 
    printf("%20.20s INFO: filling 2D slice in variables %d and %d, hist \"%s\" (%dx%d bins)\n", 
         m_name, numx, numy, hist->GetName(), xbins, ybins ); 
  
  UInt_t ix, iy;
  set_timer(); 
  for (ix=1; ix<=xbins; ix++) {
    if (printout && timer(1)) printf("%20.20s INFO: filling row %d\n", m_name, ix);
    for (iy=1; iy<=ybins; iy++) {
      Double_t xi = hist->GetXaxis()->GetBinCenter(ix); 
      Double_t yi = hist->GetYaxis()->GetBinCenter(iy); 
      point[numx] = xi; 
      point[numy] = yi; 
      Double_t d = 0;
      if (!inPhaseSpace || phaseSpace()->withinLimits(point)) 
        d = density(point); 
      if (d<0) d = 0;
      hist->SetBinContent(ix, iy, d);
    }
  }
}

void AbsDensity::project(TH1F* hist, Bool_t printout) {
  if (phaseSpace()->dimensionality() != 1) {
    printf("%20.20s ERROR: Projection of %dD density to 1D histogram is not supported\n", m_name, 
           phaseSpace()->dimensionality());
    abort();
  }
  std::vector<Double_t> x(1); 
  x[0] = 0;
  slice(x, 0, hist, printout); 
}

void AbsDensity::project(TH2F* hist, Bool_t printout, Bool_t inPhaseSpace) {
  if (phaseSpace()->dimensionality() != 2) {
    printf("%20.20s ERROR: Projection of %dD density to 2D histogram is not supported\n", m_name, 
           phaseSpace()->dimensionality());
    abort();
  }
  std::vector<Double_t> x(2); 
  x[0] = 0;
  x[1] = 0;
  slice(x, 0, 1, hist, printout, inPhaseSpace); 
}

Double_t AbsDensity::generate(std::vector<Double_t> &x) {

  Bool_t success = 0;
  Double_t d = 0;
  
  UInt_t dimensionality = phaseSpace()->dimensionality(); 
  UInt_t t; 

  Bool_t estimateMajorant = (m_majorant <= 0); 
  if (estimateMajorant) {
    printf("%20.20s INFO: Generating toy MC distribution\n", m_name);
    printf("%20.20s INFO: Majorant will be estimated with %d tries\n", m_name, m_maxTries);
  }

  for (t = 0; t < m_maxTries; t++) {

    // Generate random point
    UInt_t var;
    for (var = 0; var < dimensionality; var++) {
      Double_t lowerLimit = phaseSpace()->lowerLimit(var);
      Double_t upperLimit = phaseSpace()->upperLimit(var);
      x[var] = lowerLimit + m_rnd.Rndm()*(upperLimit-lowerLimit);
    }

    Bool_t inPhsp = phaseSpace()->withinLimits(x); 
    if (inPhsp) {
      Double_t y = m_majorant*m_rnd.Rndm();
      d = density(x); 
      if (d > m_majorant) {
        if (!estimateMajorant)
          printf("%20.20s WARNING: Updating majorant: %f -> %f\n", m_name, m_majorant, 1.1*d);
        m_majorant = 1.1*d; 
      } 
      if (d > y) {
        success = 1; 
        if (!estimateMajorant) break;
      }
    }
  }
  if (!success) {
    printf("%20.20s WARNING: failed to generate a point within phase space after %d tries\n", m_name, m_maxTries); 
    return 0;
  }
  
  if (estimateMajorant) {
    printf("%20.20s INFO: Estimated majorant = %f\n", m_name, m_majorant);
  }
  
  return d; 
}

void AbsDensity::generate(TNtuple* tree, UInt_t numEvents) {

  Float_t array[11]; 
  UInt_t dimensionality = phaseSpace()->dimensionality(); 
  
  if (dimensionality > 10) {
    printf("%20.20s ERROR: Generation is not supported for more than 10 dimensions\n", m_name);
    abort();
  }

  std::vector<Double_t> x(dimensionality); 
  UInt_t i; 
  set_timer(); 
  for (i=0; i<numEvents; i++) { 
    Double_t eff = generate(x); 

    UInt_t n; 
    for (n=0; n<dimensionality; n++) array[n] = x[n]; 
    array[dimensionality] = eff; 
    
    tree->Fill(array); 
    
    if (i % 100 == 0 && timer(2)) printf("%20.20s INFO: Ntuple event %d/%d (%f%%)\n", 
                                         m_name, i, numEvents, 100.*(Double_t)i/(Double_t)numEvents);
  } 
  tree->Write(); 
}

#endif
