#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include <iostream>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "RooRealVar.h"

#include "Roo1DBinnedDensity.hh"

#define MAX_VECTOR_SIZE 20000000

/// Constructor that reads from file
Roo1DBinnedDensity::Roo1DBinnedDensity(const char* pdfName, const char* title, RooRealVar& _x, Double_t down, Double_t up, const char* filename)
    : RooAbsPdf(pdfName, title), m_var("x", "x", this, _x) {
  strncpy(m_name, pdfName, 255);
  m_name[255] = 0;
  m_var_down = down;
  m_var_up = up;

  if (!TString(filename).EndsWith(".root")) {
    std::cout << "Invalid input filename " << filename << "! Terminating ..." << std::endl;
    std::exit(-1);
  }

  readFromRootFile(filename);
}

/// Read from ROOT file
void Roo1DBinnedDensity::readFromRootFile(const char* filename) {
  TFile file(filename);

  TTree* dimTree;
  TTree* mapTree;

  dimTree = (TTree*)gROOT->FindObject("DimTree");
  mapTree = (TTree*)gROOT->FindObject("MapTree");

  if (!dimTree) {
    printf("%20.20s ERROR: DimTree not found in file \"%s\"\n", m_name, filename);
    abort();
  }
  if (!mapTree) {
    printf("%20.20s ERROR: MapTree not found in file \"%s\"\n", m_name, filename);
    abort();
  }

  UInt_t dim = dimTree->GetEntries();

  printf("%20.20s INFO: Reading binned density over %dD phase space from ROOT file \"%s\"\n", m_name, dim, filename);

  if (dim != m_dim) {
    printf("%20.20s ERROR: Dimensionality of phase space (%d) does not match binning vector size (%d)\n", m_name, m_dim, dim);
    abort();
  }

  Int_t j;
  m_binning.resize(dim);
  Int_t nbins;
  dimTree->SetBranchAddress("bins", &nbins);
  for (j = 0; j < (Int_t)dim; j++) {
    dimTree->GetEvent(j);
    m_binning[j] = nbins;
  }

  UInt_t size = 1;
  std::vector<UInt_t>::iterator i;
  for (i = m_binning.begin(); i != m_binning.end(); i++) {
    size *= (*i);
  }

  printf("%20.20s INFO: Map size=%d\n", m_name, size);
  if (size > MAX_VECTOR_SIZE) {
    printf("%20.20s ERROR: Map size too large!\n", m_name);
    abort();
  }

  m_map.resize(size);

  // Zero iterator vector
  std::vector<Double_t> x(dim);
  std::vector<UInt_t> iter(dim);
  for (j = 0; j < (Int_t)dim; j++) {
    iter[j] = 0;
  }

  // Loop through the bins
  Float_t e;
  Bool_t inphsp;
  mapTree->SetBranchAddress("dens", &e);
  mapTree->SetBranchAddress("inphsp", &inphsp);
  if (mapTree->GetEntries() != size) {
    printf("%20.20s ERROR: Map size (%d) does not match number of entries in MapTree (%d)!\n", m_name, size, (int)mapTree->GetEntries());
    abort();
  }

  do {
    UInt_t index = 0;
    for (j = dim - 1; j >= 0; j--) {
      Double_t low = m_var_down;
      Double_t up = m_var_up;
      x[j] = low + (Double_t)iter[j] / ((Double_t)m_binning[j] - 1) * (up - low);
      if (j == (Int_t)dim - 1) {
        index = iter[j];
      } else {
        index = index * m_binning[j] + iter[j];
      }
    }

    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size);
      abort();
    } else {
      mapTree->GetEvent(index);
      m_map[index] = e;
    }

    Bool_t run = 0;
    for (j = 0; j < (Int_t)dim; j++) {
      if (iter[j] < m_binning[j] - 1) {
        iter[j]++;
        run = 1;
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;

  } while (1);

  file.Close();
}

Roo1DBinnedDensity::Roo1DBinnedDensity(const Roo1DBinnedDensity& pdf, const char* name)
    : RooAbsPdf(pdf, name),
      m_var("x", this, pdf.m_var),
      m_map(pdf.m_map),
      m_binning(pdf.m_binning),
      m_var_down(pdf.m_var_down),
      m_var_up(pdf.m_var_up) {}

Roo1DBinnedDensity::~Roo1DBinnedDensity() {}

Double_t Roo1DBinnedDensity::evaluate() const { return density(m_var); }

Double_t Roo1DBinnedDensity::density(Double_t _x) const {
  UInt_t dim = m_dim;

  Int_t j;
  std::vector<UInt_t> ivect(dim);
  std::vector<Double_t> x(dim);
  x.at(0) = _x;

  for (j = 0; j < (Int_t)dim; j++) {
    Double_t low = m_var_down;
    Double_t up = m_var_up;
    Double_t xj = x[j];
    if (xj < low || xj > up) {
      return 0.;
    }
    Int_t ij = (int)floor((xj - low) / (up - low) * (m_binning[j] - 1));

    if (ij == (Int_t)m_binning[j] - 1) ij--;

    ivect[j] = ij;
  }

  Double_t e = 0.;
  Double_t wsum = 0.;

  std::vector<UInt_t> iter(dim);

  for (j = 0; j < (Int_t)dim; j++) {
    iter[j] = 0;
  }

  // Loop through the vertices of the N-dim cube
  do {
    // Calculate weight
    Double_t weight = 1;
    for (j = 0; j < (Int_t)dim; j++) {
      Double_t low = m_var_down;
      Double_t up = m_var_up;

      Double_t xj1 = low + ((Double_t)ivect[j]) / ((Double_t)m_binning[j] - 1.) * (up - low);
      Double_t xj2 = low + ((Double_t)ivect[j] + 1.) / ((Double_t)m_binning[j] - 1.) * (up - low);

      if (x[j] < xj1 || x[j] > xj2) {
        printf("%20.20s WARNING: dim %d: x=%f, x1=%f, x2=%f\n", m_name, j, x[j], xj1, xj2);
      }

      Double_t fweight;
      if (iter[j] == 0) {
        fweight = 1. - (x[j] - xj1) / (xj2 - xj1);
      } else {
        fweight = (x[j] - xj1) / (xj2 - xj1);
      }
      weight *= fweight;

      //      printf("DEBUG:   Weight fraction: dim%d, i=%d, x=%f, x1=%f, x2=%f, fweight=%f\n", j, iter[j], x[j], xj1, xj2, fweight);
    }

    UInt_t index = 0;
    for (j = dim - 1; j >= 0; j--) {
      UInt_t ij = ivect[j] + iter[j];
      if (j == (Int_t)dim - 1) {
        index = ij;
      } else {
        index = index * m_binning[j] + ij;
      }
    }

    e += weight * m_map[index];
    wsum += weight;

    //    printf("DEBUG: Weight=%f, index=%d, density=%f\n", weight, index, m_map[index]);

    // Increment iterator
    Bool_t run = 0;
    for (j = 0; j < (Int_t)dim; j++) {
      if (iter[j] == 0) {
        iter[j]++;
        run = 1;
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;

  } while (1);

  //  printf("density=%f, wsum=%f\n", e, wsum);

  return e;
}

Int_t Roo1DBinnedDensity::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* ) const
{
  if(matchArgs(allVars, analVars, m_var)) return 1;
  return 0;
}

Double_t Roo1DBinnedDensity::analyticalIntegral(Int_t code, const char* ) const
{
  assert(code == 1);
  return (m_var_up - m_var_down);
}
