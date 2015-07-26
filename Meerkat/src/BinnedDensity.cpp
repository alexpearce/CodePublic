#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TRandom3.h"

#include "AbsPhaseSpace.hh"
#include "AbsDensity.hh"
#include "BinnedDensity.hh"

#include "Timer.hh"

#define MAX_VECTOR_SIZE 20000000

/// Constructor that fills bins from AbsDensity
BinnedDensity::BinnedDensity(const char* pdfName, 
                             AbsPhaseSpace* thePhaseSpace, 
                             AbsDensity* d, 
                             UInt_t bins1, 
                             UInt_t bins2, 
                             UInt_t bins3, 
                             UInt_t bins4, 
                             UInt_t bins5) : AbsDensity(pdfName) {
  std::vector<UInt_t> bins; 
  bins.push_back( bins1 ); 
  if (bins2 > 0) bins.push_back( bins2 ); 
  if (bins3 > 0) bins.push_back( bins3 ); 
  if (bins4 > 0) bins.push_back( bins4 ); 
  if (bins5 > 0) bins.push_back( bins5 ); 
  
  init(thePhaseSpace, bins, d); 
}

/// Constructor that fills bins from AbsDensity
BinnedDensity::BinnedDensity(const char* pdfName, 
                             AbsPhaseSpace* thePhaseSpace, 
                             std::vector<UInt_t> &binning, 
                             AbsDensity* d) : AbsDensity(pdfName) {
  init(thePhaseSpace, binning, d); 
}

/// Initialise method used by both constructors
void BinnedDensity::init(AbsPhaseSpace* thePhaseSpace, 
                    std::vector<UInt_t> &binning, 
                    AbsDensity* d) {

  m_phaseSpace = thePhaseSpace; 
  m_binning = binning; 
  m_density = d; 
  UInt_t dim = m_phaseSpace->dimensionality(); 

  printf("%20.20s INFO: Creating binned density over %dD phase space from density \"%s\"\n", m_name, dim, d->name() ); 
  
  if (m_binning.size() != m_phaseSpace->dimensionality()) {
    printf("%20.20s ERROR: Dimensionality of phase space (%d) does not match binning vector size (%d)\n", 
           m_name, m_phaseSpace->dimensionality(), (UInt_t)m_binning.size());
    abort(); 
  }
  
  UInt_t size = 1; 
  std::vector<UInt_t>::iterator i;
  for (i=m_binning.begin(); i!=m_binning.end(); i++) {
    size *= (*i);
  }
  
  printf("%20.20s INFO: Map size=%d\n", m_name, size); 
  if (size > MAX_VECTOR_SIZE) {
    printf("%20.20s ERROR: Map size (%d) too large!\n", m_name, size); 
    abort(); 
  }
  
  m_map.resize(size); 
  
  std::vector<Double_t> x(dim);
  std::vector<UInt_t> iter(dim); 
  Int_t j;
  for (j=0; j<(Int_t)dim; j++) {
    iter[j] = 0;
  }
  
  Double_t phspSum = 0.;
  UInt_t phspNum = 0;
  
  set_timer(); 
  
  do {
    UInt_t index = 0;
    for (j=dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)dim-1) {
        index = iter[j]; 
      } else {
        index = index*m_binning[j] + iter[j];
      }
//      if (iter[0] == 0)
//        printf("  Dim%d, bin%d, x=%f\n", j, iter[j], x[j]); 
    }

    Double_t e = m_density->density(x); 
    if (m_phaseSpace->withinLimits(x)) {
      phspSum += e;
      phspNum++;
    }
    
    if ((index % 100) == 0 && timer(2))
      printf("%20.20s INFO: Index %d, density=%f\n", m_name, index, e); 

    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size); 
      abort(); 
    } else {
      m_map[index] = e; 
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)dim; j++) {
      if (iter[j] < m_binning[j]-1) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;
    

  } while(1); 
  
  // Normalize map such that its average equals to 1.
  phspSum /= (Double_t)phspNum; 
  for (j = 0; j<(Int_t)size; j++) m_map[j] /= phspSum; 
}

/// Constructor that reads from file
BinnedDensity::BinnedDensity(const char* pdfName, AbsPhaseSpace* thePhaseSpace, const char* filename) : AbsDensity(pdfName) {
  m_phaseSpace = thePhaseSpace; 
  m_density = 0;

  readFromFile(filename); 
}

void BinnedDensity::readFromFile(const char* filename) {
  if (TString(filename).EndsWith(".root") ) {
    readFromRootFile(filename); 
  } else {
    readFromTextFile(filename); 
  }
}

/// Read from text file
void BinnedDensity::readFromTextFile(const char* filename) {

  FILE* file = fopen(filename, "r"); 
  if (!file) {
    printf("%20.20s ERROR: binned density file \"%s\" not found\n", m_name, filename ); 
    abort(); 
  }

  UInt_t dim; 
  if (fscanf(file, "%d", &dim) != 1) {
    printf("%20.20s ERROR: Phase space dimensionality not found in file \"%s\"\n", m_name, filename); 
    abort(); 
  } 

  printf("%20.20s INFO: Reading binned density over %dD phase space from text file \"%s\"\n", m_name, dim, filename ); 

  if (dim != m_phaseSpace->dimensionality()) {
    printf("%20.20s ERROR: Dimensionality of phase space (%d) does not match binning vector size (%d)\n", 
           m_name, m_phaseSpace->dimensionality(), dim);
    abort(); 
  }
  
  Int_t j; 
  m_binning.resize(dim); 
  for (j=0; j<(Int_t)dim; j++) {
    int nbins; 
    if (fscanf(file, "%d", &nbins) != 1) {
      printf("%20.20s ERROR: Number of bins not found in file \"%s\"\n", m_name, filename); 
      abort(); 
    }
    m_binning[j] = nbins; 
  }
  
  UInt_t size = 1; 
  std::vector<UInt_t>::iterator i;
  for (i=m_binning.begin(); i!=m_binning.end(); i++) {
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
  for (j=0; j<(Int_t)dim; j++) {
    iter[j] = 0;
  }
  
  // Loop through the bins
  
  set_timer(); 
  do {

    UInt_t index = 0;
    for (j=dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)dim-1) {
        index = iter[j]; 
      } else {
        index = index*m_binning[j] + iter[j];
      }
//      if (iter[0] == 0)
//        printf("  Dim%d, bin%d, x=%f\n", j, iter[j], x[j]); 
    }
    
    Double_t e;
    Int_t dummy; 
    if (fscanf(file, "%lf %d", &e, &dummy) != 2) {
      printf("%20.20s ERROR: Error reading map from file \"%s\", index %d\n", m_name, filename, index); 
      abort(); 
    }
    if ((index % 100) == 0 && timer(2)) printf("%20.20s INFO: Index %d, density=%f\n", m_name, index, e); 
    
    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size); 
      abort(); 
    } else {
      m_map[index] = e; 
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)dim; j++) {
      if (iter[j] < m_binning[j]-1) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;

  } while(1); 
  
  fclose(file); 
}

/// Read from ROOT file
void BinnedDensity::readFromRootFile(const char* filename) {

  TFile file(filename); 
//  if (!file) {
//    printf("%20.20s ERROR: binned density file \"%s\" not found\n", m_name, filename ); 
//    abort(); 
//  }
  
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

  printf("%20.20s INFO: Reading binned density over %dD phase space from ROOT file \"%s\"\n", m_name, dim, filename ); 

  if (dim != m_phaseSpace->dimensionality()) {
    printf("%20.20s ERROR: Dimensionality of phase space (%d) does not match binning vector size (%d)\n", 
           m_name, (int)m_phaseSpace->dimensionality(), dim);
    abort(); 
  }

  Int_t j; 
  m_binning.resize(dim); 
  Int_t nbins;
  dimTree->SetBranchAddress("bins", &nbins); 
  for (j=0; j<(Int_t)dim; j++) {
    dimTree->GetEvent(j); 
    m_binning[j] = nbins; 
  }

  UInt_t size = 1; 
  std::vector<UInt_t>::iterator i;
  for (i=m_binning.begin(); i!=m_binning.end(); i++) {
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
  for (j=0; j<(Int_t)dim; j++) {
    iter[j] = 0;
  }
  
  // Loop through the bins
  Float_t e; 
  Bool_t inphsp; 
  mapTree->SetBranchAddress("dens", &e); 
  mapTree->SetBranchAddress("inphsp", &inphsp); 
  if (mapTree->GetEntries() != size) {
    printf("%20.20s ERROR: Map size (%d) does not match number of entries in MapTree (%d)!\n", m_name, 
           size, (int)mapTree->GetEntries()); 
    abort(); 
  }
  
  set_timer(); 
  do {

    UInt_t index = 0;
    for (j=dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)dim-1) {
        index = iter[j]; 
      } else {
        index = index*m_binning[j] + iter[j];
      }
//      if (iter[0] == 0)
//        printf("  Dim%d, bin%d, x=%f\n", j, iter[j], x[j]); 
    }
    
//    if (fscanf(file, "%lf %d", &e, &dummy) != 2) {
//      printf("%20.20s ERROR: Error reading map from file \"%s\", index %d\n", m_name, filename, index); 
//      abort(); 
//    }
    if ((index % 100) == 0 && timer(2)) printf("%20.20s INFO: Index %d, density=%f\n", m_name, index, e); 
    
    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size); 
      abort(); 
    } else {
      mapTree->GetEvent(index);
      m_map[index] = e; 
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)dim; j++) {
      if (iter[j] < m_binning[j]-1) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;

  } while(1); 
  
  file.Close(); 
}


BinnedDensity::~BinnedDensity() {

}

void BinnedDensity::writeToFile(const char* filename) {
  if (TString(filename).EndsWith(".root") ) {
    writeToRootFile(filename); 
  } else {
    writeToTextFile(filename); 
  }
}

/// Write density map to file
void BinnedDensity::writeToTextFile(const char* filename) {

  printf("%20.20s INFO: Writing binned density to text file \"%s\"\n", m_name, filename ); 

  FILE* file = fopen(filename, "w+"); 
  
  UInt_t dim = m_phaseSpace->dimensionality(); 
  fprintf(file, "%d\n", dim);
  
  Int_t j; 
  for (j=0; j<(Int_t)dim; j++) {
    fprintf(file, "%d\n", m_binning[j]);
  }
  
  // Zero iterator vector
  std::vector<Double_t> x(dim);
  std::vector<UInt_t> iter(dim); 
  for (j=0; j<(Int_t)dim; j++) {
    iter[j] = 0;
  }
  
  UInt_t size = m_map.size(); 
  
  // Loop through the bins
  do {
    
    UInt_t index = 0;
    for (j=dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)dim-1) {
        index = iter[j]; 
      } else {
        index = index*m_binning[j] + iter[j];
      }
//      if (iter[0] == 0)
//        printf("  Dim%d, bin%d, x=%f\n", j, iter[j], x[j]); 
    }
    
    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size); 
      abort(); 
    } else {
      fprintf(file, "%f %d\n", m_map[index], m_phaseSpace->withinLimits(x) );
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)dim; j++) {
      if (iter[j] < m_binning[j]-1) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;

  } while(1); 
  fclose(file); 
}

/// Write density map to a ROOT file
void BinnedDensity::writeToRootFile(const char* filename) {

  printf("%20.20s INFO: Writing binned density to ROOT file \"%s\"\n", m_name, filename ); 

  TFile file(filename, "RECREATE"); 
  TTree dimTree("DimTree", "DimTree"); 

  Int_t bins; 
  dimTree.Branch("bins",&bins,"bins/I"); 
  UInt_t dim = m_phaseSpace->dimensionality(); 

  Int_t j; 
  for (j=0; j<(Int_t)dim; j++) {
    bins = m_binning[j]; 
    dimTree.Fill();
  }
  dimTree.Write(); 

  // Zero iterator vector
  std::vector<Double_t> x(dim);
  std::vector<UInt_t> iter(dim); 
  for (j=0; j<(Int_t)dim; j++) {
    iter[j] = 0;
  }

  UInt_t size = m_map.size(); 

  TTree mapTree("MapTree", "MapTree"); 
  
  Bool_t inphsp; 
  Float_t dens; 
  mapTree.Branch("dens",  &dens,"dens/F"); 
  mapTree.Branch("inphsp",&inphsp,"inphsp/B"); 
  
  // Loop through the bins
  do {
    
    UInt_t index = 0;
    for (j=dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up  = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)dim-1) {
        index = iter[j];
      } else {
        index = index*m_binning[j] + iter[j];
      }
//      if (iter[0] == 0)
//        printf("  Dim%d, bin%d, x=%f\n", j, iter[j], x[j]); 
    }

    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size); 
      abort(); 
    } else {
      dens = density(x); 
      inphsp = m_map[index];
      mapTree.Fill(); 
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)dim; j++) {
      if (iter[j] < m_binning[j]-1) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;

  } while(1); 

  mapTree.Write();
  file.Close(); 
}

Double_t BinnedDensity::density(std::vector<Double_t> &x) {
  
  UInt_t dim = m_phaseSpace->dimensionality(); 

  Int_t j;
  std::vector<UInt_t> ivect(dim); 
  
  for (j=0; j<(Int_t)dim; j++) {
    Double_t low = m_phaseSpace->lowerLimit(j);
    Double_t up = m_phaseSpace->upperLimit(j);
    Double_t xj = x[j]; 
    if (xj < low || xj > up) {
      return 0.;
    }
    Int_t ij = (int)floor((xj-low)/(up-low)*(m_binning[j]-1));

    if (ij == (Int_t)m_binning[j]-1) ij--;

    ivect[j] = ij;
  }
  
  Double_t e = 0.;
  Double_t wsum = 0.;
  
  std::vector<UInt_t> iter(dim); 
  
  for (j=0; j<(Int_t)dim; j++) {
    iter[j] = 0;
  }

  // Loop through the vertices of the N-dim cube
  do {

    // Calculate weight
    Double_t weight = 1; 
    for (j=0; j<(Int_t)dim; j++) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up = m_phaseSpace->upperLimit(j);
      
      Double_t xj1 = low + ((Double_t)ivect[j])/((Double_t)m_binning[j]-1.)*(up-low);
      Double_t xj2 = low + ((Double_t)ivect[j]+1.)/((Double_t)m_binning[j]-1.)*(up-low);

      if (x[j] < xj1 || x[j] > xj2) {
        printf("%20.20s WARNING: dim %d: x=%f, x1=%f, x2=%f\n", m_name, j, x[j], xj1, xj2); 
      }

      Double_t fweight;
      if (iter[j] == 0) {
        fweight = 1. - (x[j]-xj1)/(xj2-xj1); 
      } else {
        fweight = (x[j]-xj1)/(xj2-xj1); 
      }
      weight *= fweight; 
      
//      printf("DEBUG:   Weight fraction: dim%d, i=%d, x=%f, x1=%f, x2=%f, fweight=%f\n", j, iter[j], x[j], xj1, xj2, fweight);
      
    }

    UInt_t index = 0;
    for (j=dim-1; j>=0; j--) {
      UInt_t ij = ivect[j] + iter[j]; 
      if (j==(Int_t)dim-1) {
        index = ij; 
      } else {
        index = index*m_binning[j] + ij;
      }
    }

    e += weight*m_map[index]; 
    wsum += weight; 
    
//    printf("DEBUG: Weight=%f, index=%d, density=%f\n", weight, index, m_map[index]); 

    // Increment iterator
    Bool_t run = 0; 
    for (j=0; j<(Int_t)dim; j++) {
      if (iter[j] == 0) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = 0;
      }
    }
    if (!run) break;
    
  } while(1); 
  
//  printf("density=%f, wsum=%f\n", e, wsum); 
  
  return e; 
  
}
