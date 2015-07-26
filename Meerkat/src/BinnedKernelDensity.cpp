#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include "AbsPhaseSpace.hh"
#include "AbsDensity.hh"
#include "BinnedKernelDensity.hh"

#include "Timer.hh"

BinnedKernelDensity::BinnedKernelDensity(const char* pdfName, 
                             AbsPhaseSpace* thePhaseSpace, 
                             TTree* tree, 
                             std::vector<TString> &vars, 
                             std::vector<UInt_t> &binning, 
                             std::vector<Double_t> &width, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfName) {
  init(thePhaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfName, 
                             AbsPhaseSpace* thePhaseSpace, 
                             TTree* tree, 
                             const char *vars1, 
                             UInt_t bins1, 
                             Double_t width1, 
                             AbsDensity* d, 
                             UInt_t toyEvents,
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfName) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  binning.push_back(bins1);
  width.push_back(width1);
                           
  init(thePhaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfName, 
                             AbsPhaseSpace* thePhaseSpace, 
                             TTree* tree, 
                             const char *vars1, 
                             const char *weight, 
                             UInt_t bins1, 
                             Double_t width1, 
                             AbsDensity* d, 
                             UInt_t toyEvents,
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfName) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(weight)); 
  binning.push_back(bins1);
  width.push_back(width1);
                           
  init(thePhaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, 
                             UInt_t bins1, UInt_t bins2, 
                             Double_t width1, Double_t width2, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  width.push_back(width1);
  width.push_back(width2);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, 
                             const char *weight, 
                             UInt_t bins1, UInt_t bins2, 
                             Double_t width1, Double_t width2, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(weight)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  width.push_back(width1);
  width.push_back(width2);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, const char *vars3, 
                             UInt_t bins1, UInt_t bins2, UInt_t bins3, 
                             Double_t width1, Double_t width2, Double_t width3, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(vars3)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  binning.push_back(bins3);
  width.push_back(width1);
  width.push_back(width2);
  width.push_back(width3);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, const char *vars3, 
                             const char *weight, 
                             UInt_t bins1, UInt_t bins2, UInt_t bins3, 
                             Double_t width1, Double_t width2, Double_t width3, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(vars3)); 
  vars.push_back(TString(weight)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  binning.push_back(bins3);
  width.push_back(width1);
  width.push_back(width2);
  width.push_back(width3);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, const char *vars3, const char *vars4, 
                             UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, 
                             Double_t width1, Double_t width2, Double_t width3, Double_t width4, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(vars3)); 
  vars.push_back(TString(vars4)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  binning.push_back(bins3);
  binning.push_back(bins4);
  width.push_back(width1);
  width.push_back(width2);
  width.push_back(width3);
  width.push_back(width4);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, const char *vars3, const char *vars4, 
                             const char *weight, 
                             UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, 
                             Double_t width1, Double_t width2, Double_t width3, Double_t width4, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(vars3)); 
  vars.push_back(TString(vars4)); 
  vars.push_back(TString(weight)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  binning.push_back(bins3);
  binning.push_back(bins4);
  width.push_back(width1);
  width.push_back(width2);
  width.push_back(width3);
  width.push_back(width4);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, const char *vars3, const char *vars4, const char *vars5, 
                             UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, UInt_t bins5, 
                             Double_t width1, Double_t width2, Double_t width3, Double_t width4, Double_t width5, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(vars3)); 
  vars.push_back(TString(vars4)); 
  vars.push_back(TString(vars5)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  binning.push_back(bins3);
  binning.push_back(bins4);
  binning.push_back(bins5);
  width.push_back(width1);
  width.push_back(width2);
  width.push_back(width3);
  width.push_back(width4);
  width.push_back(width5);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::BinnedKernelDensity(const char* pdfname, 
                             AbsPhaseSpace* thephaseSpace, 
                             TTree* tree, 
                             const char *vars1, const char *vars2, const char *vars3, const char *vars4, const char *vars5, 
                             const char *weight, 
                             UInt_t bins1, UInt_t bins2, UInt_t bins3, UInt_t bins4, UInt_t bins5, 
                             Double_t width1, Double_t width2, Double_t width3, Double_t width4, Double_t width5, 
                             AbsDensity* d, 
                             UInt_t toyEvents, 
                             UInt_t maxEvents, 
                             UInt_t skipEvents
                           ) : AbsDensity(pdfname) {
                           
  std::vector<TString> vars;
  std::vector<UInt_t> binning; 
  std::vector<Double_t> width; 
  
  vars.push_back(TString(vars1)); 
  vars.push_back(TString(vars2)); 
  vars.push_back(TString(vars3)); 
  vars.push_back(TString(vars4)); 
  vars.push_back(TString(vars5)); 
  vars.push_back(TString(weight)); 
  binning.push_back(bins1);
  binning.push_back(bins2);
  binning.push_back(bins3);
  binning.push_back(bins4);
  binning.push_back(bins5);
  width.push_back(width1);
  width.push_back(width2);
  width.push_back(width3);
  width.push_back(width4);
  width.push_back(width5);
                           
  init(thephaseSpace, tree, vars, binning, width, d, toyEvents, maxEvents, skipEvents); 
}

BinnedKernelDensity::~BinnedKernelDensity() {

}

/// Initialise method used by all constructors
void BinnedKernelDensity::init(AbsPhaseSpace* thephaseSpace, 
                    TTree* tree, 
                    std::vector<TString> &vars, 
                    std::vector<UInt_t> &binning, 
                    std::vector<Double_t> &width, 
                    AbsDensity* d, 
                    UInt_t toyEvents, 
                    UInt_t maxEvents, 
                    UInt_t skipEvents
                  ) {

  m_phaseSpace = thephaseSpace; 
  m_binning = binning; 
  m_width = width; 
  m_approxDensity = d; 
  m_dim = m_phaseSpace->dimensionality(); 
  
  m_fractionalMode = false; 

  printf("%20.20s INFO: Creating binned kernel density over %dD phase space\n", m_name, m_dim ); 
  
  if (m_binning.size() != m_dim) {
    printf("%20.20s ERROR: Dimensionality of phase space (%d) does not match binning vector size (%d)\n", 
           m_name, m_dim, (UInt_t)m_binning.size());
    abort(); 
  }
  
  if (m_width.size() != m_dim) {
    printf("%20.20s ERROR: Dimensionality of phase space (%d) does not match kernel width vector size (%d)\n", 
           m_name, m_dim, (UInt_t)m_width.size());
    abort(); 
  }
  
  UInt_t size = 1; 
  std::vector<UInt_t>::iterator i;
  for (i=m_binning.begin(); i!=m_binning.end(); i++) {
    size *= (*i);
  }

  printf("%20.20s INFO: Map size=%d\n", m_name, size);
  if (size > 20000000) {
    printf("%20.20s ERROR: Map size (%d) too large!\n", m_name, size);
    abort();
  }

  m_map.resize(size);
  m_approxMap.resize(size); 

  fillMapFromTree(tree, vars, maxEvents, skipEvents);
  fillMapFromDensity(m_approxDensity, toyEvents);

  normalise();

}

/// Calculate map index for a given iterator vector
UInt_t BinnedKernelDensity::iterToIndex( std::vector<UInt_t> &iter ) {
  UInt_t index = 0;
  Int_t n;
  for (n=m_dim-1; n>=0; n--) {
    if (n==(Int_t)m_dim-1) {
      index = iter[n]; 
    } else {
      index = index*m_binning[n] + iter[n];
    }
  }
  return index;
}

void BinnedKernelDensity::addToMap(std::vector<Double_t> &map, std::vector<Double_t> &point, Double_t weight) {

  // Fill the map
  std::vector<UInt_t> initBin(m_dim); 
  std::vector<UInt_t> finalBin(m_dim); 
  std::vector<Double_t> lowLimit(m_dim);
  std::vector<Double_t> coeff(m_dim); 

  UInt_t size = map.size(); 

  // Calculate the initial and final N-dim bins
  Int_t n;
  for (n=0; n<(Int_t)m_dim; n++) {
    Double_t x1 = point[n] - m_width[n]; 
    Double_t x2 = point[n] + m_width[n];

    Double_t low = m_phaseSpace->lowerLimit(n);
    Double_t up  = m_phaseSpace->upperLimit(n);
    Int_t i1 = (int)TMath::Ceil( (x1-low)/(up-low)*(m_binning[n]-1) );
    if (i1 < 0) i1 = 0; 
    if (i1 >= (Int_t)m_binning[n]) i1 = m_binning[n] - 1;
    Int_t i2 = (int)TMath::Floor( (x2-low)/(up-low)*(m_binning[n]-1) );
    if (i2 < 0) i2 = 0; 
    if (i2 >= (Int_t)m_binning[n]) i2 = m_binning[n] - 1;

    if (i1 > i2) {
      printf("%20.20s ERROR: i1(%d)>i2(%d) for n=%d, bin size is larger than kernel width!\n", m_name, i1, i2, n); 
      abort(); 
    }

    initBin[n] = i1;
    finalBin[n] = i2;
    
    lowLimit[n] = phaseSpace()->lowerLimit(n); 
    coeff[n] = ( phaseSpace()->upperLimit(n) - lowLimit[n] ) / ((Double_t)m_binning[n]-1); 
    lowLimit[n] -= point[n]; 
    
    coeff[n] /= m_width[n]; 
    lowLimit[n] /= m_width[n]; 
  }

  // Loop through the map and fill bins
  std::vector<UInt_t> iter(m_dim); 
  for (n=0; n<(Int_t)m_dim; n++) iter[n] = initBin[n];

  do {

    UInt_t index = iterToIndex( iter ); 

    if (index >= size) {
      printf("%20.20s ERROR: index (%d) is larger than array size (%d)\n", m_name, index, size); 
      abort(); 
    } else {
      Double_t sqsum = 0; 
      for (n=0; n<(Int_t)m_dim; n++) {
//        Double_t low = m_phaseSpace->lowerLimit(n);
//        Double_t up  = m_phaseSpace->upperLimit(n);
//        Double_t binx = low + (Double_t)iter[n]/((Double_t)m_binning[n]-1)*(up-low);
//        Double_t dx = binx/m_width[n]; 
        Double_t dx = lowLimit[n] + (Double_t)iter[n]*coeff[n];
        if (fabs(dx) < 1.) sqsum += dx*dx; 
      }
      if (sqsum < 1.) map[index] += weight*(1.-sqsum); 

    }

    // Increase iterator
    Bool_t run = 0; 
    for (n=0; n<(Int_t)m_dim; n++) {
      if (iter[n] < finalBin[n]) {
        iter[n]++; 
        run = 1; 
        break;
      } else {
        iter[n] = initBin[n];
      }
    }
    if (!run) break;

  } while(1); 

}

void BinnedKernelDensity::fillMapFromTree( TTree* tree, std::vector<TString> &vars, UInt_t maxEvents, UInt_t skipEvents) {

  if (vars.size() != m_dim && vars.size() != m_dim + 1) {
    printf("%20.20s ERROR: Number of TTree variables (%d) in tree \"%s\" does not correspond to phase space dimensionality (%d)\n", 
           m_name, (UInt_t)vars.size(), tree->GetName(), m_phaseSpace->dimensionality() ); 
    abort(); 
  }  

  if (vars.size() == m_dim + 1) {
    printf("%20.20s INFO: Using variable \"%s\" as weight\n", 
           m_name, vars[m_dim].Data() ); 
  }
  
  UInt_t nvars = vars.size(); 

  tree->ResetBranchAddresses();

  Long64_t nentries = tree->GetEntries();
  if (maxEvents > 0 && skipEvents + maxEvents < nentries) nentries = skipEvents + maxEvents; 

  printf("%20.20s INFO: Will read %lld events (skipping first %d)\n", 
           m_name, nentries-skipEvents, skipEvents ); 

  Long64_t i;

  std::vector<Float_t> varArray(nvars); 

  Int_t n; 
  for (n=0; n < (Int_t)nvars; n++) {
    printf("%20.20s INFO: Will read branch \"%s\" from tree \"%s\"\n", m_name, vars[n].Data(), tree->GetName()); 
    Int_t status = tree->SetBranchAddress(vars[n], &( varArray[n] ));
    if (status < 0) {
      printf("%20.20s WARNING: Error setting branch, status=%d\n", m_name, status); 
      abort(); 
    }
  }

  std::vector<Double_t> point(m_dim); 

  Long64_t nout = 0;
  
  set_timer(); 
  
  for(i=skipEvents; i<nentries; i++) {
    tree->GetEntry(i);
    for (n=0; n<(Int_t)m_dim; n++) {
      point[n] = varArray[n]; 
    }
    Double_t weight = 1.; 
    if (m_dim + 1 == nvars) 
      weight = varArray[m_dim]; 

    if (!m_phaseSpace->withinLimits( point )) {
      nout ++; 
//      printf("%20.20s WARNING: Ntuple point (", m_name); 
//      for (n=0; n<(Int_t)nvars; n++) {
//        printf("%f ", point[n]);
//      }
//      printf(") outside phase space\n");
    } else {
      addToMap(m_map, point, weight); 
    }
    
    if (i % 100 == 0 && timer(2)) {
      printf("%20.20s INFO: Read %lld/%lld events (%f%%), %lld out\n", m_name, i-skipEvents, nentries-skipEvents, 
             100.*float(i-skipEvents)/float(nentries-skipEvents), nout);
    }
  }

  printf("%20.20s INFO: %lld events read in from \"%s\", %lld out\n", m_name, nentries-nout, tree->GetName(), nout ); 

}

void BinnedKernelDensity::fillMapFromDensity(AbsDensity* theDensity, UInt_t toyEvents) {

  std::vector<Double_t> x(m_dim);
//  UInt_t size = m_approxMap.size(); 

  if (theDensity == 0) {
    printf("%20.20s INFO: Will use uniform density for approximation\n", m_name); 
  } else {
    printf("%20.20s INFO: Will use density \"%s\" for approximation\n", m_name, theDensity->name()); 
  }

  if (toyEvents == 0) {
    // Fill map in nodes of the binning
    printf("%20.20s INFO: Convolution of approx. density using rectangular grid\n", m_name); 

    // Create iterator vector
    std::vector<UInt_t> iter(m_dim); 
    Int_t j;
    for (j=0; j<(Int_t)m_dim; j++) {
      iter[j] = 0;
    }

    // Loop through the nodes
    Int_t index = 0; 

    set_timer(); 
    do {
      index++; 
      for (j=m_dim-1; j>=0; j--) {
        Double_t low = m_phaseSpace->lowerLimit(j);
        Double_t up  = m_phaseSpace->upperLimit(j);
        x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      }

      Double_t e = 1.; 
      if (theDensity) e = theDensity->density(x); 

      if ((index % 100) == 0 && timer(2))
        printf("%20.20s INFO: Index %d, density=%f\n", m_name, index, e); 

      if (m_phaseSpace->withinLimits(x)) addToMap(m_approxMap, x, e);

      // Increase iterator
      Bool_t run = 0; 
      for (j=0; j<(Int_t)m_dim; j++) {
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

  } else {
  
    // Fill map from random points
    
    printf("%20.20s INFO: Convolution of approx. density using MC with %d events\n", m_name, toyEvents); 
    Int_t i; 
    std::vector<Double_t> lower(m_dim); 
    std::vector<Double_t> coeff(m_dim); 
    for (i=0; i<(Int_t)m_dim; i++) {
      lower[i] = m_phaseSpace->lowerLimit(i);
      coeff[i] = m_phaseSpace->upperLimit(i) - lower[i];
    }
    
    set_timer(); 
    for (i=0; i<(Int_t)toyEvents; i++) {

      Bool_t success = 0;
      UInt_t t; 
      for (t = 0; t < m_maxTries; t++) {

        // Generate random point
        UInt_t var;
        for (var = 0; var < m_dim; var++) {
          x[var] = lower[var] + m_rnd.Rndm()*coeff[var];
        }

        Bool_t inPhsp = m_phaseSpace->withinLimits(x); 
        if (inPhsp) {
          success = 1;
          break;
        }
      }
      if (!success) {
        printf("%20.20s WARNING: failed to generate a point within phase space after %d tries\n", m_name, m_maxTries); 
        continue; 
      }

      Double_t e = 1; 
      if (theDensity) e = theDensity->density(x);

      if ((i % 100) == 0 && timer(1))
        printf("%20.20s INFO: toy event %d/%d (%f%%), density=%f\n", m_name, i, toyEvents, 100*float(i)/float(toyEvents), e); 

      addToMap(m_approxMap, x, e);

    }
  }
}

/// Write density map to a file, ROOT or text depending on extension
void BinnedKernelDensity::writeToFile(const char* filename) {
  if (TString(filename).EndsWith(".root") ) {
    writeToRootFile(filename); 
  } else {
    writeToTextFile(filename); 
  }
}

/// Write density map to a text file
void BinnedKernelDensity::writeToTextFile(const char* filename) {

  printf("%20.20s INFO: Writing binned density to text file \"%s\"\n", m_name, filename ); 

  FILE* file = fopen(filename, "w+"); 
  fprintf(file, "%d\n", m_dim);
  
  Int_t j; 
  for (j=0; j<(Int_t)m_dim; j++) {
    fprintf(file, "%d\n", m_binning[j]);
  }
  
  // Zero iterator vector
  std::vector<Double_t> x(m_dim);
  std::vector<UInt_t> iter(m_dim); 
  for (j=0; j<(Int_t)m_dim; j++) {
    iter[j] = 0;
  }

  UInt_t size = m_map.size(); 
  
  // Loop through the bins
  do {
    
    UInt_t index = 0;
    for (j=m_dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up  = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)m_dim-1) {
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
      fprintf(file, "%f %d\n", density(x), m_phaseSpace->withinLimits(x) );
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)m_dim; j++) {
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
void BinnedKernelDensity::writeToRootFile(const char* filename) {

  printf("%20.20s INFO: Writing binned density to ROOT file \"%s\"\n", m_name, filename ); 
  
  TDirectory* curr_dir = gDirectory;

  TFile file(filename, "RECREATE"); 
  TTree dimTree("DimTree", "DimTree"); 

  Int_t bins; 
  dimTree.Branch("bins",&bins,"bins/I"); 

  Int_t j; 
  for (j=0; j<(Int_t)m_dim; j++) {
    bins = m_binning[j]; 
    dimTree.Fill();
  }
  dimTree.Write(); 

  // Zero iterator vector
  std::vector<Double_t> x(m_dim);
  std::vector<UInt_t> iter(m_dim); 
  for (j=0; j<(Int_t)m_dim; j++) {
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
    for (j=m_dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up  = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)m_dim-1) {
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
      inphsp = m_phaseSpace->withinLimits(x);
      mapTree.Fill(); 
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)m_dim; j++) {
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

  gDirectory = curr_dir; 
}

/// Normalise density to have the average PDF value over the allowed phase space equal to 1
void BinnedKernelDensity::normalise(void) {

  printf("%20.20s INFO: Normalising density\n", m_name); 

  Double_t sum = 0; 
  UInt_t num = 0;

  // Zero iterator vector
  std::vector<Double_t> x(m_dim);
  std::vector<UInt_t> iter(m_dim); 
  UInt_t size = m_map.size(); 

  Int_t j;

  for (j=0; j<(Int_t)m_dim; j++) {
    iter[j] = 0;
  }

  // Loop through the bins
  do {
    
    UInt_t index = 0;
    for (j=m_dim-1; j>=0; j--) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up  = m_phaseSpace->upperLimit(j);
      x[j] = low + (Double_t)iter[j]/((Double_t)m_binning[j]-1)*(up-low);
      if (j==(Int_t)m_dim-1) {
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
      if (m_phaseSpace->withinLimits(x)) {
        sum += density(x); 
        num ++; 
      }
    }

    Bool_t run = 0; 
    for (j=0; j<(Int_t)m_dim; j++) {
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
  
  sum /= (Double_t)num; 
  
  printf("%20.20s INFO: Average PDF value before normalisation is %f\n", m_name, sum); 

  // Loop through the map and scale its entries
  for (j=0; j<(Int_t)size; j++) m_map[j] /= sum; 

}


Double_t BinnedKernelDensity::mapDensity(std::vector<Double_t> &map, std::vector<Double_t> &x) {

  Int_t j;
  std::vector<UInt_t> ivect(m_dim); 

  for (j=0; j<(Int_t)m_dim; j++) {
    Double_t low = m_phaseSpace->lowerLimit(j);
    Double_t up  = m_phaseSpace->upperLimit(j);
    Double_t xj = x[j]; 
    if (xj < low || xj > up) {
      return 0.;
    }
    Int_t ij = (Int_t)floor((xj-low)/(up-low)*(m_binning[j]-1));

    if (ij == (Int_t)m_binning[j]-1) ij--;

    ivect[j] = ij;
  }

  Double_t e = 0.;
  Double_t wsum = 0.;

  std::vector<UInt_t> iter(m_dim); 

  for (j=0; j<(Int_t)m_dim; j++) {
    iter[j] = 0;
  }

  // Loop through the vertices of the N-dim cube
  do {

    // Calculate weight
    Double_t weight = 1; 
    for (j=0; j<(Int_t)m_dim; j++) {
      Double_t low = m_phaseSpace->lowerLimit(j);
      Double_t up  = m_phaseSpace->upperLimit(j);

      Double_t xj1 = low + ((Double_t)ivect[j]   )/((Double_t)m_binning[j]-1.)*(up-low);
      Double_t xj2 = low + ((Double_t)ivect[j]+1.)/((Double_t)m_binning[j]-1.)*(up-low);

//      if (x[j] < xj1 || x[j] > xj2) {
//        printf("%20.20s WARNING: dim %d: x=%f, x1=%f, x2=%f\n", m_name, j, x[j], xj1, xj2); 
//      }

      Double_t fweight;
      if (iter[j] == 0) {
        fweight = 1. - (x[j]-xj1)/(xj2-xj1); 
      } else {
        fweight = (x[j]-xj1)/(xj2-xj1); 
      }
      if (fweight < 0.) {
        printf("%20.20s WARNING: dim %d: x=%f, weight=%f\n", m_name, j, x[j], fweight); 
        fweight = 0.;
      }
      if (fweight > 1.) {
        printf("%20.20s WARNING: dim %d: x=%f, weight=%f\n", m_name, j, x[j], fweight); 
        fweight = 1.;
      }
      
      weight *= fweight; 

//      printf("DEBUG:   Weight fraction: dim%d, i=%d, x=%f, x1=%f, x2=%f, fweight=%f\n", j, iter[j], x[j], xj1, xj2, fweight);

    }

    UInt_t index = 0;
    for (j=m_dim-1; j>=0; j--) {
      Int_t ij = ivect[j] + iter[j]; 
      if (j==(Int_t)m_dim-1) {
        index = ij; 
      } else {
        index = index*m_binning[j] + ij;
      }
    }

    e += weight*map[index]; 
    wsum += weight; 

//    printf("DEBUG: Weight=%f, index=%d, density=%f\n", weight, index, m_map[index]); 

    // Increment iterator
    Bool_t run = 0; 
    for (j=0; j<(Int_t)m_dim; j++) {
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

Double_t BinnedKernelDensity::density(std::vector<Double_t> &x) {
  Double_t a = mapDensity(m_approxMap, x); 
  if (a>0.) {
    if (m_approxDensity && !m_fractionalMode) { 
      return mapDensity(m_map, x)/a*m_approxDensity->density(x); 
    } else { 
      return mapDensity(m_map, x)/a; 
    }
  } else {
    return 0.; 
  }
}
