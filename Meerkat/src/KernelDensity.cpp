#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TMath.h"
#include "TRandom3.h"

#include "AbsPhaseSpace.hh"
#include "AbsDensity.hh"
#include "KernelDensity.hh"

KernelDensity::KernelDensity(const char* pdfname, 
                  AbsPhaseSpace* thephaseSpace, 
                  UInt_t apprSize, 
                  Double_t width1, 
                  Double_t width2, 
                  Double_t width3, 
                  Double_t width4, 
                  Double_t width5) : AbsDensity(pdfname) {
  std::vector<Double_t> width_list; 

  width_list.push_back( width1 ); 
  if (width2 > 0) width_list.push_back( width2 ); 
  if (width3 > 0) width_list.push_back( width3 ); 
  if (width4 > 0) width_list.push_back( width4 ); 
  if (width5 > 0) width_list.push_back( width5 ); 

  m_phaseSpace = thephaseSpace;
  m_approxDensity = 0;
  m_maxTries = 10000;

  if (width_list.size() != m_phaseSpace->dimensionality()) {
    printf("%20.20s WARNING: Number of non-zero elements in width list (%d) does not match phase space dimensionality (%d)\n", 
           m_name, (UInt_t)width_list.size(), m_phaseSpace->dimensionality());
    abort(); 
  }

  setWidth(width_list);
  generateApproximation(apprSize);

}

KernelDensity::KernelDensity(const char* pdfname, 
                  AbsPhaseSpace* thephaseSpace, 
                  UInt_t apprSize, 
                  AbsDensity* approxDensity, 
                  Double_t width1, 
                  Double_t width2, 
                  Double_t width3, 
                  Double_t width4, 
                  Double_t width5) : AbsDensity(pdfname)  {
  std::vector<Double_t> width_list; 

  width_list.push_back( width1 ); 
  if (width2 > 0) width_list.push_back( width2 ); 
  if (width3 > 0) width_list.push_back( width3 ); 
  if (width4 > 0) width_list.push_back( width4 ); 
  if (width5 > 0) width_list.push_back( width5 ); 

  m_phaseSpace = thephaseSpace;
  m_approxDensity = approxDensity;
  m_maxTries = 10000;

  if (width_list.size() != m_phaseSpace->dimensionality()) {
    printf("%20.20s ERROR: Number of non-zero elements in width list (%d) does not match phase space dimensionality (%d)\n", 
           m_name, (UInt_t)width_list.size(), m_phaseSpace->dimensionality());
    abort(); 
  }
  
  if (approxDensity && approxDensity->phaseSpace() != m_phaseSpace) {
    printf("%20.20s ERROR: Phase space definitions for Kernel density and its approximation differ\n", m_name); 
    abort(); 
  }

  setWidth(width_list);
  generateApproximation(apprSize);

}

KernelDensity::KernelDensity(const char* pdfName, 
                             AbsPhaseSpace* thePhaseSpace, 
                             std::vector<Double_t> &width, 
                             UInt_t apprSize, 
                             AbsDensity* approxDensity) : AbsDensity(pdfName)  {

  m_phaseSpace = thePhaseSpace;
  m_approxDensity = approxDensity;
  m_maxTries = 10000;

  if (width.size() != m_phaseSpace->dimensionality()) {
    printf("%20.20s ERROR: Number of elements in width list (%d) does not match phase space dimensionality (%d)\n", 
           m_name, (UInt_t)width.size(), m_phaseSpace->dimensionality());
    abort(); 
  }
  
  if (approxDensity && approxDensity->phaseSpace() != m_phaseSpace) {
    printf("%20.20s ERROR: Phase space definitions for Kernel density and its approximation differ\n", m_name); 
    abort(); 
  }

  setWidth(width);
  generateApproximation(apprSize);
}

KernelDensity::~KernelDensity() {

}

/// Set kernel width
void KernelDensity::setWidth(std::vector<Double_t> &width) {

  UInt_t n = 0; 
  std::vector<Double_t>::const_iterator i; 
  for (i=width.begin(); i != width.end(); i++) {
    n++; 
    if (*i <= 0.) {
      printf("%20.20s ERROR: Kernel width for variable %d less or equal to zero\n", m_name, n);
      abort(); 
    }
  }

  m_width = width;

}

/// Create normalisation vector
Bool_t KernelDensity::generateApproximation(UInt_t apprSize) {
  
  UInt_t dimensionality = m_phaseSpace->dimensionality(); 
  
  UInt_t cells = numCells(); 
  
  m_apprVector.resize(cells); 
  
  Int_t cell;
  for (cell = 0; cell < (Int_t)cells; cell++)
    m_apprVector[cell].reserve(apprSize/cells); 

  std::vector<Double_t> point(dimensionality); 

  printf("%20.20s INFO: Generating approximation sample for %d-dim distribution, %d points, %d cells\n", 
         m_name, dimensionality, apprSize, cells);

  UInt_t i; 
  for (i=0; i<apprSize; i++) {

    if (!m_approxDensity) {
      UInt_t t; 
      Bool_t success = 0;
      for (t = 0; t < m_maxTries; t++) {

        // Generate random point
        UInt_t var;
        for (var = 0; var < dimensionality; var++) {
          Double_t lowerLimit = m_phaseSpace->lowerLimit(var);
          Double_t upperLimit = m_phaseSpace->upperLimit(var);
          point[var] = lowerLimit + m_rnd.Rndm()*(upperLimit-lowerLimit);
        }

        Bool_t inPhsp = m_phaseSpace->withinLimits(point); 
        if (inPhsp) {
          success = 1;
          break;
        }
      }
      if (!success) {
        printf("%20.20s WARNING: failed to generate a point within phase space after %d tries\n", m_name, m_maxTries); 
        return 0;
      }

    } else {
      m_approxDensity->generate(point); 
    }

    cell = cellIndex(point); 
    if (cell < (Int_t)cells && cell >= 0)
      m_apprVector[cell].push_back(point); 
    else if (cell > 0) {
      printf("%20.20s WARNING: cell number %d exceeds vector size %d\n", m_name, cell, cells); 
      abort(); 
    }

    if (i % 1000 == 0) printf("%20.20s INFO: %d%% done (%d/%d)\n", m_name, (100*i/apprSize), i, apprSize); 

  }

  return 1;
}

Bool_t KernelDensity::readTuple(TTree* tree, const char* var1, const char* var2, 
                                const char* var3, const char* var4, 
                                const char* var5, const char* var6, UInt_t maxEvents) {

  UInt_t dim = phaseSpace()->dimensionality(); 
  std::vector<TString> varList(dim); 
  
  if (dim>=1 && var1!=0) varList[0] = TString(var1); 
  if (dim>=2 && var2!=0) varList[1] = TString(var2); 
  if (dim>=3 && var3!=0) varList[2] = TString(var3); 
  if (dim>=4 && var4!=0) varList[3] = TString(var4); 
  if (dim>=5 && var5!=0) varList[4] = TString(var5); 
  if (dim>=6 && var6!=0) varList[5] = TString(var6);
  
  return readTuple(tree, varList, maxEvents); 

}

Bool_t KernelDensity::readTuple(TTree* tree, const char* var1, const char* var2, 
                                const char* var3, const char* var4, 
                                const char* var5, UInt_t maxEvents) {

  UInt_t dim = phaseSpace()->dimensionality(); 
  std::vector<TString> varList(dim); 
  
  if (dim>=1 && var1!=0) varList[0] = TString(var1); 
  if (dim>=2 && var2!=0) varList[1] = TString(var2); 
  if (dim>=3 && var3!=0) varList[2] = TString(var3); 
  if (dim>=4 && var4!=0) varList[3] = TString(var4); 
  if (dim>=5 && var5!=0) varList[4] = TString(var5); 
  
  return readTuple(tree, varList, maxEvents); 

}

Bool_t KernelDensity::readTuple(TTree* tree, const char* var1, const char* var2, 
                                const char* var3, const char* var4, 
                                UInt_t maxEvents) {

  UInt_t dim = phaseSpace()->dimensionality(); 
  std::vector<TString> varList(dim); 
  
  if (dim>=1 && var1!=0) varList[0] = TString(var1); 
  if (dim>=2 && var2!=0) varList[1] = TString(var2); 
  if (dim>=3 && var3!=0) varList[2] = TString(var3); 
  if (dim>=4 && var4!=0) varList[3] = TString(var4); 
  
  return readTuple(tree, varList, maxEvents); 

}

Bool_t KernelDensity::readTuple(TTree* tree, const char* var1, const char* var2, 
                                const char* var3, UInt_t maxEvents) {

  UInt_t dim = phaseSpace()->dimensionality(); 
  std::vector<TString> varList(dim); 
  
  if (dim>=1 && var1!=0) varList[0] = TString(var1); 
  if (dim>=2 && var2!=0) varList[1] = TString(var2); 
  if (dim>=3 && var3!=0) varList[2] = TString(var3); 
  
  return readTuple(tree, varList, maxEvents); 

}

Bool_t KernelDensity::readTuple(TTree* tree, const char* var1, const char* var2, 
                                UInt_t maxEvents) {

  UInt_t dim = phaseSpace()->dimensionality(); 
  std::vector<TString> varList(dim); 
  
  if (dim>=1 && var1!=0) varList[0] = TString(var1); 
  if (dim>=2 && var2!=0) varList[1] = TString(var2); 
  
  return readTuple(tree, varList, maxEvents); 

}

Bool_t KernelDensity::readTuple(TTree* tree, const char* var1, UInt_t maxEvents) {

  UInt_t dim = phaseSpace()->dimensionality(); 
  std::vector<TString> varList(dim); 
  
  if (dim>=1 && var1!=0) varList[0] = TString(var1); 
  
  return readTuple(tree, varList, maxEvents); 

}

UInt_t KernelDensity::numCells(void) {
  UInt_t dim = phaseSpace()->dimensionality(); 

  UInt_t i; 
  UInt_t cells = 1;
  for (i = 0; i<dim; i++) {
    Double_t lower = phaseSpace()->lowerLimit(i); 
    Double_t upper = phaseSpace()->upperLimit(i); 

    Int_t dimCells = TMath::Ceil( (upper-lower)/m_width[i] ); 

    cells *= dimCells; 
  }

  return cells; 
}

Int_t KernelDensity::cellIndex(std::vector<Double_t>& x) {

  UInt_t dim = phaseSpace()->dimensionality(); 

  UInt_t i; 
  UInt_t cell = 0;
  for (i = 0; i<dim; i++) {
    Double_t lower = phaseSpace()->lowerLimit(i); 
    Double_t upper = phaseSpace()->upperLimit(i); 

    Double_t xi = x[i];
    if (xi < lower || xi > upper) {
      return -1;
    }

    if (i<dim-1) { 
      Int_t dimCells = TMath::Ceil( (upper-lower)/m_width[i+1] ); 
      cell *= dimCells;
    }

    Int_t dimCells = TMath::Ceil( (upper-lower)/m_width[i] ); 
    Int_t dimCell  = TMath::Floor( (xi-lower)/m_width[i] ); 
    if (dimCell == dimCells) dimCell--;
    
    if (dimCell < 0 || dimCell >= dimCells) {
      printf("%20.20s ERROR: Something wrong! upper=%f, lower=%f, x=%f, width=%f, cell=%d, ncells=%d\n", m_name, 
             upper, lower, xi, m_width[i], dimCell, dimCells); 
      abort(); 
      return -1;
    }
//    if (dimCell == 20 && dim==1)
//      printf("DEBUG: CellIndex = %d, lower=%f, upper=%f, x=%f, width=%f\n", dimCell, lower, upper, xi, m_width[i]); 

    cell += dimCell; 
  }

  return cell; 

}

Bool_t KernelDensity::readTuple(TTree* tree, std::vector<TString> &vars, UInt_t maxEvents) {

  if (vars.size() != m_phaseSpace->dimensionality() ) {
    printf("%20.20s ERROR: Number of TTree variables (%d) in tree \"%s\" does not correspond to phase space dimensionality (%d)\n", 
           m_name, (UInt_t)vars.size(), tree->GetName(), m_phaseSpace->dimensionality() ); 
    abort(); 
  }  
  
  UInt_t nvars = vars.size(); 

  tree->ResetBranchAddresses();
  
  Long64_t nentries = tree->GetEntries();
  
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents; 
  
  Long64_t i;

  std::vector<Float_t> varArray(nvars); 

  UInt_t n; 
  for (n=0; n < nvars; n++) {
    printf("%20.20s INFO: Will read branch \"%s\" from tree \"%s\"\n", m_name, vars[n].Data(), tree->GetName()); 
    Int_t status = tree->SetBranchAddress(vars[n], &( varArray[n] ));
    if (status < 0) {
      printf("%20.20s WARNING: Error setting branch, status=%d\n", m_name, status); 
      abort(); 
    }
  }

  std::vector<Double_t> point(nvars); 

  UInt_t cells = numCells(); 

  m_dataVector.resize(cells); 
  
  Int_t cell; 
  for (cell = 0; cell < (Int_t)cells; cell++) {
    m_dataVector[cell].reserve(nentries/cells); 
  }

  UInt_t nout = 0;
  for(i=0; i<nentries; i++) {
//    printf("DEBUG: before: %f\n", varArray[0]);
    tree->GetEntry(i);
//    printf("DEBUG: after: %f\n", varArray[0]);
    
    for (n=0; n<nvars; n++) {
      point[n] = varArray[n]; 
    }
    
    if (!m_phaseSpace->withinLimits( point )) {
      nout ++; 
      printf("%20.20s WARNING: Ntuple point (", m_name); 
      for (n=0; n<nvars; n++) {
        printf("%f ", point[n]);
      }
      printf(", %f%%) outside phase space\n", 100.*float(nout)/float(i));
    } else {
      cell = cellIndex( point ); 
      if (cell>=0) m_dataVector[cell].push_back(point); 
    }
    
    if (i % 10000 == 0) {
      printf("%20.20s INFO: Read %lld/%lld events (%f%%)\n", m_name, i, nentries, 100.*float(i)/float(nentries));
    }
  }
  
  printf("%20.20s INFO: %lld events read in from \"%s\"\n", m_name, nentries-nout, tree->GetName() ); 

  return 1;
}

Double_t KernelDensity::rawDensity(std::vector<Double_t> &x, std::vector<TCell> &vector) {

//  UInt_t cells = numCells(); 

  UInt_t dim = m_phaseSpace->dimensionality(); 

  std::vector<UInt_t> iter(dim); 
  
  UInt_t j;
  for (j=0; j<dim; j++) {
    iter[j] = -1;
  }

  Double_t d = 0.;

  // Loop through the vertices of the N-dim 3^N cube
  do {
  
//    printf("DEBUG: Iterator "); 
//    for (j=0; j<dim; j++) {
//      printf(" %d", iter[j]); 
//    }
//    printf("\n"); 

    std::vector<Double_t> point(dim); 
    for (j=0; j<dim; j++) {
      point[j] = x[j] + iter[j]*m_width[j]; 
    }

    Int_t cell = cellIndex(point); 

    //// !!!!!! Check if the cell is unique! 

//    printf("DEBUG: Cell %d/%d\n", cell, cells); 

    if (cell >= 0) {

      UInt_t size = vector[cell].size(); 

      for (UInt_t i=0; i<size; i++) {
        Double_t sqsum = 0.;
        for (UInt_t var=0; var < dim; var++) {
          Double_t dx = (vector[cell][i][var] - x[var])/m_width[var];
          if (TMath::Abs(dx) < 1.) {
            sqsum += dx*dx;
          } else {
            sqsum = 1.;
            break;
          }
        }
        if (sqsum < 1.) d += (1.-sqsum); 
      }

    } // if cell >= 0

    // Increment iterator
    Bool_t run = 0; 
    for (j=0; j<dim; j++) {
      if (iter[j] <= 0) {
        iter[j]++; 
        run = 1; 
        break;
      } else {
        iter[j] = -1;
      }
    }
    if (!run) break;

  } while(1); 


  return d;
}

Double_t KernelDensity::density(std::vector<Double_t> &x) {

  Double_t rawData = rawDensity(x, m_dataVector);
  Double_t rawNorm = rawDensity(x, m_apprVector);
  
  if (rawNorm <= 0) {
//    printf("WARNING: Normalisation density <= 0!\n"); 
    return 0.;
  }

  Double_t corrEff = rawData/rawNorm; 
  if (m_approxDensity) {
    corrEff *= m_approxDensity->density(x); 
  }

//  printf("DEBUG: %f %f %f %f\n", x[0], rawData, rawNorm, corrEff);
  
  return corrEff;
}
