#include "PolynomialDensity.hh"
#include "AbsPhaseSpace.hh"
#include "OneDimPhaseSpace.hh"
#include "UniformDensity.hh"

#include "TMath.h"
#include "TVirtualFitter.h"
#include "TTree.h"

#include "Timer.hh"

static std::vector<std::vector<Double_t> > gData; 
static Int_t gPower; 
static AbsPhaseSpace* gPhsp; 
static std::vector<Double_t> gMiddle; 
static Double_t gMaximum; 
static std::vector<std::pair<Double_t, Double_t> > gIntegVector2D; 
static Int_t gIter;
static const char* gName; 

Double_t gPoly1D(Double_t* pars, Double_t x, Int_t pow) {
  // Number of parameters = pow
  Int_t psum, pnum; 
  Double_t sum = 1.;
  pnum = 0;
  for (psum=1; psum<=pow; psum++) {
    sum += pars[pnum]*TMath::Power(x - gMiddle[0], psum);
    pnum ++; 
  }

  return sum; 

}

Double_t gPoly1DInt(Double_t* pars, Int_t pow, AbsPhaseSpace* phsp) {
  Double_t low = phsp->lowerLimit(0); 
  Double_t up = phsp->upperLimit(0); 

  Int_t psum, pnum; 
  Double_t sum = up-low;
  
  pnum = 0;
  for (psum=1; psum<=pow; psum++) {
    sum += pars[pnum]*(TMath::Power(up - gMiddle[0], psum+1) - TMath::Power(low - gMiddle[0], psum+1))/(psum+1.);
    pnum ++; 
  }

  return sum; 

}

void gPoly1DLH(Int_t &, Double_t *, Double_t &f, Double_t *pars, Int_t) {
  std::vector<std::vector<Double_t> >::iterator i; 
  Double_t lhsum = 0.;
  Int_t error = 0; 
  gIter++;
  
  Double_t integ = gPoly1DInt(pars, gPower, gPhsp);
  for (i=gData.begin(); i != gData.end(); i++) {
    Double_t x = (*i)[0]; 
    Double_t pdf = gPoly1D(pars, x, gPower); 
    if (pdf < 0) {
      error = 1; 
    } else {
      lhsum += TMath::Log( pdf/integ ); 
    }
  }
  if (error) {
    f = gMaximum; 
  } else {
    f = -2.*lhsum; 
    if (f > gMaximum) gMaximum = f; 
  }

  if (timer(4)) {
    printf("%20.20s INFO: Iteration %d\n", gName, gIter); 
    for (Int_t v=0; v<gPower; v++)
      printf("%20.20s INFO:   P%d = %f\n", gName, v+1, pars[v]); 
  }
}

Double_t gPoly2D(Double_t* pars, Double_t x, Double_t y, Int_t pow) {
  // Number of parameters: 2 + 3 + ... +  (pow+1) = (pow+1)*pow/2-1

  Int_t psum, pnum; 
  Double_t sum = 1.;

  pnum = 0;
  for (psum=1; psum<=pow; psum++) {
    Int_t px; 
    for (px = 0; px <= psum; px++) {
      Int_t py = psum-px; 
      sum += pars[pnum]*TMath::Power(x - gMiddle[0], px)*TMath::Power(y - gMiddle[1], py); 
      pnum ++; 
    }
  }

  return sum; 
}

Double_t gPoly2DInt(Double_t* pars, Int_t pow, __attribute__((unused)) AbsPhaseSpace* phsp) {
  Double_t sum = 0; 

  std::vector<std::pair<Double_t, Double_t> >::iterator i;
  for (i = gIntegVector2D.begin(); i != gIntegVector2D.end(); i++) {
    Double_t x = (*i).first;
    Double_t y = (*i).second;

    Double_t pdf = gPoly2D(pars, x, y, pow);
    sum += pdf; 
  }
//  sum /= (Double_t) gIntegVector2D.size(); 

  return sum; 
}

void gPoly2DLH(Int_t &, Double_t *, Double_t &f, Double_t *pars, Int_t) {
  std::vector<std::vector<Double_t> >::iterator i; 
  Double_t lhsum = 0.;
  Int_t error = 0; 
  gIter++;
  
  Double_t integ = gPoly2DInt(pars, gPower, gPhsp);
  for (i=gData.begin(); i != gData.end(); i++) {
    Double_t x = (*i)[0];
    Double_t y = (*i)[1];
    Double_t pdf = gPoly2D(pars, x, y, gPower); 
    if (pdf < 0) {
      error = 1; 
    } else {
      lhsum += TMath::Log( pdf/integ ); 
    }
  }
  if (error) {
    f = gMaximum; 
  } else {
    f = -2.*lhsum; 
    if (f > gMaximum) gMaximum = f; 
  }

  if (timer(4)) {
    printf("%20.20s INFO: Iteration %d\n", gName, gIter); 
    Int_t pnum = 0;
    Int_t psum; 
    for (psum=1; psum<=(Int_t)gPower; psum++) {
      Int_t px; 
      for (px = 0; px <= psum; px++) {
        Int_t py = psum-px; 
        pnum++;
        printf("%20.20s INFO:   P%d_%d = %f\n", gName, px, py, pars[pnum-1]); 
      }
    }
  }
}

PolynomialDensity::PolynomialDensity(const char* pdfName, 
                      OneDimPhaseSpace* thePhaseSpace, 
                      UInt_t maxPower, 
                      TTree* tree, 
                      const char* var, 
                      UInt_t maxEvents) : AbsDensity(pdfName) {
  m_power = maxPower; 
  m_phaseSpace = thePhaseSpace; 
  m_dim = m_phaseSpace->dimensionality(); 
  gPhsp = thePhaseSpace; 
  gPower = maxPower; 
  gMaximum = -1.e50;
  gIter = 0;
  gName = pdfName; 

  if (m_dim != 1) {
    printf("%20.20s INFO: Dimensionality of the phase space is %d for constructor with 1 parameter\n", m_name, m_dim); 
    abort(); 
  }


  gMiddle.resize(1); 
  gMiddle[0] = (gPhsp->upperLimit(0) + gPhsp->lowerLimit(0))/2.;

  std::vector<TString> vars(1); 
  vars[0] = TString(var);
  readTuple(tree, vars, maxEvents); 

  // Create fitter and parameters
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, gPower);
  minuit->SetFCN(gPoly1DLH);

  double arglist[100];

  Double_t halfSize = (gPhsp->upperLimit(0) - gPhsp->lowerLimit(0))/2.; 
  for (Int_t v=1; v<=gPower; v++) {
    char parname[32]; 
    sprintf(parname, "P%d", v); 
    minuit->SetParameter(v-1, parname, 0., 0.01*TMath::Power(halfSize, 1./(Double_t)v), 0, 0);
  }

  set_timer(); 
  
  // Fit polynomial
  arglist[0] = 1;
  minuit->ExecuteCommand("MIGRAD", arglist, 0);

  for (Int_t v=0; v<gPower; v++) {
    m_par[v] = minuit->GetParameter(v); 
    printf("%20.20s INFO: Parameter %d = %f\n", m_name, v, m_par[v]); 
  }

}

PolynomialDensity::PolynomialDensity(const char* pdfName, 
                      AbsPhaseSpace* thePhaseSpace, 
                      UInt_t maxPower, 
                      TTree* tree, 
                      const char* var1, 
                      const char* var2, 
                      UInt_t integEvents, 
                      UInt_t maxEvents) : AbsDensity(pdfName) {
  gPhsp = thePhaseSpace; 
  gPower = maxPower; 
  m_power = maxPower; 
  m_phaseSpace = thePhaseSpace; 
  m_dim = m_phaseSpace->dimensionality(); 
  gMaximum = -1.e50;
  gIter = 0;
  gName = pdfName; 
  
  if (m_dim != 2) {
    printf("%20.20s INFO: Dimensionality of the phase space is %d for constructor with 2 parameters\n", m_name, m_dim); 
    abort(); 
  }

  gMiddle.resize(2); 
  gMiddle[0] = (gPhsp->upperLimit(0) + gPhsp->lowerLimit(0))/2.;
  gMiddle[1] = (gPhsp->upperLimit(1) + gPhsp->lowerLimit(1))/2.;
  
  std::vector<TString> vars(2); 
  vars[0] = TString(var1);
  vars[1] = TString(var2);
  readTuple(tree, vars, maxEvents); 
  
  UInt_t i; 
  
  UniformDensity uniform(pdfName, m_phaseSpace); 
  std::vector<Double_t> genVector(2); 
  gIntegVector2D.clear(); 
  for (i=0; i<integEvents; i++) {
    uniform.generate(genVector); 
    Double_t x = genVector[0]; 
    Double_t y = genVector[1]; 
    std::pair<Double_t, Double_t> point(x, y); 
    gIntegVector2D.push_back(point);
  } 
  
  // Create fitter and parameters
  
  Int_t npars = (gPower+2)*(gPower+1)/2 - 1; 
  
  TVirtualFitter *minuit = TVirtualFitter::Fitter(0, npars);
  minuit->SetFCN(gPoly2DLH);
  
  double arglist[100];

  Int_t psum, pnum;
  pnum = 0;
  Double_t halfSizeX = (gPhsp->upperLimit(0) - gPhsp->lowerLimit(0))/2.; 
  Double_t halfSizeY = (gPhsp->upperLimit(1) - gPhsp->lowerLimit(1))/2.; 

  for (psum=1; psum<=(Int_t)m_power; psum++) {
    Int_t px; 
    for (px = 0; px <= psum; px++) {
      Int_t py = psum-px; 
      pnum ++; 
      char parname[32]; 
      sprintf(parname, "P%d_%d", px, py); 
      Double_t err = 0.01;
      if (px > 0) err /= TMath::Power(halfSizeX, 1./(Double_t)px);
      if (py > 0) err /= TMath::Power(halfSizeY, 1./(Double_t)py);
      minuit->SetParameter(pnum-1, parname, 0., err, 0, 0);
    }
  }

  set_timer(); 

  // Fit polynomial
  arglist[0] = 1;
  minuit->ExecuteCommand("MIGRAD", arglist, 0);

  for (Int_t v=0; v<npars; v++) {
    m_par[v] = minuit->GetParameter(v); 
    printf("%20.20s INFO: Parameter %d = %f\n", m_name, v, m_par[v]); 
  }

}

PolynomialDensity::~PolynomialDensity() {

}

void PolynomialDensity::readTuple(TTree* tree, std::vector<TString> &vars, UInt_t maxEvents) {

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

  gData.clear();
  gData.reserve(nentries); 
  
  set_timer(); 

  UInt_t nout = 0;
  for(i=0; i<nentries; i++) {
    tree->GetEntry(i);
    
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
      for (UInt_t v = 0; v < nvars; v++) {
        gData.push_back(point); 
      }
    }
    
    if (i % 100 == 0 && timer(2)) {
      printf("%20.20s INFO: Read %lld/%lld events (%f%%)\n", m_name, i, nentries, 100.*float(i)/float(nentries));
    }
  }
  
  printf("%20.20s INFO: %lld events read in from \"%s\"\n", m_name, nentries-nout, tree->GetName() ); 

}

Double_t PolynomialDensity::density(std::vector<Double_t> &x) {
  gPhsp = m_phaseSpace; 
  gPower = m_power; 
  gMaximum = -1.e50;
  if (m_dim == 1) {
    gMiddle.resize(1); 
    gMiddle[0] = (gPhsp->upperLimit(0) + gPhsp->lowerLimit(0))/2.;
    return gPoly1D(m_par, x[0], m_power); 
  } else if (m_dim == 2) {
    gMiddle.resize(2); 
    gMiddle[0] = (gPhsp->upperLimit(0) + gPhsp->lowerLimit(0))/2.;
    gMiddle[1] = (gPhsp->upperLimit(1) + gPhsp->lowerLimit(1))/2.;
    return gPoly2D(m_par, x[0], x[1], m_power); 
  }
  return 0.;
}
