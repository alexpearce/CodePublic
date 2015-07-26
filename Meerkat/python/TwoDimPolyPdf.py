# This is the example of using 2D polynomial fit class. 
#
#  Generate a distribution in the circular phase space x**2+y**2<1
#  and use PolynomialDensity class to fit it. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import OneDimPhaseSpace, ParametricPhaseSpace, FormulaDensity, PolynomialDensity
from ROOT import TFile, TNtuple, TCanvas, TH2F, TText

# First create the 1D phase space for variable x
xphsp = OneDimPhaseSpace("PhspX", -1., 1.) 
  
# Now create the parametric phase space for (x,y) where limits on variable y are functions of x
phsp = ParametricPhaseSpace("PhspParam", xphsp, "-sqrt(1-x^2)", "sqrt(1-x^2)", -1., 1.)

outfile = TFile("TwoDimPolyPdf.root", "RECREATE")
ntuple = TNtuple("ntuple", "2D NTuple", "x:y") 

# True PDF
true_pdf = FormulaDensity("TruePDF", phsp, "1.+0.1*x-0.4*x^2-0.8*y^2") 

# Generate 10000 toys according to the true PDF
true_pdf.generate(ntuple, 50000) 

poly = PolynomialDensity("PolyPDF", phsp, 4, ntuple, "x", "y", 100000) 

poly_hist = TH2F("poly", "Polynomial PDF", 200, -1.1, 1.1, 200, -1.1, 1.1) 

poly.project(poly_hist) 
poly_hist.Write()

gStyle.SetOptStat(0)
  
canvas = TCanvas("canvas", "TwoDimPolyPdf", 400, 400) 
poly_hist.Draw("zcol") 

canvas.Print("TwoDimPolyPdf.png") 

#outfile.Close() 

