#  An example of using 2D Dalitz plot phase space class. 
#  Generate points across the Dalitz plot and use kernel 
#  estimation to estimate efficiency shape. 
#
#  Demonstrates the use of 2D projection and 1D slices. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import DalitzPhaseSpace, BinnedKernelDensity, FactorisedDensity, FormulaDensity
from ROOT import BinnedDensity, AdaptiveKernelDensity, PolynomialDensity
from ROOT import TFile, TNtuple, TCanvas, TH1F, TH2F, TText

# Define masses of initial ("d") and final state particles ("a,b,c") for the Dalitz phase space. 
md = 1.8646  # D0
ma = 0.497   # KS0
mb = 0.139   # Pi
mc = 0.139   # Pi

# Define Dalitz phase space for D0->KsPiPi decay
# The phase space is defined in x=m^2(ab), y=m^2(ac) variables
phsp = DalitzPhaseSpace("PhspDalitz", md, ma, mb, mc)

# Create a true PDF over this phase space, polynomial of power 4 in x and y
truepdf = FormulaDensity("TruePDF", phsp, "sqrt(1.-0.1*(x-1.3)^4)*(1.+2.*exp(-y))")

outfile = TFile.Open("DalitzPdf.root", "RECREATE")
ntuple = TNtuple("ntuple", "2D NTuple", "x:y")

# Generate 50000 toys according to the "true" PDF
truepdf.generate(ntuple, 25000)

# Create an approximation PDF, polynomial of power 3
poly = PolynomialDensity("PolyPDF", 
                         phsp,    # Phase space
                         2,       # Power of the polynomial
                         ntuple,  # Input ntuple
                         "x","y", # Ntuple variables
                         50000    # Sample size for MC integration
                         )

# Create kernel PDF from the generated sample
kde = BinnedKernelDensity("KernelPDF", 
                          phsp,    # Phase space
                          ntuple,  # Input ntuple
                          "x","y", # Variables to use
                          200,200, # Numbers of bins
                          0.4, 0.4,# Kernel widths
                          poly,    # Approximation PDF (0 for flat approximation)
                          50000    # Sample size for MC convolution (0 for binned convolution)
                         )

# Write binned PDF into a file
kde.writeToFile("DalitzPdf_bins.root")

# That's it. Now fill some histograms and show the results. 

true_hist = TH2F("true", "True PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                     100, phsp.lowerLimit(1), phsp.upperLimit(1))
poly_hist = TH2F("poly", "Polynomial PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                     100, phsp.lowerLimit(1), phsp.upperLimit(1))
kernel_hist = TH2F("kernel", "Kernel PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                           100, phsp.lowerLimit(1), phsp.upperLimit(1))

true_kpi = TH1F("true_kpi", "M^{2}(K_{S}#pi) slice", 100, phsp.lowerLimit(0), phsp.upperLimit(0) )
poly_kpi = TH1F("poly_kpi", "M^{2}(K_{S}#pi) slice", 100, phsp.lowerLimit(0), phsp.upperLimit(0) )
kernel_kpi = TH1F("kernel_kpi", "M^{2}(K_{S}#pi) slice", 100, phsp.lowerLimit(0), phsp.upperLimit(0) )
true_pipi = TH1F("true_pipi", "M^{2}(#pi#pi) slice", 100, phsp.lowerLimit(1), phsp.upperLimit(1) )
poly_pipi = TH1F("poly_pipi", "M^{2}(#pi#pi) slice", 100, phsp.lowerLimit(1), phsp.upperLimit(1) )
kernel_pipi = TH1F("kernel_pipi", "M^{2}(#pi#pi) slice", 100, phsp.lowerLimit(1), phsp.upperLimit(1) )

# Coordinates in the phase space where we will create the 1D slices.
from ROOT import std, Double
v = std.vector(Double)(2)
v[0] = 1.0;
v[1] = 0.3;

# Fill histogram with the result of kernel density estimation
truepdf.project(true_hist)
true_hist.Write()
poly.project(poly_hist)
poly_hist.Write()
kde.project(kernel_hist)
kernel_hist.Write()

# Fill 1D slice in x projection, at the point y=v[1]
truepdf.slice(v, 0, true_kpi) 
true_kpi.Write() 
poly.slice(v, 0, poly_kpi)
poly_kpi.Write()
kde.slice(v, 0, kernel_kpi)
kernel_kpi.Write()

# Fill 1D slice in y projection, at the point x=v[0]
truepdf.slice(v, 1, true_pipi)
true_pipi.Write()
poly.slice(v, 1, poly_pipi)
poly_pipi.Write()
kde.slice(v, 1, kernel_pipi)
kernel_pipi.Write()

# Plot 2D and 1D histograms

gStyle.SetOptStat(0)

canvas = TCanvas("canvas", "DalitzPdf", 600, 600)
canvas.Divide(2,2)

#canvas.cd(1)
#ntuple.Draw("y:x", "", "")

canvas.cd(1)
true_hist.Draw("zcol")
true_hist.GetXaxis().SetTitle("M^{2}(K_{S}#pi), GeV")
true_hist.GetYaxis().SetTitle("M^{2}(#pi#pi), GeV")

canvas.cd(2)
poly_hist.Scale( true_hist.GetSumOfWeights() / poly_hist.GetSumOfWeights() )
poly_hist.Draw("zcol")
poly_hist.GetXaxis().SetTitle("M^{2}(K_{S}#pi), GeV")
poly_hist.GetYaxis().SetTitle("M^{2}(#pi#pi), GeV")

canvas.cd(3)
kernel_hist.Scale( true_hist.GetSumOfWeights() / kernel_hist.GetSumOfWeights() )
kernel_hist.Draw("zcol")
kernel_hist.GetXaxis().SetTitle("M^{2}(K_{S}#pi), GeV")
kernel_hist.GetYaxis().SetTitle("M^{2}(#pi#pi), GeV")

canvas.cd(4)
kernel_kpi.Scale( true_kpi.GetSumOfWeights() / kernel_kpi.GetSumOfWeights() )
poly_kpi.Scale( true_kpi.GetSumOfWeights() / poly_kpi.GetSumOfWeights() )
true_kpi.Draw()
kernel_kpi.SetLineColor(2)
kernel_kpi.Draw("same")
poly_kpi.SetLineColor(3)
poly_kpi.Draw("same")

#canvas.cd(5)
#kernel_pipi.Scale( true_pipi.GetSumOfWeights() / kernel_pipi.GetSumOfWeights() )
#poly_pipi.Scale( true_pipi.GetSumOfWeights() / poly_pipi.GetSumOfWeights() )
#true_pipi.Draw()
#kernel_pipi.SetLineColor(2)
#kernel_pipi.Draw("same")
#poly_pipi.SetLineColor(3)
#poly_pipi.Draw("same")

canvas.Print("DalitzPdf.png")

#outfile.Close()
