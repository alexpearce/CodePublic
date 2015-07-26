#  An example of using 2D parametric phase space class. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import OneDimPhaseSpace, ParametricPhaseSpace, BinnedKernelDensity, FormulaDensity
from ROOT import TFile, TNtuple, TCanvas, TH1F, TH2F, TText

# First create the 1D phase space for variable x
xphsp = OneDimPhaseSpace("PhspX", -1., 1.)
  
# Now create the parametric phase space for (x,y) where limits on variable y are functions of x
phsp = ParametricPhaseSpace("PhspParam", xphsp, "-sqrt(1-x^2)", "sqrt(1-x^2)", -1., 1.)

# Create a uniform PDF over this phase space
# Sum of a polynomial shape and a Gaussian peak with correlation between x and y
truepdf = FormulaDensity("TruePDF", phsp, "1.-0.8*x^2-0.8*y^2+exp(-(x-y)^2/2./(0.4)^2-(x+y)^2/2./(0.8)^2)") 

# Create an approximation PDF over this phase space
# Just a polymonial shape. It will mostly be used to approximate the PDF 
# at the edges of the phase space.
approxpdf = FormulaDensity("TruePDF", phsp, "1.-0.8*x^2-0.8*y^2") 

outfile = TFile("ParamPhsp.root", "RECREATE")
ntuple = TNtuple("ntuple", "2D NTuple", "x:y") 

# Generate 50000 toys according to the "true" PDF
truepdf.generate(ntuple, 50000) 

# Create kernel PDF from the generated sample. Use polynomial shape as an approximation PDF 
kde = BinnedKernelDensity("KernelPDF", 
                          phsp,      # Phase space
                          ntuple,    # Input ntuple
                          "x","y",   # Variables to use
                          200,200,   # Numbers of bins
                          0.2, 0.2,  # Kernel widths
                          approxpdf, # Approximation PDF (0 for flat approximation)
                          100000     # Sample size for MC convolution (0 for binned convolution)
                         ) 

# Write binned PDF into a file
kde.writeToFile("ParamPhsp_bins.root") 

# That's it. Now fill some histograms and show the results. 

true_hist = TH2F("true", "True PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                     100, phsp.lowerLimit(1), phsp.upperLimit(1)) 
kernel_hist = TH2F("kernel", "Kernel PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                           100, phsp.lowerLimit(1), phsp.upperLimit(1)) 

true_kpi = TH1F("true_kpi", "X slice", 100, phsp.lowerLimit(0), phsp.upperLimit(0) ) 
kernel_kpi = TH1F("kernel_kpi", "X slice", 100, phsp.lowerLimit(0), phsp.upperLimit(0) ) 
true_pipi = TH1F("true_pipi", "Y slice", 100, phsp.lowerLimit(1), phsp.upperLimit(1) ) 
kernel_pipi = TH1F("kernel_pipi", "Y slice", 100, phsp.lowerLimit(1), phsp.upperLimit(1) ) 

# Coordinates in the phase space where we will create the 1D slices.
from ROOT import std, Double
v = std.vector(Double)(2) 
v[0] = 0.
v[1] = 0.

# Fill histogram with the result of kernel density estimation
# In this particular case, since we are plotting a 2D slice of a 2D distribution, 
# the value of "v" is ignored. 
truepdf.slice(v, 0, 1, true_hist) 
true_hist.Write()
kde.slice(v, 0, 1, kernel_hist) 
kernel_hist.Write()

# Fill 1D slice in x projection, at the point y=v[1]
truepdf.slice(v, 0, true_kpi) 
true_kpi.Write() 
kde.slice(v, 0, kernel_kpi) 
kernel_kpi.Write() 

# Fill 1D slice in y projection, at the point x=v[0]
truepdf.slice(v, 1, true_pipi) 
true_pipi.Write() 
kde.slice(v, 1, kernel_pipi) 
kernel_pipi.Write() 

# Plot 2D and 1D histograms

gStyle.SetOptStat(0)

canvas = TCanvas("canvas", "ParamPhsp", 500, 500) 
canvas.Divide(2,2) 

canvas.cd(1) 
true_hist.Draw("zcol") 
#ntuple.Draw("y:x", "", "same") 
canvas.cd(2)
kernel_hist.Scale( true_hist.GetSumOfWeights() / kernel_hist.GetSumOfWeights() ) 
kernel_hist.Draw("zcol") 
canvas.cd(3) 
kernel_kpi.Scale( true_kpi.GetSumOfWeights() / kernel_kpi.GetSumOfWeights() ) 
true_kpi.Draw() 
kernel_kpi.SetLineColor(2) 
kernel_kpi.Draw("same") 
canvas.cd(4) 
kernel_pipi.Scale( true_pipi.GetSumOfWeights() / kernel_pipi.GetSumOfWeights() ) 
true_pipi.Draw() 
kernel_pipi.SetLineColor(2) 
kernel_pipi.Draw("same") 

canvas.Print("ParamPhsp.png") 

#outfile.Close() 
