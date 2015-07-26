#  This is the simplest possible example of using Meerkat package. 
#
#  Generate a uniform distribution in the range (-1,1)
#  and apply a Meerkat procedure to estimate PDF from it. 
#
#  The application creates the ROOT file with the generated ntuple, 
#  true and estimated distributions. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import OneDimPhaseSpace, UniformDensity, BinnedKernelDensity
from ROOT import TFile, TNtuple, TCanvas, TH1F

# Define 1D phase space for a variable in range (-1, 1)
phsp = OneDimPhaseSpace("Phsp1D", -1., 1.)

# Create a uniform PDF over this phase space
uniform = UniformDensity("UniformPDF", phsp)

outfile = TFile("OneDimPdf.root", "RECREATE")
ntuple = TNtuple("ntuple", "1D NTuple", "x")

# Generate 10000 toys according to the uniform PDF
uniform.generate(ntuple, 10000)

# Create kernel PDF from the generate sample
kde = BinnedKernelDensity("KernelPDF", 
                          phsp, 
                          ntuple, # Input ntuple
                          "x",    # Variable to use
                          1000,   # Number of bins
                          0.2,    # Kernel width
                          0,      # Approximation PDF (0 for flat approximation)
                          100000  # Sample size for MC convolution (0 for binned convolution)
                         )

# Write binned PDF into a file
kde.writeToFile("OneDimPdfBins.root")

uniform_hist = TH1F("unform", "PDF", 200, -1.5, 1.5)
kernel_hist = TH1F("kernel", "Kernel PDF", 200, -1.5, 1.5)

uniform.project(uniform_hist) 
kde.project(kernel_hist)

gStyle.SetOptStat(0)

canvas = TCanvas("canvas", "OneDimPdf", 400, 400)

uniform_hist.Draw()
kernel_hist.Scale( uniform_hist.GetSumOfWeights() / kernel_hist.GetSumOfWeights() )
kernel_hist.SetLineColor(2)
kernel_hist.Draw("same")
uniform_hist.GetXaxis().SetTitle("x")

canvas.Update()
canvas.Print("OneDimPdf.png")

#outfile.Close()
