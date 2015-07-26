#  This example demonstrates the use of weighted ntuple. 
#
#  Generate a uniform distribution of variable x in the range (-1,1) 
#  with the weight which is a function of x, 
#  and apply a Meerkat procedure to estimate PDF from it. 
#
#  The application creates the ROOT file with the generated ntuple, 
#  true and estimated distributions. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import OneDimPhaseSpace, BinnedKernelDensity
from ROOT import TFile, TNtuple, TCanvas, TH1F, TText, TRandom3

# Define 1D phase space for a variable in range (-1, 1)
phsp = OneDimPhaseSpace("Phsp1D", -1., 1.)

outfile = TFile("WeightedTuple.root", "RECREATE")
ntuple = TNtuple("ntuple", "1D NTuple with weight", "x:w")
  
# Generate 10000 toys according to the uniform PDF
# with weight=x**2
uniform_hist = TH1F("unform", "PDF", 200, -1.5, 1.5)
rnd = TRandom3()
for i in range(0, 10000) :
  x = rnd.Rndm()*2.-1.
  w = x*x 
  ntuple.Fill(x, w) 
  uniform_hist.Fill(x, w)
ntuple.Write() 

# Create kernel PDF from the generated sample
# If the number of variables passed to the constructor is larger by 1 than the phase 
# space dimensionality, the last variables is considered as the weight
kde = BinnedKernelDensity("KernelPDF",
                          phsp, 
                          ntuple, # Input ntuple
                          "x",     # Variable to use
                          "w",     # Weight variable
                          1000,    # Number of bins
                          0.1,     # Kernel width
                          0,       # Approximation PDF (0 for flat approximation)
                          100000   # Sample size for MC convolution (0 for binned convolution)
                         ) 

# Write binned PDF into a file
kde.writeToFile("WeightedTuple_bins.root") 

# That's it. Now fill some control histograms and plot them. 
kernel_hist = TH1F("kernel", "Kernel PDF", 200, -1.5, 1.5) 

kde.project(kernel_hist) 

uniform_hist.Write()
kernel_hist.Write()

gStyle.SetOptStat(0)
  
canvas = TCanvas("canvas", "WeightedTuple", 400, 400) 
uniform_hist.Draw() 
kernel_hist.Scale( uniform_hist.GetSumOfWeights() / kernel_hist.GetSumOfWeights() ) 
kernel_hist.SetLineColor(2) 
kernel_hist.Draw("same") 

canvas.Print("WeightedTuple.png") 

#outfile.Close() 
