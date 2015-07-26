#  This is an example of using adaptive kernel class for estimation of 1D PDF 
#  with a narrow peak. 
#
#  It also uses FormulaDensity to describe the true PDF in a parametrised way. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import OneDimPhaseSpace, FormulaDensity
from ROOT import BinnedKernelDensity, BinnedDensity, AdaptiveKernelDensity
from ROOT import TFile, TNtuple, TCanvas, TH1F, TText

# Define 1D phase space for a variable in range (-1, 1)
phsp = OneDimPhaseSpace("Phsp1D", -1., 1.)

# Create a "true" PDF over this phase space in a parametrised way
formula = FormulaDensity("FormulaPDF", phsp, "(1. + 4.*exp(-x*x/2/0.1/0.1))*( abs(x)<1 )")

outfile = TFile("OneDimAdaptiveKernel.root", "RECREATE")
ntuple = TNtuple("ntuple", "1D NTuple", "x")

# Generate toys according to the "true" PDF
formula.generate(ntuple, 10000)
  
# Create kernel PDF from the generate sample
kde = BinnedKernelDensity("KernelPDF", 
                        phsp,   # Phase space
                        ntuple, # Input ntuple
                        "x",     # Variable to use
                        1000,    # Number of bins
                        0.1,     # Kernel width
                        0,       # Approximation PDF (0 for flat approximation)
                        100000   # Sample size for MC convolution (0 for binned convolution)
                       )

# Create adaptive kernel PDF with the kernel width depending on the binned PDF 
kde_adaptive = AdaptiveKernelDensity("AdaptivePDF", 
                          phsp,        # Phase space
                          ntuple,      # Input ntuple
                          "x",         # Variable to use
                          1000,        # Number of bins
                          0.1,         # Initial kernel width
                          kde,         # PDF for width scaling
                          0,           # Approximation PDF (0 for flat approximation)
                          100000       # Sample size for MC convolution (0 for binned convolution)
                         )

# Write binned PDF into a file
kde_adaptive.writeToFile("OneDimAdaptiveKernel_bins.root")
  
true_hist = TH1F("true_pdf", "PDF", 200, -1.1, 1.1)
kernel_hist = TH1F("kernel_pdf", "Kernel PDF", 200, -1.1, 1.1)
adaptive_hist = TH1F("adaptive_pdf", "Adaptive kernel PDF", 200, -1.1, 1.1)

kernel_diff = TH1F("kernel_diff", "PDF bias", 200, -1.1, 1.1)
adaptive_diff = TH1F("adaptive_diff", "PDF bias", 200, -1.1, 1.1)

# Fill histogram with true PDF
formula.project(true_hist); 
true_hist.Write();

# Fill histogram with the result of kernel density estimation
kde.project(kernel_hist); 
kernel_hist.Write();

kde_adaptive.project(adaptive_hist)
adaptive_hist.Write()

gStyle.SetOptStat(0)

canvas = TCanvas("canvas", "OneDimAdaptiveKernel", 400, 600)
canvas.Divide(1, 2)
canvas.cd(1)
true_hist.SetLineColor(1)
true_hist.Draw()
kernel_hist.Scale( true_hist.GetSumOfWeights() / kernel_hist.GetSumOfWeights() )
kernel_hist.SetLineColor(2)
kernel_hist.Draw("same")
adaptive_hist.Scale( true_hist.GetSumOfWeights() / adaptive_hist.GetSumOfWeights() )
adaptive_hist.SetLineColor(4)
adaptive_hist.Draw("same")

t = TText()
t.SetTextSize(0.04)
t.SetTextColor(1)
t.DrawText(0.2, 4.5, "True PDF")
t.SetTextColor(2)
t.DrawText(0.2, 4.2, "Fixed kernel PDF")
t.SetTextColor(4)
t.DrawText(0.2, 3.9, "Adaptive kernel PDF")

kernel_diff.Add(kernel_hist, true_hist, 1., -1.)
adaptive_diff.Add(adaptive_hist, true_hist, 1., -1.)
canvas.cd(2)
kernel_diff.SetMinimum(-1.)
kernel_diff.SetMaximum( 1.)
kernel_diff.SetLineColor(2)
kernel_diff.Draw()
adaptive_diff.SetLineColor(4)
adaptive_diff.Draw("same")

t.SetTextColor(2)
t.DrawText(0.2, 0.8, "Fixed kernel PDF")
t.SetTextColor(4)
t.DrawText(0.2, 0.6, "Adaptive kernel PDF")

canvas.Print("OneDimAdaptiveKernel.png")

#outfile.Close()
