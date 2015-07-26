#   This example demonstrated the use of combined phase spaces and factorised density. 
#   Here we create a 2D phase space for two variables in range (-1,1), and estimate 
#   the density which has a slight correlation between the two variables. 
#
#   We compare three approaches to estimate the true PDF: 
#     * Factorised PDF, which is simply a product of kernel PDFs in each variable
#     * Kernel PDF with the flat approximation
#     * Kernel PDF with factorised approximation. 
#   The latter gives better results at the edges of the phase space. 

from ROOT import gSystem, gStyle

gSystem.Load("../lib/libMeerkat.so")

from ROOT import OneDimPhaseSpace, CombinedPhaseSpace, BinnedKernelDensity, FactorisedDensity, FormulaDensity
from ROOT import TFile, TNtuple, TCanvas, TH2F, TText

# Define 1D phase spaces for each of the two variables
phsp_x = OneDimPhaseSpace("PhspX", -1., 1.)
phsp_y = OneDimPhaseSpace("PhspY", -1., 1.)

# Define combined phase space for the two variables (AKA direct product)
phsp = CombinedPhaseSpace("PhspCombined", phsp_x, phsp_y)

outfile = TFile("CombinedPdf.root", "RECREATE")
ntuple = TNtuple("ntuple", "2D NTuple", "x:y")

truepdf = FormulaDensity("TruePDF", phsp, "1.-0.05*(x+y)^4-0.2*(x-y)^2")
truepdf.generate(ntuple, 50000)

# Create kernel PDF from the generated sample for projection to variable x
kde_x = BinnedKernelDensity("KernelPDF_X",
                          phsp_x,  # Phase space
                          ntuple,  # Input ntuple
                          "x",     # Variables to use
                          200,     # Numbers of bins
                          0.2,     # Kernel width
                          0,       # Approximation PDF (0 for flat approximation)
                          100000   # Sample size for MC convolution (0 for binned convolution)
                         )

# Create kernel PDF from the generated sample for projection to variable y
kde_y = BinnedKernelDensity("KernelPDF_Y",
                          phsp_y,  # Phase space
                          ntuple,  # Input ntuple
                          "y",     # Variables to use
                          200,     # Numbers of bins
                          0.2,     # Kernel width
                          0,       # Approximation PDF (0 for flat approximation)
                          100000   # Sample size for MC convolution (0 for binned convolution)
                         )

# Create factorised PDF which is the product of PDFs in x and y
fact = FactorisedDensity("FactPDF", phsp, kde_x, kde_y)

# Create kernel PDF from the generated sample, using flat approximation
kde = BinnedKernelDensity("KernelPDF",
                          phsp,     # Phase space
                          ntuple,   # Input ntuple
                          "x","y",  # Variables to use
                          100,100,  # Numbers of bins
                          0.4, 0.4, # Kernel widths
                          0,        # Approximation PDF (0 for flat approximation)
                          100000    # Sample size for MC convolution (0 for binned convolution)
                         )

# Create kernel PDF from the generated sample, using factorised PDF as approximation
kde_factappr = BinnedKernelDensity("KernelPDFWithFactApprox",
                          phsp,     # Phase space
                          ntuple,   # Input ntuple
                          "x","y",  # Variables to use
                          100,100,  # Numbers of bins
                          0.4, 0.4, # Kernel widths
                          fact,     # Approximation PDF (0 for flat approximation)
                          100000    # Sample size for MC convolution (0 for binned convolution)
                         )

# That's it. Now fill some histograms and show the results. 

true_hist = TH2F("true", "True PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                     100, phsp.lowerLimit(1), phsp.upperLimit(1))
fact_hist = TH2F("fact", "Factorised PDF", 100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                           100, phsp.lowerLimit(1), phsp.upperLimit(1))
kernel_hist = TH2F("kernel", "Kernel PDF with flat approximation", 
                                           100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                           100, phsp.lowerLimit(1), phsp.upperLimit(1))
kernel_factappr_hist = TH2F("kernel_fa", "Kernel PDF with factorised approximation", 
                                           100, phsp.lowerLimit(0), phsp.upperLimit(0), 
                                           100, phsp.lowerLimit(1), phsp.upperLimit(1))

truepdf.project(true_hist)
true_hist.Write()

# Fill histogram with the result of kernel density estimation
kde.project(kernel_hist)
kernel_hist.Write()

kde_factappr.project(kernel_factappr_hist)
kernel_factappr_hist.Write()

fact.project(fact_hist)
fact_hist.Write()

gStyle.SetOptStat(0)

canvas = TCanvas("canvas", "CombinedPdf", 500, 500)
canvas.Divide(2,2)

canvas.cd(1)
true_hist.SetMaximum(1.1)
true_hist.SetMinimum(0.0)
true_hist.Draw("zcol")
#ntuple.Draw("y:x", "", "same")
canvas.cd(2)
fact_hist.Scale( true_hist.GetSumOfWeights() / fact_hist.GetSumOfWeights() )
fact_hist.SetMaximum(1.1)
fact_hist.SetMinimum(0.0)
fact_hist.Draw("zcol")
canvas.cd(3)
kernel_hist.Scale( true_hist.GetSumOfWeights() / kernel_hist.GetSumOfWeights() )
kernel_hist.SetMaximum(1.1)
kernel_hist.SetMinimum(0.0)
kernel_hist.Draw("zcol")
canvas.cd(4)
kernel_factappr_hist.Scale( true_hist.GetSumOfWeights() / kernel_factappr_hist.GetSumOfWeights() )
kernel_factappr_hist.SetMaximum(1.1)
kernel_factappr_hist.SetMinimum(0.0)
kernel_factappr_hist.Draw("zcol")

canvas.Print("CombinedPdf.png")

#outfile.Close()
