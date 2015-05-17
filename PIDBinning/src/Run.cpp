#include "EfficiencyBin.h"
#include "EfficiencyBinCollection.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

#include "TH2Poly.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"

int main() {
  TStyle* atlasStyle = new TStyle("ATLAS", "Atlas style");
  atlasStyle->SetHistMinimumZero();

  Int_t icol = 0;  // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  Double_t red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
  Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
  Double_t blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  atlasStyle->SetNumberContours(NCont);
  atlasStyle->SetPadTopMargin(0.06);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.16);
  atlasStyle->SetTitleXOffset(1.4);
  atlasStyle->SetTitleYOffset(1.4);
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2.);
  atlasStyle->SetLineStyleString(2, "[12 12]");  // postscript dashes
  atlasStyle->SetEndErrorSize(0.);
  atlasStyle->SetOptTitle(0);
  atlasStyle->SetOptStat(0);
  atlasStyle->SetOptFit(0);
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  // Let's create a collection object
  auto col = new EfficiencyBinCollection();

  // Put in a few squares
  Int_t nbins = 50;
  for (int iX = 0; iX < nbins; iX++) {
    for (int iY = 0; iY < nbins; iY++) {
      Double_t eta_low = 2. + (5. - 2.)/nbins * iY;
      Double_t eta_high = 2. + (5. - 2.)/nbins * (iY + 1);
      Double_t p_low = 3000. + (100000. - 3000.)/nbins * iX;
      Double_t p_high = 3000. + (100000. - 3000.)/nbins * (iX + 1);
      col->AddBin(EfficiencyBin::RectangularBin(p_low, eta_low, p_high, eta_high));
    }
  }

  // for (int iX = 0; iX < 50; iX++) {
  // for (int iY = 0; iY < 50; iY++) {
  // col->AddBin(EfficiencyBin::SquareBin(0.1 * iX, 0.1 * iY, 0.1));
  //}
  //}

  auto file = TFile::Open("./out_K.root");
  auto tree = dynamic_cast<TTree*>(file->Get("RSDStCalib"));
  Double_t P = 0.;
  Double_t ETA = 0.;
  Double_t PID = 0.;
  Double_t weight = 0.;
  Double_t run_number = 0.;
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("K_P", 1);
  tree->SetBranchStatus("K_Eta", 1);
  tree->SetBranchStatus("K_CombDLLK", 1);
  tree->SetBranchStatus("runNumber", 1);
  tree->SetBranchStatus("nsig_sw", 1);
  tree->SetBranchAddress("K_P", &P);
  tree->SetBranchAddress("K_Eta", &ETA);
  tree->SetBranchAddress("K_CombDLLK", &PID);
  tree->SetBranchAddress("runNumber", &run_number);
  tree->SetBranchAddress("nsig_sw", &weight);

  for (int i = 0; i < tree->GetEntries(); i++) {
    if (i % 10000 == 0) {
      std::cout << "Loaded " << i << " / " << tree->GetEntries() << std::endl;
    }
    tree->GetEntry(i);
    //if(run_number<126915|| run_number>126933){continue;}
    //if(run_number<125951|| run_number>125951){continue;}

    col->Fill(P, ETA, PID > 5., weight);
  }

  // for (int i = 0; i < 100000; ++i) {
  // double x = gRandom->Gaus(2.5, 1.0);
  // double y = gRandom->Exp(2.);
  // bool pass = (gRandom->Uniform(0., TMath::Gaus(2.5, 2.5, 1.) * TMath::Gaus(2.5, 2.5, 2.)) > TMath::Gaus(x, 2.5, 1.) * TMath::Gaus(y, 2.5, 2.));
  // col->Fill(x, y, pass);
  //}

  col->UpdateEfficiencies();

  auto poly = col->MakePH2Poly();
  auto canvas = new TCanvas("", "", 800, 600);
  poly->Draw("COLZ");
  poly->GetXaxis()->SetTitle("#text{$p$ [MeV]}");
  poly->GetYaxis()->SetTitle("#text{$#eta$}");
  poly->GetZaxis()->SetRangeUser(0., 1.);
  canvas->Print("before.tex");
  canvas->Clear();

  col->BuiltNeighbourhood();
  col->PrintBins();

  col->MergeBins(5.);

  col->UpdateEfficiencies();

  col->PrintBins();

  auto outfile = TFile::Open("output.root","RECREATE");

  auto poly_after = col->MakePH2Poly();
  poly_after->Draw("COLZ");
  poly_after->GetXaxis()->SetTitle("#text{$p$ [MeV]}");
  poly_after->GetYaxis()->SetTitle("#text{$#eta$}");
  poly_after->GetZaxis()->SetRangeUser(0., 1.);
  canvas->Print("after.tex");
  canvas->Clear();

  outfile->WriteTObject(poly_after, "poly");
  outfile->Close();

  auto poly_map = col->MakePH2Poly(true);
  poly_map->Draw("COLZ");
  poly_map->Draw("TEXTsame");
  poly_map->GetXaxis()->SetTitle("#text{$p$ [MeV]}");
  poly_map->GetYaxis()->SetTitle("#text{$#eta$}");
  canvas->Print("map.tex");
  canvas->Clear();

  return 0;
}
