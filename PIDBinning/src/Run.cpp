#include "EfficiencyBin.h"
#include "EfficiencyBinCollection.h"
#include "TCanvas.h"
#include "TRandom3.h"

#include "TH2Poly.h"

int main() {
  // Let's create a collection object
  auto col = new EfficiencyBinCollection();

  // Put in a few squares
  for (int iX = 0; iX < 50; iX++) {
    for (int iY = 0; iY < 50; iY++) {
      col->AddBin(EfficiencyBin::SquareBin(0.1*iX, 0.1*iY, 0.1 ));
    }
  }


  //for (int i = 0; i < 100000; ++i) {
    //double x = gRandom->Uniform(0., 5.);
    //double y = gRandom->Uniform(0., 5.);
    //bool pass = (gRandom->Uniform(0., 1.) > (0.1*x+0.1*y));
    //col->Fill(x, y, pass);
  //}

  //col->PrintBins();

  //col->UpdateEfficiencies();
  //col->PrintBins();

  //auto poly = col->MakePH2Poly();
  //auto canvas = new TCanvas("", "", 800, 600);
  //poly->Draw("COLZ");
  //canvas->Print("test.pdf");

  return 0;
}
