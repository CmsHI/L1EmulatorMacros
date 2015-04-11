#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>

void plotTurnOn(int option=1)
{
  const Int_t INPUTFILES = 8;

  TString outFileTag = "Hydjet502_500Hz";
  TString inFileName[INPUTFILES] = {"hist_Hydjet502_zeroWalls.root",
                    "hist_Hydjet502_sigmaSubtractionzeroWalls.root",
  				    "hist_Hydjet502_twoByTwoANDzeroWallsANDsigmaSubtraction.root",
  				    "hist_Hydjet502_twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf.root",
  				    "hist_Hydjet502_oneByOneANDzeroWallsANDsigmaSubtraction.root",
  				    "hist_Hydjet502_slidingSubtractionTwoRegionsNoGapANDzeroWalls.root",
  				    "hist_Hydjet502_slidingSubtractionTwoRegionsGapANDzeroWalls.root",
  				    "hist_Hydjet502_slidingSubtractionDoubleTwoRegionsNoGapANDzeroWalls.root"
  };
  TString labels[INPUTFILES] = {"zeroWalls, L1_48, 525Hz",
                "zeroWalls sigma subtraction, L1_24, 370 Hz",
  				"2x2 zeroWalls sigma subtraction, L1_20, 632Hz",
                "2x2 zeroWalls 1.5 sigma subtraction, L1_16, 348 Hz",
  				"1x1 zeroWalls sigma subtraction, L1_16, 423Hz",
  				"zeroWalls sliding bkg no gap, L1_40, 653 Hz",
  				"zeroWalls sliding bkg gap, L1_56, 426 Hz",
  				"zeroWalls sliding DOUBLE bkg no gap , L1_48, 539 Hz"};
  const Double_t L1_THRESHOLD[INPUTFILES] = {48,24,20,16,16,40,56,48};

  TFile *inFile[INPUTFILES];
  const Int_t COLORS[INPUTFILES] = {kBlack, kRed, kRed+2, kBlue, kBlue-2, kGreen, kGreen-2, kMagenta};//, kMagenta+3};
  TGraphAsymmErrors *asymm[INPUTFILES];//[2];

  for(int i = 0; i < INPUTFILES; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      inFile[i] = TFile::Open(inFileName[i]);
      asymm[i] = (TGraphAsymmErrors*)inFile[i]->Get(Form("asymm_pt_%d_0",(int)L1_THRESHOLD[i]));
      asymm[i]->SetMarkerColor(COLORS[i]);
      asymm[i]->SetLineColor(COLORS[i]);
      asymm[i]->SetLineWidth(2);
    }
    //asymm[i][1]->SetMarkerStyle(25);
  }

  // these values MUST MATCH those used in makeTurnOn.C
  const int nBins = 75;
  const double maxPt = 300;

  TH1D *hEmpty = new TH1D("hEmpty",Form(";Jet p_{T} (GeV);Efficiency"),nBins,0,maxPt);

  TCanvas *c1 = new TCanvas();
  hEmpty->SetMinimum(0);
  hEmpty->SetMaximum(1.2);
  hEmpty->Draw();
  //c1->SetLogy();

  TLine *line = new TLine(0,1,maxPt,1);
  line->Draw();

  for(int i = 0; i < INPUTFILES; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      if(option==0 && (i==0 || i==2 ||i==5 || i==6 || i==7)) asymm[i]->Draw("p");
      if(option==1 && (i==1 || i==2 || i==3 || i==4)) asymm[i]->Draw("p");
    }
  }

  TLegend *leg = new TLegend(0.55,0.2,0.9,0.5,"|#eta| < 2");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(18);

  for(int i = 0; i < INPUTFILES; i++)
  {
    //for(int j = 0; j < 2; j++)
    {
      if(option==0 && (i==0 || i==2 ||i==5 || i==6 || i==7)) leg->AddEntry(asymm[i], Form("%s", labels[i].Data()), "lp");
      if(option==1 && (i==1 || i==2 || i==3 || i==4)) leg->AddEntry(asymm[i], Form("%s", labels[i].Data()), "lp");
    }
  }

  leg->Draw();

  c1->SaveAs(Form("%s_turnon.pdf",outFileTag.Data()));

}

int main(int argc, char **argv)
{
  if(argc == 1)
  {
    plotTurnOn();
    return 0;
  }
  else
  {
    std::cout << "Usage:\nplotTurnOn_multiMethods.exe <output_tag>" << std::cout;
    std::cout << "An output pdf will be named <output_tag>_turnon.pdf" << std::cout;
    return 1;
  }
}
