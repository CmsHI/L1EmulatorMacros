#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TString.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>

#include <vector>
#include <iostream>

//#include "EventMatchingCMS.h"

const int MAXL1JETS    = 8;
const int MAXJETS      = 500;
const Int_t THRESHOLDS = 30;
const Double_t L1_THRESHOLD[THRESHOLDS] = {0, 4, 8, 12, 16, 20, 24,
					   28, 32, 36, 40, 44, 48,
					   52, 56, 60, 64, 68, 72,
					   76, 80, 84, 88, 92, 96,
					   100, 104, 108, 112, 116};


void makeTurnOn_fromSameFile_Scatter(const char* inL1Name = "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v5.root", 
		const char* inHiForestFileName = "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_HiForestAndEmulator_v5.root", 
		const char* outFileName        = "jetL1RecoPlots", 
		bool montecarlo = false, //only to decide if to plot vs gen or vs reco pt
		int etaRegionCut = 0) // 0: |eta|<2; 1: |eta|<1.4; 2: 1.4<|eta|<2
{
  // gROOT->Macro("/Users/eusmartass/Documents/cmswrite/hin-10-006/logon.C+");

  float etaLow  = 0.;
  float etaHigh = 2.;
  switch (etaRegionCut) 
    {
    case 1:
      etaHigh = 1.4;
      break;
    case 2: 
      etaLow = 1.4;
      break;
    };
 
  TFile *inFile    = TFile::Open(Form("%s",inHiForestFileName));
  TTree *f1Tree    = (TTree*)inFile->Get("akPu3CaloJetAnalyzer/t"); 
  TTree *fEvtTree  = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree"); 
  TTree *fSkimTree = (TTree*)inFile->Get("skimanalysis/HltTree");
  
  TFile *inL1File  = TFile::Open(Form("%s",inL1Name));
  TTree *l1Tree    = (TTree*)inL1File->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  //f1Tree->AddFriend(l1Tree);

  Int_t l1_event, l1_run, l1_lumi;
  Int_t l1_hwPt[MAXL1JETS], l1_hwEta[MAXL1JETS], l1_hwPhi[MAXL1JETS];
  Double_t l1_pt[MAXL1JETS];

  l1Tree->SetBranchStatus("*"           , 0);
  l1Tree->SetBranchStatus("event",     1);
  l1Tree->SetBranchStatus("lumi",      1);
  l1Tree->SetBranchStatus("jet_hwPt",  1);
  l1Tree->SetBranchStatus("jet_hwEta", 1);
  l1Tree->SetBranchStatus("jet_hwPhi",  1);
  l1Tree->SetBranchStatus("jet_pt",     1);

  l1Tree->SetBranchAddress("event",    &l1_event);
  l1Tree->SetBranchAddress("run",      &l1_run);
  l1Tree->SetBranchAddress("lumi",     &l1_lumi);
  l1Tree->SetBranchAddress("jet_hwPt", l1_hwPt);
  l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  l1Tree->SetBranchAddress("jet_pt",   l1_pt);
  
  // f1Tree->SetBranchAddress("event",&l1_event);
  // f1Tree->SetBranchAddress("run",&l1_run);
  // f1Tree->SetBranchAddress("lumi",&l1_lumi);
  // f1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  // f1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  // f1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  // f1Tree->SetBranchAddress("jet_pt",l1_pt);

  Int_t f_evt, f_run, f_lumi;
  Float_t vz;
  Int_t hiBin;
  fEvtTree->SetBranchStatus("*",  0);
  fEvtTree->SetBranchStatus("evt",1);
  fEvtTree->SetBranchStatus("run",1);
  fEvtTree->SetBranchStatus("lumi",1);
  fEvtTree->SetBranchStatus("vz",  1);
  fEvtTree->SetBranchStatus("hiBin",1);

  fEvtTree->SetBranchAddress("evt",  &f_evt);
  fEvtTree->SetBranchAddress("run" , &f_run);
  fEvtTree->SetBranchAddress("lumi", &f_lumi);
  fEvtTree->SetBranchAddress("vz",   &vz);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t f_num;
  Float_t f_pt[MAXJETS];
  Float_t f_eta[MAXJETS];
  Float_t f_phi[MAXJETS];
  Float_t f_rawpt[MAXJETS];

  Int_t num_gen;
  Float_t genpt[MAXJETS], geneta[MAXJETS], genphi[MAXJETS];
  f1Tree->SetBranchStatus("*", 0);
  f1Tree->SetBranchStatus("nref",1);
  f1Tree->SetBranchStatus("jtpt", 1);
  f1Tree->SetBranchStatus("jteta",1);
  f1Tree->SetBranchStatus("jtphi",1);
  f1Tree->SetBranchStatus("rawpt",1);

  f1Tree->SetBranchAddress("nref", &f_num);
  f1Tree->SetBranchAddress("jtpt", f_pt);
  f1Tree->SetBranchAddress("jteta",f_eta);
  f1Tree->SetBranchAddress("jtphi",f_phi);
  f1Tree->SetBranchAddress("rawpt",f_rawpt);

  if(montecarlo)
  {
    f1Tree->SetBranchStatus("*",  0);
    f1Tree->SetBranchStatus("ngen", 1);
    f1Tree->SetBranchStatus("genpt",1);
    f1Tree->SetBranchStatus("geneta",1);
    f1Tree->SetBranchStatus("genphi",1);

    f1Tree->SetBranchAddress("ngen",  &num_gen);
    f1Tree->SetBranchAddress("genpt", genpt);
    f1Tree->SetBranchAddress("geneta",geneta);
    f1Tree->SetBranchAddress("genphi",genphi);
  }

  const int nBins    = 75;
  const double maxPt = 300;

  TH1D *l1Pt = new TH1D("l1Pt",";L1 p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fPt[3];
  fPt[0] = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_1");
  fPt[2] = (TH1D*)fPt[0]->Clone("fPt_2");

  TH1D *fPtRaw[3];
  fPtRaw[0] = new TH1D("fPtRaw_0",";offline raw p_{T} (GeV)",nBins,0,maxPt);
  fPtRaw[1] = (TH1D*)fPtRaw[0]->Clone("fPtRaw_1");
  fPtRaw[2] = (TH1D*)fPtRaw[0]->Clone("fPtRaw_2");

  TH2D *corr[3];
  corr[0]  = new TH2D("corr_0",";offline p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);
  corr[1]  = (TH2D*)corr[0]->Clone("corr_1");
  corr[2]  = (TH2D*)corr[0]->Clone("corr_2");

  TH2D *corr_raw[3];
  corr_raw[0]  = new TH2D("corr_raw_0",";offline raw p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);
  corr_raw[1]  = (TH2D*)corr_raw[0]->Clone("corr_raw_1");
  corr_raw[2]  = (TH2D*)corr_raw[0]->Clone("corr_raw_2");

  // turn-on efficiencies curves
  TH1D *accepted[THRESHOLDS][3];
  for(int i = 0; i < THRESHOLDS; ++i)
  {
    for(int t = 0; t < 3; ++t)
    {
      accepted[i][t] = new TH1D(Form("accepted_pt%d_%d",(int)L1_THRESHOLD[i],t),";offline p_{T}",nBins,0,maxPt);
    }
  }
  
  Long64_t entries = f1Tree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    // Only use good collision events ********
    fEvtTree->GetEntry(j);
    fSkimTree->GetEntry(j);
    bool goodEvent = false;
    // 5.02 TeV Hydjet missing pcollisionEventSelection for now
    if( (pcollisionEventSelection == 1) && (TMath::Abs(vz) < 15) )
      {
      goodEvent = true;
    }
    if(!goodEvent) continue;

    //std::cout << accepted[0][0]->GetName() << std::endl;
    l1Tree->GetEntry(j); // l1 tree
    f1Tree->GetEntry(j); // offline jet tree

    //std::cout << accepted[0][0]->GetName() << std::endl;

    // search for leading jets in each trees
    // L1 tree
    double maxl1pt = 0.;
    for(int i = 0; i < MAXL1JETS; ++i)
    {
      // if(TMath::Abs(l1_hwEta[i]) > etaHigh || TMath::Abs(l1_hwEta[i]) < etaLow) continue;
      if(l1_pt[i] > maxl1pt)
	maxl1pt = l1_pt[i];
    }
    // cout<<"Event "<<j<<"leading jet en "<<maxl1pt<<endl;
    l1Pt->Fill(maxl1pt);

    // jet tree (reco or gen)
    double maxfpt     = 0;
    double maxfpt_raw = 0;
    if(!montecarlo)
    {
      for(int i = 0; i < f_num; ++i)
      {
       	if(TMath::Abs(f_eta[i]) > etaHigh || TMath::Abs(f_eta[i]) < etaLow) continue;
    	if(f_pt[i] > maxfpt)        { maxfpt = f_pt[i];} // JEC applied on reco
    	if(f_rawpt[i] > maxfpt_raw) { maxfpt_raw = f_rawpt[i];} // raw jet energy
      }
    }
    else
    {
      for(int i = 0; i < num_gen; ++i)
      {
    	if(TMath::Abs(geneta[i]) > etaHigh || TMath::Abs(geneta[i]) < etaLow) continue;
    	if(genpt[i] > maxfpt)
    	  maxfpt = genpt[i];
      }
    }

    // fill the 0-100% histogram with the leading jet
    fPt[0]->Fill(maxfpt);
    fPtRaw[0]->Fill(maxfpt_raw);
    // fill 2D pt_l1 vs pt_reco histograms
    // 0-100%
    corr[0]->Fill(maxfpt,maxl1pt);
    corr_raw[0]->Fill(maxfpt_raw,maxl1pt);

    // same for 2 different centrality regions
    if(hiBin < 60) // centrality 0-30%
      {
    	fPt[1]->Fill(maxfpt); 
    	fPtRaw[1]->Fill(maxfpt_raw);    
    	corr[1]->Fill(maxfpt,maxl1pt);
    	corr_raw[1]->Fill(maxfpt_raw,maxl1pt);
      }
    else if (hiBin >= 100) // centrality 50-100%
      {
    	fPt[2]->Fill(maxfpt); 
    	fPtRaw[2]->Fill(maxfpt_raw);
    	corr[2]->Fill(maxfpt,maxl1pt);
    	corr_raw[2]->Fill(maxfpt_raw,maxl1pt);
      }

    // fill numerator 
    for(int i = 0; i < THRESHOLDS; ++i)
    {
      if(maxl1pt>L1_THRESHOLD[i])
      {
    	accepted[i][0]->Fill(maxfpt);// cent 0-100%
    	if(hiBin < 60) // cent 0-30%
    	  accepted[i][1]->Fill(maxfpt);
    	else if (hiBin >= 100) // cent 50-100%
    	  accepted[i][2]->Fill(maxfpt);
      }
    }

  }// tree entries
  
  //############### calcualte the efficiencies turn-on curves
  // denominator--all reco jets with pt>threshold; 
  // numerator: all L1 jets with pt>threshold   
  TGraphAsymmErrors *a[THRESHOLDS][3];
  for(int k = 0; k < THRESHOLDS; ++k)
    {
      for(int l = 0; l < 3; ++l)
	{
	  a[k][l] = new TGraphAsymmErrors();
	  a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
	  a[k][l]->SetName(Form("asymm_pt_%d_%d",(int)L1_THRESHOLD[k],l));
	}
    }

  // some drawing: 
  TLatex *tex2 = new TLatex(0.25,0.85,Form("%.1f<|#eta|<%.1f",etaLow,etaHigh));
  tex2->SetNDC();
  tex2->SetTextAlign(13);
  tex2->SetTextFont(43);
  tex2->SetTextSize(17);
  tex2->SetLineWidth(2);

  TLatex *lt1 = new TLatex(); lt1->SetNDC();
 
  // ---------------------------------------------------------------------------
  TCanvas *pc = new TCanvas("pc","pc",1200,600);
  pc->Divide(3,1);
  pc->cd(1);
  corr[0]->Draw("colz");
  tex2->Draw();
  lt1->DrawLatex(0.25,0.75,Form("Cent.0-100 %%"));

  pc->cd(2);
  corr[1]->Draw("colz");
  lt1->DrawLatex(0.25,0.75,Form("Cent.0-30 %%"));

  pc->cd(3);
  corr[2]->Draw("colz");
  lt1->DrawLatex(0.25,0.75,Form("Cent.50-100 %%"));

  pc->SaveAs(Form("pcPtL1Reco_etaRegion%d.png",etaRegionCut));
  // ---------------------------------------------------------------------------

  // ############### write out all the histograms
  TFile *outFile   = new TFile(Form("%s_etaRegion%d_mc%d.root",outFileName,etaRegionCut,montecarlo),"RECREATE");
  outFile->cd();
  l1Pt->Write();
  for (int i=0; i<3; i++)
    {
      fPt[i]->Write();
      corr[i]->Write();
      fPtRaw[i]->Write();
      corr_raw[i]->Write();
    }

  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      accepted[k][l]->Write();
      a[k][l]->Write();
    }
  }

  outFile->Close();
  
}

/// / ------------- main piece
// int main(int argc, char **argv)
// {
//   if(argc == 4)
//   {
//     makeTurnOn(argv[1], argv[2], argv[3]);
//     return 0;
//   }
//   else if(argc == 5)
//   {
//     makeTurnOn(argv[1], argv[2], argv[3], atoi(argv[4]));
//   }
//   else
//   {
//     std::cout << "Usage: \nmakeTurnOn_fromSameFile.exe <input_HiForest_file> <output_file>" << std::endl;
//     return 1;
//   }
