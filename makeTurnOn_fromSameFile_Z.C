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
//#include <string>

//#include "L1EmulatorSimulator.h"

//using namespace L1EmulatorSimulator;

const int MAXL1EMCANDS = 144;
const int MAXL1REGIONS = 396;
const int MAXL1JETS = 8;
const int MAXPHOTONS = 500;
const Int_t THRESHOLDS = 6;
const Int_t MAXGEN = 200000;

const TString THRESHNAMES[THRESHOLDS] = {"10_10",
					 "15_10",
					 "15_15",
					 "20_10",
					 "20_15",
					 "20_20"};

void makeTurnOn(/*TString inHiForestFileName, TString outFileName*/)
{
  const TString inHiForestFileName = "/mnt/hadoop/cms/store/user/richard/L1Emulator/twang-ZEE_5TeV_HiForest_v2/HiForest_54_1_634.root";
  const TString outFileName = "hist_Zee_Z_emcand.root";
  TFile *outFile = new TFile(outFileName,"RECREATE");

  TFile *inFile = TFile::Open(inHiForestFileName);
  TTree *f1Tree = (TTree*)inFile->Get("HiGenParticleAna/hi");
  TTree *fEvtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");
  TTree *fSkimTree = (TTree*)inFile->Get("skimanalysis/HltTree");
  TTree *l1Tree = (TTree*)inFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  Int_t l1_event, l1_run, l1_lumi;
  Int_t l1_hwPt[MAXL1JETS], l1_hwEta[MAXL1JETS], l1_hwPhi[MAXL1JETS];
  Double_t l1_pt[MAXL1JETS];

  l1Tree->SetBranchAddress("event",&l1_event);
  l1Tree->SetBranchAddress("run",&l1_run);
  l1Tree->SetBranchAddress("lumi",&l1_lumi);
  l1Tree->SetBranchAddress("jet_hwPt",l1_hwPt);
  l1Tree->SetBranchAddress("jet_hwEta",l1_hwEta);
  l1Tree->SetBranchAddress("jet_hwPhi",l1_hwPhi);
  l1Tree->SetBranchAddress("jet_pt",l1_pt);

  Int_t emcand_hwPt[MAXL1EMCANDS], emcand_hwEta[MAXL1EMCANDS], emcand_hwPhi[MAXL1EMCANDS];
  Int_t region_hwPt[MAXL1REGIONS], region_hwEta[MAXL1REGIONS], region_hwPhi[MAXL1REGIONS];

  l1Tree->SetBranchAddress("emcand_hwPt",emcand_hwPt);
  l1Tree->SetBranchAddress("emcand_hwEta",emcand_hwEta);
  l1Tree->SetBranchAddress("emcand_hwPhi",emcand_hwPhi);
  l1Tree->SetBranchAddress("region_hwPt",region_hwPt);
  l1Tree->SetBranchAddress("region_hwEta",region_hwEta);
  l1Tree->SetBranchAddress("region_hwPhi",region_hwPhi);

  Int_t f_evt, f_run, f_lumi;
  Float_t vz;
  Int_t hiBin;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  fEvtTree->SetBranchAddress("vz",&vz);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fSkimTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fSkimTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  Int_t mult;
  Float_t pt[MAXGEN], eta[MAXGEN], phi[MAXGEN];
  Int_t pdg[MAXGEN], sta[MAXGEN];
  Int_t nDaughters[MAXGEN];
  Int_t daughterIDx[MAXGEN][200];

  f1Tree->SetBranchAddress("mult",&mult);
  f1Tree->SetBranchAddress("pt",pt);
  f1Tree->SetBranchAddress("eta",eta);
  f1Tree->SetBranchAddress("phi",phi);
  f1Tree->SetBranchAddress("pdg",pdg);
  f1Tree->SetBranchAddress("sta",sta);
  f1Tree->SetBranchAddress("nDaughters",nDaughters);
  f1Tree->SetBranchAddress("daughterIDx",daughterIDx);

  const int nBins = 100;
  const int maxPt = 100;

  TH1D *l1Pt = new TH1D("l1Pt",";L1 p_{T} (GeV)",nBins,0,maxPt);
  TH1D *fPt[3];
  fPt[0] = new TH1D("fPt_0",";offline p_{T} (GeV)",nBins,0,maxPt);
  fPt[1] = (TH1D*)fPt[0]->Clone("fPt_1");
  fPt[2] = (TH1D*)fPt[0]->Clone("fPt_2");
  TH1D *accepted[THRESHOLDS][3];

  for(int i = 0; i < THRESHOLDS; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      accepted[i][j] = new TH1D(Form("accepted_pt%s_%d",THRESHNAMES[i].Data(),j),";offline p_{T}",nBins,0,maxPt);
    }
  }

  TH2D *corr = new TH2D("corr",";offline p_{T};l1 p_{T}",nBins,0,maxPt,nBins,0,maxPt);
  TH2D *deltaMap = new TH2D("deltamap",";#Delta #eta;#Delta #phi",100, 0, 10, 100, 0, 3.15);

  Long64_t entries = f1Tree->GetEntries();
  for(Long64_t j = 0; j < entries; ++j)
  {
    if(j % 10000 == 0)
      printf("%lld / %lld\n",j,entries);

    fEvtTree->GetEntry(j);
    fSkimTree->GetEntry(j);

    l1Tree->GetEntry(j);
    f1Tree->GetEntry(j);

    double maxl1pt = 0;
    double maxl1pt2 = 0;
    double maxl1eta = -10;
    double maxl1phi = -10;
    // switch(seed){
    // case emcands:
    // {
      for(int i = 0; i < MAXL1EMCANDS; ++i)
      {
	//if(emcand_hwEta[i] < 7 || emcand_hwEta[i] > 14) continue;
	if(emcand_hwPt[i] > maxl1pt)
	{
	  maxl1pt2 = maxl1pt;
	  maxl1pt = emcand_hwPt[i];
	  //maxl1eta = L1EmulatorSimulator::physicalEta(emcand_hwEta[i]);
	  //maxl1phi = L1EmulatorSimulator::physicalPhi(emcand_hwPhi[i]);
	} else if (emcand_hwPt[i] > maxl1pt2)
	{
	  maxl1pt2 = emcand_hwPt[i];
	}
      }
    // }
    // break;
    // case regions:
    // {
    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	if(region_hwPt[i] * 0.5 > maxl1pt)
    // 	{
    // 	  if(region_hwEta[i] < 6 || region_hwEta[i] > 15) continue;
    // 	  maxl1pt = region_hwPt[i] * 0.5;
    // 	  maxl1eta = physicalEta(region_hwEta[i]);
    // 	  maxl1phi = physicalPhi(region_hwPhi[i]);
    // 	}
    //   }
    // }
    // break;
    // case subRegions:
    // {
    //   L1EmulatorSimulator::cand regions[396];
    //   L1EmulatorSimulator::cand subregions[396];
    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	regions[i].pt = region_hwPt[i];
    // 	regions[i].eta = region_hwEta[i];
    // 	regions[i].phi = region_hwPhi[i];
    //   }
    //   L1EmulatorSimulator::CaloRingBackgroundSubtraction(regions, subregions);

    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	if(subregions[i].eta < 6 || subregions[i].eta > 15) continue;
    // 	if(subregions[i].pt*0.5 > maxl1pt)
    // 	{
    // 	  maxl1pt = subregions[i].pt * 0.5;
    // 	  maxl1eta = physicalEta(subregions[i].eta);
    // 	  maxl1phi = physicalPhi(subregions[i].phi);
    // 	}
    //   }
    // }
    // break;
    // case twoByTwoJets:
    // {
    //   L1EmulatorSimulator::cand regions[396];
    //   L1EmulatorSimulator::cand subregions[396];
    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	regions[i].pt = region_hwPt[i];
    // 	regions[i].eta = region_hwEta[i];
    // 	regions[i].phi = region_hwPhi[i];
    //   }
    //   L1EmulatorSimulator::CaloRingBackgroundSubtraction(regions, subregions);
    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	if(subregions[i].eta < 6 || subregions[i].eta > 15)
    // 	  subregions[i].pt = 0;
    //   }
    //   L1EmulatorSimulator::cand jets[8];
    //   L1EmulatorSimulator::TwoByTwoFinder(subregions, jets);

    //   for(int i = 0; i < MAXL1JETS; ++i)
    //   {
    // 	if(jets[i].pt*4 > maxl1pt)
    // 	{
    // 	  maxl1pt = jets[i].pt * 4;
    // 	  maxl1eta = physicalEta(jets[i].eta);
    // 	  maxl1phi = physicalPhi(jets[i].phi);
    // 	}
    //   }
    // }
    // break;
    // case threeByThreeJets:
    // {
    //   L1EmulatorSimulator::cand regions[396];
    //   L1EmulatorSimulator::cand subregions[396];
    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	regions[i].pt = region_hwPt[i];
    // 	regions[i].eta = region_hwEta[i];
    // 	regions[i].phi = region_hwPhi[i];
    //   }
    //   L1EmulatorSimulator::CaloRingBackgroundSubtraction(regions, subregions);
    //   for(int i = 0; i < MAXL1REGIONS; ++i)
    //   {
    // 	if(subregions[i].eta < 6 || subregions[i].eta > 15)
    // 	  subregions[i].pt = 0;
    //   }
    //   L1EmulatorSimulator::cand jets[8];
    //   L1EmulatorSimulator::SlidingWindowJetFinder(subregions, jets);

    //   for(int i = 0; i < MAXL1JETS; ++i)
    //   {
    // 	if(jets[i].pt*4 > maxl1pt)
    // 	{
    // 	  maxl1pt = jets[i].pt * 4;
    // 	  maxl1eta = physicalEta(jets[i].eta);
    // 	  maxl1phi = physicalPhi(jets[i].phi);
    // 	}
    //   }
    // }
    // break;
    // default:
    //   std::cout << "Invalid seed type specified." << std::endl;
    //   return;
    // }

    double maxfpt = 0;
    double maxfeta = -10;
    double maxfphi = -10;
    for(int i = 0; i < mult; ++i)
    {
      if(pdg[i] != 23) continue;
      if(abs(pdg[daughterIDx[i][0]]) != 11) continue;
      if(pt[i] > maxfpt)
      {
	maxfpt = pt[i];
	maxfeta = eta[i];
	maxfphi = phi[i];
      }
    }
    l1Pt->Fill(maxl1pt);

    bool goodEvent = false;
    if((pcollisionEventSelection == 1) && (TMath::Abs(vz) < 15))
    {
      goodEvent = true;
    }
    if(!goodEvent) continue;

    fPt[0]->Fill(maxfpt);
    if(hiBin < 60)
      fPt[1]->Fill(maxfpt);
    else if (hiBin >= 100)
      fPt[2]->Fill(maxfpt);

    corr->Fill(maxfpt,maxl1pt);
    deltaMap->Fill(TMath::Abs(maxl1eta-maxfeta), TMath::ACos(TMath::Cos(maxl1phi - maxfphi)));

    // for(int i = 0; i < THRESHOLDS; ++i)
    // {
    //   if(maxl1pt >= i)
    //   {
    // 	accepted[i][0]->Fill(maxfpt);
    // 	if(hiBin < 60)
    // 	  accepted[i][1]->Fill(maxfpt);
    // 	else if (hiBin >= 100)
    // 	  accepted[i][2]->Fill(maxfpt);
    //   }
    // }
    if(maxl1pt >= 10 && maxl1pt2 >= 10)
    {
      accepted[0][0]->Fill(maxfpt);
      if(hiBin < 60)
	accepted[0][1]->Fill(maxfpt);
      else if (hiBin >= 100)
	accepted[0][2]->Fill(maxfpt);
    }
    if(maxl1pt >= 15 && maxl1pt2 >= 10)
    {
      accepted[1][0]->Fill(maxfpt);
      if(hiBin < 60)
	accepted[1][1]->Fill(maxfpt);
      else if (hiBin >= 100)
	accepted[2][2]->Fill(maxfpt);
    }
    if(maxl1pt >= 15 && maxl1pt2 >= 15)
    {
      accepted[2][0]->Fill(maxfpt);
      if(hiBin < 60)
	accepted[2][1]->Fill(maxfpt);
      else if (hiBin >= 100)
	accepted[2][2]->Fill(maxfpt);
    }
    if(maxl1pt >= 20 && maxl1pt2 >= 10)
    {
      accepted[3][0]->Fill(maxfpt);
      if(hiBin < 60)
	accepted[3][1]->Fill(maxfpt);
      else if (hiBin >= 100)
	accepted[3][2]->Fill(maxfpt);
    }
    if(maxl1pt >= 20 && maxl1pt2 >= 15)
    {
      accepted[4][0]->Fill(maxfpt);
      if(hiBin < 60)
	accepted[4][1]->Fill(maxfpt);
      else if (hiBin >= 100)
	accepted[4][2]->Fill(maxfpt);
    }
    if(maxl1pt >= 20 && maxl1pt2 >= 20)
    {
      accepted[5][0]->Fill(maxfpt);
      if(hiBin < 60)
	accepted[5][1]->Fill(maxfpt);
      else if (hiBin >= 100)
	accepted[5][2]->Fill(maxfpt);
    }

  }

  TGraphAsymmErrors *a[THRESHOLDS][3];
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      a[k][l] = new TGraphAsymmErrors();
      a[k][l]->BayesDivide(accepted[k][l],fPt[l]);
      a[k][l]->SetName(Form("asymm_pt_%s_%d",THRESHNAMES[k].Data(),l));
    }
  }

  outFile->cd();

  l1Pt->Write();
  fPt[0]->Write();
  fPt[1]->Write();
  fPt[2]->Write();
  corr->Write();
  deltaMap->Write();
  for(int k = 0; k < THRESHOLDS; ++k){
    for(int l = 0; l < 3; ++l)
    {
      accepted[k][l]->Write();
      a[k][l]->Write();
    }
  }

  outFile->Close();
}

int main(/*int argc, char **argv*/)
{
  // if(argc == 3)
  // {
  //   makeTurnOn(argv[1], argv[2]);
  //   return 0;
  // }
  // else
  // {
    makeTurnOn();
    return 0;
    //}
}
