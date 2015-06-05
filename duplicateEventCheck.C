#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>

#include "EventMatchingCMS.h"

void makeTurnOn(TString inHiForestFileName)
{
  TFile *inFile = TFile::Open(inHiForestFileName);
  TTree *fEvtTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");

  Int_t f_evt, f_run, f_lumi;
  Float_t vz;
  Int_t hiBin;
  fEvtTree->SetBranchAddress("evt",&f_evt);
  fEvtTree->SetBranchAddress("run",&f_run);
  fEvtTree->SetBranchAddress("lumi",&f_lumi);
  fEvtTree->SetBranchAddress("vz",&vz);
  fEvtTree->SetBranchAddress("hiBin",&hiBin);

  // Make the event-matching map ************
  EventMatchingCMS *matcher = new EventMatchingCMS();
  int duplicates = 0;
  std::cout << "Begin making map." << std::endl;
  Long64_t f_entries = fEvtTree->GetEntries();
  for(Long64_t j = 0; j < f_entries; ++j)
  {
    fEvtTree->GetEntry(j);
    bool status = matcher->addEvent(f_evt, f_lumi, f_run, j);
    if(status == false)
      duplicates++;
  }
  std::cout << "Finished making map." << std::endl;
  std::cout << "Duplicate entries: " << duplicates << std::endl;
  // **********************

  std::cout << "Total entries: " << f_entries << std::endl;
}

int main(int argc, char **argv)
{
  if(argc == 2)
  {
    makeTurnOn(argv[1]);
    return 0;
  }
  else
  {
    return 1;
  }
}
