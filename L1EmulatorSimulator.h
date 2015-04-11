#include <iostream>
#include <vector>
#include <algorithm>

#include "TMath.h"

#ifndef L1EMULATORSIMULATOR
#define L1EMULATORSIMULATOR

namespace L1EmulatorSimulator {
  struct cand{
    int pt;
    int eta;
    int phi;
  };

  enum algoVariation {
    nominal,
    zeroWalls,
    doubleSubtraction,
    sigmaSubtraction,
    sigmaSubtractionzeroWalls,
    barrelOnly,
    oneByOne,
    twoByTwo,
    oneByOneANDzeroWalls,
    oneByOneANDzeroWallsANDsigmaSubtraction,
    twoByTwoANDzeroWalls,
    twoByTwoANDzeroWallsANDsigmaSubtraction,
    twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf,
    slidingSubtractionTwoRegionsNoGapANDzeroWalls,
    slidingSubtractionTwoRegionsGapANDzeroWalls,
    slidingSubtractionDoubleTwoRegionsNoGapANDzeroWalls
  };

  int deltaGctPhi(int phi1, int phi2);
  void CaloRingBackgroundSubtraction(cand region[396], cand subregion[396]);
  void CaloSlidingBackgroundSubtractionTwoRegionsNoGap(cand region[396], cand subregion[396]);
  void CaloSlidingBackgroundSubtractionDoubleTwoRegionsNoGap(cand region[396], cand subregion[396]);
  void CaloSlidingBackgroundSubtractionTwoRegionsGap(cand region[396], cand subregion[396]);
  void CaloRingSigmaBackgroundSubtraction(cand region[396], cand subregion[396]);
  void CaloRingSigmaBackgroundSubtraction1half(cand region[396], cand subregion[396]);
  void SlidingWindowJetFinder(cand region[396], cand output[8], algoVariation algo);
  void OneByOneFinder(cand region[396], cand output[8], algoVariation algo);
  void TwoByTwoFinder(cand region[396], cand output[8], algoVariation algo);

  int deltaGctPhi(int phi1, int phi2)
  {
    int diff = phi1 - phi2;
    if (std::abs(phi1 - phi2) == 17) { //18 regions in phi
      diff = -diff/std::abs(diff);
    }
    return diff;
  }
  
  void CaloRingBackgroundSubtraction(cand region[396], cand subregion[396])
  {
    int etaCount[22];
    int puLevelHI[22];
    float r_puLevelHI[22];

    // 22 values of eta
    for(unsigned i = 0; i < 22; ++i)
    {
      puLevelHI[i] = 0;
      r_puLevelHI[i] = 0.0;
      etaCount[i] = 0;
    }

    for(int i = 0; i < 396; ++i){
      r_puLevelHI[region[i].eta] += region[i].pt;
      etaCount[region[i].eta]++;
    }

    for(unsigned i = 0; i < 22; ++i)
    {
      if(etaCount[i] != 18)
	std::cout << "ERROR: wrong number of regions in phi ring." << std::endl;
      puLevelHI[i] = floor(r_puLevelHI[i]/18. + 0.5); // this floating point operation should probably be replaced
    }

    for(int i = 0; i < 396; ++i){
      subregion[i].pt = std::max(0, region[i].pt - puLevelHI[region[i].eta]);
      subregion[i].eta = region[i].eta;
      subregion[i].phi = region[i].phi;
    }
  }


  void CaloSlidingBackgroundSubtractionTwoRegionsNoGap(cand region[396], cand subregion[396])
  {

    int puLevelHI[396];  
    float r_puLevelHI[396]; 
    int array[22][18];
  
    for(unsigned i = 0; i < 396; ++i){
      puLevelHI[i]=0;
      r_puLevelHI[i]=0.;
      array[region[i].eta][region[i].phi] = i;
    }
  
    for(unsigned i = 0; i < 396; ++i){
    
	  int regionEta = region[i].eta; 
	  int regionPhi = region[i].phi;
	  
	  int indexPhill=((regionPhi-2)%18+18)%18;
	  int indexPhil=((regionPhi-1)%18+18)%18;
	  int indexPhir=((regionPhi+1)%18+18)%18;
	  int indexPhirr=((regionPhi+2)%18+18)%18; 

	  int globalindll=array[regionEta][indexPhill];
	  int globalindl=array[regionEta][indexPhil];
	  int globalindr=array[regionEta][indexPhir];
	  int globalindrr=array[regionEta][indexPhirr];
	  
	  r_puLevelHI[i]=(region[globalindll].pt+region[globalindl].pt+region[globalindr].pt+region[globalindrr].pt);
	  puLevelHI[i] = floor(r_puLevelHI[i]/4. + 0.5); 
      subregion[i].pt = std::max(0, region[i].pt - puLevelHI[i]);
      subregion[i].eta = region[i].eta;
      subregion[i].phi = region[i].phi;
    }
  }

  void CaloSlidingBackgroundSubtractionDoubleTwoRegionsNoGap(cand region[396], cand subregion[396])
  {

    int puLevelHI[396];  
    float r_puLevelHI[396]; 
    int array[22][18];
  
    for(unsigned i = 0; i < 396; ++i){
      puLevelHI[i]=0;
      r_puLevelHI[i]=0.;
      array[region[i].eta][region[i].phi] = i;
    }
  
    for(unsigned i = 0; i < 396; ++i){
    
	  int regionEta = region[i].eta; 
	  int regionPhi = region[i].phi;
	  
	  if(region[i].phi%2==0){
	
	    int indexPartner=((regionPhi+1)%18+18)%18;
	    int indexPhill=((regionPhi-2)%18+18)%18;
	    int indexPhil=((regionPhi-1)%18+18)%18;
	    int indexPhir=((regionPhi+2)%18+18)%18;
	    int indexPhirr=((regionPhi+3)%18+18)%18; 

	    int globalindPartner=array[regionEta][indexPartner];
	    int globalindll=array[regionEta][indexPhill];
	    int globalindl=array[regionEta][indexPhil];
	    int globalindr=array[regionEta][indexPhir];
	    int globalindrr=array[regionEta][indexPhirr];
	    r_puLevelHI[i]=(region[globalindll].pt+region[globalindl].pt+region[globalindr].pt+region[globalindrr].pt);
	    r_puLevelHI[globalindPartner]=(region[globalindll].pt+region[globalindl].pt+region[globalindr].pt+region[globalindrr].pt);
	  }
	}
	  
    for(unsigned i = 0; i < 396; ++i){
	  puLevelHI[i] = floor(r_puLevelHI[i]/4. + 0.5); 
      subregion[i].pt = std::max(0, region[i].pt - puLevelHI[i]);
      subregion[i].eta = region[i].eta;
      subregion[i].phi = region[i].phi;
    }
  }

  
  void CaloSlidingBackgroundSubtractionTwoRegionsGap(cand region[396], cand subregion[396])
  {

    int puLevelHI[396];  
    float r_puLevelHI[396]; 
    int array[22][18];
  
    for(unsigned i = 0; i < 396; ++i){
      puLevelHI[i]=0;
      r_puLevelHI[i]=0.;
      array[region[i].eta][region[i].phi] = i;
    }
  
    for(unsigned i = 0; i < 396; ++i){
    
	  int regionEta = region[i].eta; 
	  int regionPhi = region[i].phi;
	  
	  int indexPhill=((regionPhi-3)%18+18)%18;
	  int indexPhil=((regionPhi-2)%18+18)%18;
	  int indexPhir=((regionPhi+2)%18+18)%18;
	  int indexPhirr=((regionPhi+3)%18+18)%18; 

	  int globalindll=array[regionEta][indexPhill];
	  int globalindl=array[regionEta][indexPhil];
	  int globalindr=array[regionEta][indexPhir];
	  int globalindrr=array[regionEta][indexPhirr];
	  
	  r_puLevelHI[i]=(region[globalindll].pt+region[globalindl].pt+region[globalindr].pt+region[globalindrr].pt);
	  puLevelHI[i] = floor(r_puLevelHI[i]/4. + 0.5); 
      subregion[i].pt = std::max(0, region[i].pt - puLevelHI[i]);
      subregion[i].eta = region[i].eta;
      subregion[i].phi = region[i].phi;
    }
  }


  

  void CaloRingSigmaBackgroundSubtraction(cand region[396], cand subregion[396])
  {
    int etaCount[22];
    int puLevelHI[22];
    float r_puLevelHI[22];
    float r_puLevelHI2[22];

    // 22 values of eta
    for(unsigned i = 0; i < 22; ++i)
    {
      puLevelHI[i] = 0;
      r_puLevelHI[i] = 0.0;
      r_puLevelHI2[i] = 0.0;
      etaCount[i] = 0;
    }

    for(int i = 0; i < 396; ++i){
      r_puLevelHI[region[i].eta] += region[i].pt;
      r_puLevelHI2[region[i].eta] += (region[i].pt * region[i].pt);
      etaCount[region[i].eta]++;
    }

    for(unsigned i = 0; i < 22; ++i)
    {
      if(etaCount[i] != 18)
	std::cout << "ERROR: wrong number of regions in phi ring." << std::endl;
      puLevelHI[i] = floor(r_puLevelHI[i]/18. + 0.5); // this floating point operation should probably be replaced
      // also subtract an extra sigma
      puLevelHI[i] += floor(TMath::Sqrt( r_puLevelHI2[i]/18. - (r_puLevelHI[i]/18.)*(r_puLevelHI[i]/18.)) + 0.5);
    }

    for(int i = 0; i < 396; ++i){
      subregion[i].pt = std::max(0, region[i].pt - puLevelHI[region[i].eta]);
      subregion[i].eta = region[i].eta;
      subregion[i].phi = region[i].phi;
    }
  }
  
  void CaloRingSigmaBackgroundSubtraction1half(cand region[396], cand subregion[396])
  {
    int etaCount[22];
    int puLevelHI[22];
    float r_puLevelHI[22];
    float r_puLevelHI2[22];

    // 22 values of eta
    for(unsigned i = 0; i < 22; ++i)
    {
      puLevelHI[i] = 0;
      r_puLevelHI[i] = 0.0;
      r_puLevelHI2[i] = 0.0;
      etaCount[i] = 0;
    }

    for(int i = 0; i < 396; ++i){
      r_puLevelHI[region[i].eta] += region[i].pt;
      r_puLevelHI2[region[i].eta] += (region[i].pt * region[i].pt);
      etaCount[region[i].eta]++;
    }

    for(unsigned i = 0; i < 22; ++i)
    {
      if(etaCount[i] != 18)
	std::cout << "ERROR: wrong number of regions in phi ring." << std::endl;
      puLevelHI[i] = floor(r_puLevelHI[i]/18. + 0.5); // this floating point operation should probably be replaced
      // also subtract an extra sigma
      puLevelHI[i] += floor(1.5 *TMath::Sqrt( r_puLevelHI2[i]/18. - (r_puLevelHI[i]/18.)*(r_puLevelHI[i]/18.)) + 0.5);
    }

    for(int i = 0; i < 396; ++i){
      subregion[i].pt = std::max(0, region[i].pt - puLevelHI[region[i].eta]);
      subregion[i].eta = region[i].eta;
      subregion[i].phi = region[i].phi;
    }
  }


  void SlidingWindowJetFinder(cand region[396], cand output[8], algoVariation algo = nominal)
  {
    if((algo == oneByOne) || (algo == oneByOneANDzeroWalls) || (algo == oneByOneANDzeroWallsANDsigmaSubtraction))
    {
      OneByOneFinder(region, output, algo);
      return;
    }
    else if((algo == twoByTwo) || (algo == twoByTwoANDzeroWalls) || (algo == twoByTwoANDzeroWallsANDsigmaSubtraction) ||(algo == twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf))
    {
      TwoByTwoFinder(region, output, algo);
      return;
    }
    else
    {
      std::vector<cand> forjets;
      std::vector<cand> cenjets;

      for(int i = 0; i < 396; i++) {
	int regionET = region[i].pt;
	int regionEta = region[i].eta;
	int regionPhi = region[i].phi;
	if(algo == zeroWalls || algo == sigmaSubtractionzeroWalls || algo == slidingSubtractionTwoRegionsNoGapANDzeroWalls || algo == slidingSubtractionTwoRegionsGapANDzeroWalls || algo == slidingSubtractionDoubleTwoRegionsNoGapANDzeroWalls)
	{
	  if(regionEta == 4 || regionEta == 17) regionET =0;
	}
	else if (algo == barrelOnly)
	{
	  if(regionEta <=5 || regionEta >= 16) regionET=0;
	}
	int neighborN_et = 0;
	int neighborS_et = 0;
	int neighborE_et = 0;
	int neighborW_et = 0;
	int neighborNE_et = 0;
	int neighborSW_et = 0;
	int neighborNW_et = 0;
	int neighborSE_et = 0;
	unsigned int nNeighbors = 0;
	for(int j = 0; j < 396; j++) {
	  int neighborET = region[j].pt;
	  int neighborEta = region[j].eta;
	  if(algo == zeroWalls || algo == sigmaSubtractionzeroWalls || algo == slidingSubtractionTwoRegionsNoGapANDzeroWalls || algo == slidingSubtractionTwoRegionsGapANDzeroWalls || algo == slidingSubtractionDoubleTwoRegionsNoGapANDzeroWalls)
	  {
	    if(neighborEta == 4 || neighborEta == 17) neighborET =0;
	  }
	  else if(algo == barrelOnly)
	  {
	    if(neighborEta <= 5 || neighborEta >= 16) neighborET =0;
	  }
	  int neighborPhi = region[j].phi;
	  if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
	     (regionEta ) == neighborEta) {
	    neighborN_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
		  (regionEta    ) == neighborEta) {
	    neighborS_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
		  (regionEta + 1) == neighborEta) {
	    neighborE_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
		  (regionEta - 1) == neighborEta) {
	    neighborW_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
		  (regionEta + 1) == neighborEta) {
	    neighborNE_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
		  (regionEta - 1) == neighborEta) {
	    neighborSW_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == 1 &&
		  (regionEta - 1) == neighborEta) {
	    neighborNW_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	  else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
		  (regionEta + 1) == neighborEta) {
	    neighborSE_et = neighborET;
	    nNeighbors++;
	    continue;
	  }
	}
	if(regionET > neighborN_et &&
	   regionET > neighborNW_et &&
	   regionET > neighborW_et &&
	   regionET > neighborSW_et &&
	   regionET >= neighborNE_et &&
	   regionET >= neighborE_et &&
	   regionET >= neighborSE_et &&
	   regionET >= neighborS_et) {
	  unsigned int jetET = regionET +
	    neighborN_et + neighborS_et + neighborE_et + neighborW_et +
	    neighborNE_et + neighborSW_et + neighborSE_et + neighborNW_et;

	  int jetPhi = regionPhi;
	  int jetEta = regionEta;

	  bool neighborCheck = (nNeighbors == 8);
	  // On the eta edge we only expect 5 neighbors
	  if (!neighborCheck && (jetEta == 0 || jetEta == 21) && nNeighbors == 5)
	    neighborCheck = true;

	  if (!neighborCheck) {
	    std::cout << "ERROR (jet finder): Wrong number of neighbor regions." << std::endl;
	    std::cout << "phi: " << jetPhi << " eta: " << jetEta << " n: " << nNeighbors << std::endl;
	  }

	  cand theJet;
	  theJet.pt = jetET / 8; // factor of 8 comes from hardware scale change
	  theJet.eta = jetEta;
	  theJet.phi = jetPhi;

	  const bool forward = (jetEta < 4 || jetEta > 17);

	  if(forward)
	    forjets.push_back(theJet);
	  else
	    cenjets.push_back(theJet);
	}
      }

      auto comp = [&](cand i, cand j)-> bool {
	return (i.pt > j.pt );
      };


      // sort the jet collections and only take the largest 4
      // emulator only outputs top 4 in each category.
      std::sort(forjets.begin(), forjets.end(), comp);
      std::sort(cenjets.begin(), cenjets.end(), comp);
      forjets.resize(4);
      cenjets.resize(4);

      output[0]=cenjets.at(0);
      output[1]=cenjets.at(1);
      output[2]=cenjets.at(2);
      output[3]=cenjets.at(3);
      output[4]=forjets.at(0);
      output[5]=forjets.at(1);
      output[6]=forjets.at(2);
      output[7]=forjets.at(3);
    }
  }

  void OneByOneFinder(cand region[396], cand output[8], algoVariation algo)
  {
    std::vector<cand> forjets;
    std::vector<cand> cenjets;

    for(int i = 0; i < 396; i++) {
      region[i].pt = region[i].pt /8;
      if((algo == oneByOneANDzeroWalls) || (algo == oneByOneANDzeroWallsANDsigmaSubtraction))
      {
	if(region[i].eta == 4 || region[i].eta == 17)
	  region[i].pt = 0;
      }
      if(region[i].eta < 4 || region[i].eta > 17) {
	forjets.push_back(region[i]);
      } else {
	cenjets.push_back(region[i]);
      }
    }

    auto comp = [&](cand i, cand j)-> bool {
      return (i.pt > j.pt );
    };

    // sort the jet collections and only take the largest 4
    // emulator only outputs top 4 in each category.
    std::sort(forjets.begin(), forjets.end(), comp);
    std::sort(cenjets.begin(), cenjets.end(), comp);
    forjets.resize(4);
    cenjets.resize(4);

    output[0]=cenjets.at(0);
    output[1]=cenjets.at(1);
    output[2]=cenjets.at(2);
    output[3]=cenjets.at(3);
    output[4]=forjets.at(0);
    output[5]=forjets.at(1);
    output[6]=forjets.at(2);
    output[7]=forjets.at(3);

  }

  void TwoByTwoFinder(cand region[396], cand output[8], algoVariation algo)
  {
    std::vector<cand> forjets;
    std::vector<cand> cenjets;

    for(int i = 0; i < 396; i++) {
      int regionET = region[i].pt;
      if((algo == twoByTwoANDzeroWalls) || (algo == twoByTwoANDzeroWallsANDsigmaSubtraction) || (algo == twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf))
      {
	if(region[i].eta == 4 || region[i].eta == 17)
	  regionET = 0;
      }
      int regionEta = region[i].eta;
      int regionPhi = region[i].phi;
      int neighborS_et = 0;
      int neighborE_et = 0;
      int neighborSE_et = 0;
      unsigned int nNeighbors = 0;
      for(int j = 0; j < 396; j++) {
	int neighborET = region[j].pt;
	if((algo == twoByTwoANDzeroWalls) || (algo == twoByTwoANDzeroWallsANDsigmaSubtraction) || (algo == twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf))
	{
	  if(region[j].eta == 4 || region[j].eta == 17)
	    neighborET = 0;
	}
	int neighborEta = region[j].eta;
	int neighborPhi = region[j].phi;
	if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
	   (regionEta    ) == neighborEta) {
	  neighborS_et = neighborET;
	  nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(regionPhi, neighborPhi) == 0 &&
		(regionEta + 1) == neighborEta) {
	  neighborE_et = neighborET;
	  nNeighbors++;
	  continue;
	}
	else if(deltaGctPhi(regionPhi, neighborPhi) == -1 &&
		(regionEta + 1) == neighborEta) {
	  neighborSE_et = neighborET;
	  nNeighbors++;
	  continue;
	}
      }

      unsigned int jetET = regionET +
	neighborS_et + neighborE_et + neighborSE_et;

      int jetPhi = regionPhi;
      int jetEta = regionEta;

      cand theJet;
      theJet.pt = jetET / 8; // factor of 8 comes from hardware scale change
      theJet.eta = jetEta;
      theJet.phi = jetPhi;

      const bool forward = (jetEta < 4 || jetEta > 17);

      if(forward)
	forjets.push_back(theJet);
      else
	cenjets.push_back(theJet);
    }


    auto comp = [&](cand i, cand j)-> bool {
      return (i.pt > j.pt );
    };


    // sort the jet collections and only take the largest 4
    // emulator only outputs top 4 in each category.
    std::sort(forjets.begin(), forjets.end(), comp);
    std::sort(cenjets.begin(), cenjets.end(), comp);
    forjets.resize(4);
    cenjets.resize(4);

    output[0]=cenjets.at(0);
    output[1]=cenjets.at(1);
    output[2]=cenjets.at(2);
    output[3]=cenjets.at(3);
    output[4]=forjets.at(0);
    output[5]=forjets.at(1);
    output[6]=forjets.at(2);
    output[7]=forjets.at(3);
  }
};

#endif
