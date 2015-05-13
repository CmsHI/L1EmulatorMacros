#!/bin/bash

set -x

  # enum seedObject {
  #   emcands,
  #   regions,
  #   subRegions,
  #   twoByTwoJets,
  #   threeByThreeJets
  # };


InputHiForest=("/mnt/hadoop/cms/store/user/luck/L1Emulator/minbiasForest_merged_v2/HiForest_PbPb_Data_minbias_fromSkim_v3.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/HydjetMB_502TeV_740pre8_MCHI2_74_V3_rctconfigNoCuts_HiForestAndEmulatorAndHLT_v7.root" "/mnt/hadoop/cms/store/user/luck/L1Emulator/AllQCDPhoton30_PhotonFilter20GeV_eta3_TuneZ2_PbPb_5020GeV_embedded_rctconfigNoCuts_HiForest.root")

SampleType=(MBData Hydjet502 Photon502)

EtaCuts=(1.44 2.0)

OfflineIsolation=(0 1)

# compile the macros with g++
g++ makeTurnOn_fromSameFile_photons.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_fromSameFile_photons.exe || exit 1
g++ makeTurnOn_photons.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o makeTurnOn_photons.exe || exit 1
g++ findthes_photons.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o findthes_photons.exe || exit 1


for SAMPLE in 1 2
do
    for ETA in 0
    do
	for ISO in 0 1
	do
	    isolabel=""
	    if [[ $ISO == 1 ]]
	    then
		isolabel="iso"
	    fi
	    output="hist_${SampleType[SAMPLE]}_${isolabel}photons_eta${EtaCuts[ETA]}_emcandsBarrel.root"
	    rateoutput="rate_${SampleType[SAMPLE]}_${isolabel}photons_eta${EtaCuts[ETA]}_emcandsBarrel.root"
	    if [[ $SAMPLE == 0 ]]
	    then
		./makeTurnOn_photons.exe
	    else
		./makeTurnOn_fromSameFile_photons.exe "${InputHiForest[SAMPLE]}" "$output" 0 "${EtaCuts[ETA]}" "$ISO"
	    fi
	    ./findthes_photons.exe "${InputHiForest[1]}" "$output" "$rateoutput" 0 || exit 1
	done
    done
done
