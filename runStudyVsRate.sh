#!/bin/bash

set -x

# This is a master file to run the emulator macro and then print out rates on multiple datasets at once
# You can find the up-to-date list of samples here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TValidationSamples

  # enum algoVariation {
  #0   nominal,
  #1   zeroWalls,
  #2   doubleSubtraction,
  #3   sigmaSubtraction,
  #4   sigmaSubtractionzeroWalls,
  #5   barrelOnly,
  #6   oneByOne,
  #7   twoByTwo,
  #8   oneByOneANDzeroWalls,
  #9   oneByOneANDzeroWallsANDsigmaSubtraction,
  #10  twoByTwoANDzeroWalls,
  #11  twoByTwoANDzeroWallsANDsigmaSubtraction
  #12  twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf
  #13  slidingSubtractionTwoRegionsNoGapANDzeroWalls
  #14  slidingSubtractionTwoRegionsGapANDzeroWalls
  #15  slidingSubtractionDoubleTwoRegionsNoGapANDzeroWalls
  # };

InputType=(MBData Hydjet276 Hydjet502 Hydjet502Dijet30)

AlgoVariations=(nominal zeroWalls doubleSubtraction sigmaSubtraction sigmaSubtractionzeroWalls barrelOnly oneByOne twoByTwo oneByOneAndzeroWalls oneByOneANDzeroWallsANDsigmaSubtraction twoByTwoANDzeroWalls twoByTwoANDzeroWallsANDsigmaSubtraction twoByTwoANDzeroWallsANDsigmaSubtraction1sigmahalf slidingSubtractionTwoRegionsNoGapANDzeroWalls slidingSubtractionTwoRegionsGapANDzeroWalls slidingSubtractionDoubleTwoRegionsNoGapANDzeroWalls)

# compile the macros with g++
g++ findthes.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o findthes.exe || exit 1


for sampleNum in 2
do
    for algo in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
    do
	L1Output="~/scratch/EmulatorResults/${InputType[sampleNum]}_JetResults_${AlgoVariations[algo]}.root"
	HistOutput="hist_${InputType[sampleNum]}_${AlgoVariations[algo]}.root"
	Output="results/filerate_${InputType[sampleNum]}_${AlgoVariations[algo]}.root"
	PlotOutputTag="${InputType[sampleNum]}_${AlgoVariations[algo]}"
	./findthes.exe "$L1Output" "$HistOutput" "$Output"|| exit 1
    done
done

