#!/bin/bash

ALGOS="twoByTwoANDzeroWalls twoByTwoANDzeroWallsANDsigmaSubtraction"
BINSIZES=("" "_2GeVBin")
SAMPLES="Hydjet502 Hydjet502Dijet30 Hydjet502Dijet80"

BINSIZE=0

g++ findthes.C $(root-config --cflags --libs) -Werror -Wall -Wextra -O2 -o findthes.exe || exit 1

for sample in $SAMPLES
do
    for algo in $ALGOS
    do
	RATEFILE="~/scratch/EmulatorResults/${sample}_JetResults_${algo}${BINSIZES[BINSIZE]}.root"
	HISTFILE="hist_${sample}_${algo}${BINSIZES[BINSIZE]}.root"
	OUTFILE="rate_${sample}_${algo}${BINSIZES[BINSIZE]}"

	echo "Analyzing ${sample}_${algo}${BINSIZES[BINSIZE]}"
	./findthes.exe "$RATEFILE" "$HISTFILE" "$OUTFILE" 0 || exit 1
	echo ""
    done
done
