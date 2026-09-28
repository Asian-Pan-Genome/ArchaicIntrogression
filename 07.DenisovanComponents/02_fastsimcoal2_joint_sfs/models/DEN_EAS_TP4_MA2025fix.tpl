// TP4-FX: four ordered Denisovan pulses into EAS from one shared ghost donor; no Neanderthal population
4 samples to simulate :
// Demes: 0 AFR, 1 EAS, 2 DENI shared ghost donor, 3 DENS Altai Denisovan
86954 sfspool 0
189407 sfspool 1
5069 sfspool -1
5069 sfspool 2
// Haploid sample sizes and sampling ages; must match the observed folded MSFS
20 0
40 0
0 0
2 2203
// Growth rates
0
-0.0022539867
0
0
// No continuous migration
0
// Historical events, backward in time
9 historical event
$T_PULSE1$ 1 2 $ADMIX1$ 1 0 0
$T_PULSE2$ 1 2 $ADMIX_EVENT2$ 1 0 0
$T_PULSE3$ 1 2 $ADMIX_EVENT3$ 1 0 0
$T_PULSE4$ 1 2 $ADMIX_EVENT4$ 1 0 0
1765 1 1 0 1453 0 0 absoluteResize
2151 1 0 1 1 0 0
3640 0 0 0 38616 0 0 absoluteResize
10599 3 2 1 5069 0 0 absoluteResize
20395 2 0 1 38616 0 0 absoluteResize
// One expected multidimensional SFS locus
1 0
// Per chromosome: number of contiguous linkage blocks
1
// Per block: data type, number of loci, recombination and mutation rates
FREQ 1 0 1.25e-8 OUTEXP
