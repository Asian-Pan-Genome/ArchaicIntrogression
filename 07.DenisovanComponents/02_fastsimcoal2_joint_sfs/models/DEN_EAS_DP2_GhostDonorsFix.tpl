// DP2-GD: two pulses from two divergently related ghost Denisovan donors; no Neanderthal
5 samples to simulate :
// Demes: 0 AFR, 1 EAS, 2 D1 close donor, 3 D2 intermediate donor, 4 DENS Altai Denisovan
86954 sfspool 0
189407 sfspool 1
5069 sfspool -1
5069 sfspool -1
5069 sfspool 2
// Haploid sample sizes and sampling ages
20 0
40 0
0 0
0 0
2 2203
// Growth rates
0
-0.0022539867
0
0
0
// No continuous migration
0
// D1-Altai split ~283 kya; D2 from their ancestor ~307 kya
8 historical event
$T_PULSE1$ 1 2 $ADMIX1$ 1 0 0
$T_PULSE2$ 1 3 $ADMIX_EVENT2$ 1 0 0
1765 1 1 0 1453 0 0 absoluteResize
2151 1 0 1 1 0 0
3640 0 0 0 38616 0 0 absoluteResize
9759 4 2 1 5069 0 0 absoluteResize
10599 3 2 1 5069 0 0 absoluteResize
20395 2 0 1 38616 0 0 absoluteResize
// One expected multidimensional SFS locus
1 0
// Per chromosome: number of contiguous linkage blocks
1
// Per block: data type, number of loci, recombination and mutation rates
FREQ 1 0 1.25e-8 OUTEXP
