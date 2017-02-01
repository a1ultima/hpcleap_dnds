
Accuracy & Performance benchmarks vs. MATLAB's Bioinformatics Toolbox dnds() function for (i) whole sequence and (ii) sliding window implementations...

# Accuracy benchmarks

**Concluding remarks:** Whole sequence implementations are comporable w/ MATLAB, slightly less comparable in sliding window implementation, but it is likely due to: (i) the gaussian smoothing heuristic we used to reduce missing values, and (ii) inflating floating point errors when computing many windows as opposed to one single whole-seq window (we could have used double precision, but we used single).

**Input data:** 
 - Query: *AGAP010815_RA*
 - Reference: *AAEL001802_RA* 

## Accuracy  benchmark - Whole sequence dN/dS:

**MATLAB's whole dnds vs. hpcleap_dnds's whole dnds (msCorrect='approximate'):** 

1.0958 vs. 1.04939841193 (Elapsed time: 1.952491045s)

## Accuracy  benchmark - Sliding window dN/dSL:

**MATLAB mean(sliding dnds) vs. andy-wenping mean(sliding dnds): Mean dN/dS over all 50-codon long windows:**

1.7254581782566112 (smoothing?) vs. 1.49715496621 (incl. smoothing, stepLength=1) (Elapsed time: 2.29854488373s) (when taking the mean of all 50 window intervals in MATLAB's sliding dnds())

**Sliding window dN/dS curves**

![alt-text](https://github.com/a1ultima/hpcleap_dnds/blob/master/py/data/matlab_benchmark_sliding.png "sliding window curves vs. MATLAB")
 
 # Performance benchmarks
 
 Coming soon... @TODO
    


