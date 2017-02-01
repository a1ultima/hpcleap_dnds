# Sliding window dN/dS web aplication

## What?

Web service to do the following for available VectorBase (VB) gene IDs:

 1. **dN/dS analysis**: Reveal selective pressures along a protein sequence belonging to the user-specified VB gene id query (e.g. AGAP010815), using an orthologous VB gene id's protein as a reference (e.g. AAEL001802). Selective pressure is measured in dN/dS, as calculated using the pairwise [Nei-Gojobori algorithm][1], which we modified in a way that allows for Gaussian-smoothed dN/dS sliding window output (see image below). Responsiveness was achieved by pickling pre-computed input statistics for the dN/dS outputs, for all possible codon, and amino acid pairs: ~900 x speedup achieved. [Accuracy benchmarks vs. MATLAB's flagship Bioinformatics toolbox sliding window dnds() function.][2] 
 
 2. **Functional Domains Overlay**: Show functional protein domain annotations along the Query sequence, available from VB (see image below). Responsiveness was achieved using REST api calls to VB.
 
 3. **Aggregate data**: Bioinformatics data available for the user-specified VB gene id query (e.g. AGAP010815): General info, Transcript sequences (DNA) and Protein sequences (AA). Responsiveness was achived using REST api calls to VB.

4. **Codon alignment (coming soon)**: Codon alingments of Query and Reference sequences were necessary as input for dN/dS computations, but alignments are not yet shown in their own panel in the web app. Alignments were not optimised for responsiveness, this is currently a performance bottleneck, to be fixed soon. 

## Why?

Combining the two types of information (1. dN/dS analysis, 2. functional domains annotations) could aid in exploring hypotheses concerning ancestral evolutionary selective pressures acting on a protein and it's functional domains. One example insight that can be drawn from the combination of information: 1. and 2., is as follows: if we observe that the majority of functional domains annotated onto a protein sequence, overlap well with dN/dS values ~ 0, it is likely that the protein as a whole is in the process of becoming a pseudogene, and so we can conclude that the corresponding functional domains are for some reason no longer essential to that particular species' survival; thus elucidating evolutionary history of the species-in-question's ancestors.

## Example

- **Query id**: AGAP010815 
- **Reference id**: AAEL001802

![alt-text](https://github.com/a1ultima/hpcleap_dnds/blob/master/py/data/webapp_demo_dnds-and-domains.PNG "demo of dnds and domain panels")





## REQUIREMENTS (tested on):
 - Python 2.7.3
 - Biopython 1.66
 - Firefox 50.1.0
 - Bash
 
## USAGE (Bash terminal):

**1. Start the Python server** 
 - In its own terminal: `bash ./run_py_server.sh`

**2. Open the web-page:**
 - In its own terminal: `xdg-open ./vb-genes.html`

**3. Aggregate sequence data from VectorBase:**
 - Enter a valid VectorBase gene id into the text field (Query)
 - Click "Go!"
 
**4. Retrieve protein functional domains:**
 - Click the "Protein domain annotation (from VectorBase)" panel
 
**5. Compute sliding window dN/dS analysis curve:**
 - Click the "Sliding Window dN/dS analysis" panel
 - Enter a valid VectorBase gene id into the text field (Reference), must be an orthologue to the Query (see: [Nei-Gojobori][1] for help)
 - Click "Calculate dN/dS!"

[1]: https://www.ncbi.nlm.nih.gov/pubmed/3444411
[2]: https://github.com/a1ultima/hpcleap_dnds/blob/master/py/data/benchmarks.txt

# Aknowledgements:
 - **Boilerplate web code, responsive aggregator**: John Kirmitzoglou, Robert MacCallum
 - **Responsive functional domains overlay**: Wenping Lyu
 - **Responsive dN/dS computation, slow alignment:** Andrew Brockman
 
