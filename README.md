# Sliding window dN/dS web aplication

## What?

Web service to do the following for available VectorBase (VB) gene IDs:

 1. **dN/dS analysis**: Reveal selective pressures along a protein sequence belonging to the user-specified VB gene id query (e.g. AGAP010815), using an orthologous VB gene id's protein as a reference (e.g. AAEL001802). Selective pressure is measured in dN/dS, as calculated using the pairwise [Nei-Gojobori algorithm][1], which we modified in a way that allows for Gaussian-smoothed dN/dS sliding window output (see image below).
 
 2. **Functional Domains Overlay**: Show functional protein domain annotations along the Query sequence, available from VB (see image below).

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

`bash ./run_py_server.sh &&
xdg-open ./vb-genes.html`

[1]: https://www.ncbi.nlm.nih.gov/pubmed/3444411

# Aknowledgements:
 - **Boilerplate web code, responsive aggregator**: John Kirmitzoglou, Robert MacCallum
 - **Responsive functional domains overlay**: Wenping Lyu
 - **Responsive dN/dS computation, slow alignment:** Andrew Brockman
 
