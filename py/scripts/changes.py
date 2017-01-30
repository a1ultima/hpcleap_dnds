#!/usr/bin/python

"""

@todo: descriptions etc.

"""

# IMPORTS:
import numpy as np 
import pickle
from  itertools import permutations
import copy # @todo: remove the deepcopying?
import pdb

# FUNCTION DEFS:
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    """ Check if two values are almost equal """
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def geneticCode(name):

    """ Dictionary that maps codons to amino acids """ 

    if name == 'standard':
        gc = {  'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', \
                'AGT':'S','ATA':'I','ATC':'I','ATG':'M','ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P', \
                'CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L','CTT':'L','GAA':'E','GAC':'D','GAG':'E', \
                'GAT':'D','GCA':'A','GCC':'A','GCG':'A','GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V', \
                'GTT':'V','TAA':'*','TAC':'Y','TAG':'*','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGA':'*','TGC':'C','TGG':'W', \
                'TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F'  }
    return gc

def potential_changes_dict(nt_to_aa):

    """ Generate a dictionary, with S and N pre-calculated for all 
    possible pairs of codons (key: pair of codons, value: (S,N).

    ARGS:
        nt_to_aa, a dict mapping codons (keys), e.g. 'TTA', to 
            amino-acid letter (values), e.g. 'L'
            
            e.g. geneticCode("standard")

    Notes:

        Sources of formulae:

        http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.htm



    """

    potential_changes = {   'S': {  'AAA':0.0,'AAC':0.0,'AAG':0.0,'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, \
                                    'AGG':0.0, 'AGT':0.0, 'ATA':0.0, 'ATC':0.0, 'ATG':0.0, 'ATT':0.0, 'CAA':0.0, 'CAC':0.0, 'CAG':0.0, 'CAT':0.0, \
                                    'CCA':0.0,'CCC':0.0,'CCG':0.0,'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0, \
                                    'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0,'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0, \
                                    'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0,'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0, \
                                    'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0,'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0},

                            'N': {  'AAA':0.0, 'AAC':0.0, 'AAG':0.0, 'AAT':0.0, 'ACA':0.0, 'ACC':0.0, 'ACG':0.0, 'ACT':0.0, 'AGA':0.0, 'AGC':0.0, 'AGG':0.0, \
                                    'AGT':0.0,'ATA':0.0,'ATC':0.0,'ATG':0.0,'ATT':0.0,'CAA':0.0,'CAC':0.0,'CAG':0.0,'CAT':0.0,'CCA':0.0,'CCC':0.0,'CCG':0.0, \
                                    'CCT':0.0,'CGA':0.0,'CGC':0.0,'CGG':0.0,'CGT':0.0,'CTA':0.0,'CTC':0.0,'CTG':0.0,'CTT':0.0,'GAA':0.0,'GAC':0.0,'GAG':0.0, \
                                    'GAT':0.0,'GCA':0.0,'GCC':0.0,'GCG':0.0,'GCT':0.0,'GGA':0.0,'GGC':0.0,'GGG':0.0,'GGT':0.0,'GTA':0.0,'GTC':0.0,'GTG':0.0, \
                                    'GTT':0.0,'TAA':0.0,'TAC':0.0,'TAG':0.0,'TAT':0.0,'TCA':0.0,'TCC':0.0,'TCG':0.0,'TCT':0.0,'TGA':0.0,'TGC':0.0,'TGG':0.0, \
                                    'TGT':0.0,'TTA':0.0,'TTC':0.0,'TTG':0.0,'TTT':0.0}}   


    # Mutate (substitutions) all possible codons in the given genetic code, and count proportions of mutations that are synonymous and non-synonmyous
    for codon in nt_to_aa.keys():

        # assert (codon not in codon_record)  @DONE: no duplicate entries
        # codon_record.append(codon)  

        # Calculate S and N (i.e. potential synonymous and potential
        # non-synonymous sites) ()

        # i.e. Of possible mutations that can occur on the codon, 
        # what proportion are synonymous (no aa change)?

        # for each bp position in the codon...
        for codon_p in range(0,2+1):

            nts = ['A','G','T','C']  # @DONE: refactor away, A: we can't, since the next line

            nts.remove(codon[codon_p]) # we do not consider self substitutions, e.g. A->A

            # ...and for each nucleotide that the bp can change 
            # into... 
            for nt in nts:

                codon_mutated = list(copy.deepcopy(codon))
                #codon_mutated = codon
                codon_mutated[codon_p] = nt  # mutate the basepair
                codon_mutated = ''.join(codon_mutated)
                
                # ...count how many of them are synonymous.
                if nt_to_aa[codon]==nt_to_aa[codon_mutated]:
                    potential_changes['S'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"
                else:
                    potential_changes['N'][codon]+=1.0/3.0 #@DONE: Q: but why 1/3? to me it should be 1/(3*4), A: since it represents 1/3 of a "site"

            # assert((potential_changes['N'][codon]+potential_changes['S'][codon])==3.0)

    
        # @TEST: N + S == 3.0 per codon, as expected
        #putative_total = potential_changes['S'][codon]+potential_changes['N'][codon]
        # try:
        #     assert isclose(putative_total, 3.0)
        # except AssertionError: 
        #     pdb.set_trace()


        # @DONE:DEBUG:matlabTest: I looped over "N" and "S" dicts spearately assuming keys were alphabetical, wrong, so now I deal with codons as they go
        #potential_changes['N'][codon]=3.0-potential_changes['S'][codon]
        
    # # Calculate proportion of non-synonymous changes, by subtracting away synonymous proportions         
    # for codon in potential_changes['S'].keys():
    #     potential_changes['N'][codon]=3.0-potential_changes['S'][codon]

    codons      = nt_to_aa.keys()
    codonPairs  = list(permutations(codons,2))
    selfies     = [(i,i) for i in codons]
    codonPairs  = codonPairs + selfies 
    
    codonPair_to_potential = {}

    for pair in codonPairs:

        # @DEUG:SOLVED: but still, I don't understand why the two alternatives gives such similar values (non-directional dnds, I think this is causing matlabTest bug? A: turns out not to make much difference, mind-blown)
        
        # ALTERNATIVE 1 (no-directionality seems sensible, albeit ever so slightly more different to MATLAB values) {{

        codon1 = pair[0]
        codon2 = pair[1]
        pn1 = potential_changes['N'][codon1]
        pn2 = potential_changes['N'][codon2]
        ps1 = potential_changes['S'][codon1]
        ps2 = potential_changes['S'][codon2]
        # print str((pn1+pn2)/2.)+" "+str((ps1+ps2)/2.)
        #codonPair_to_potential[pair] = {'N':(pn1+pn2)/2.,'S':(ps1+ps2)/2.}
        codonPair_to_potential[pair] = {'N':(pn1+pn2)/2.,'S':(ps1+ps2)/2.}

        # }} 1 ALTERNATIVE 2 (on one hand this alternative gives closer answer to MATLAB, esp. dnds.dnds( ... msCorrect = 'approximate', ... ), but on the other hand it leads to "directionality", i.e. different dnds values if seq1 and seq2 are swapped in the args {{

        # codon1 = pair[0]
        # #codon2 = pair[1]
        # pn1 = potential_changes['N'][codon1]
        # #pn2 = potential_changes['N'][codon2]
        # ps1 = potential_changes['S'][codon1]
        # #ps2 = potential_changes['S'][codon2]
        # codonPair_to_potential[pair] = {'N':pn1,'S':ps1}

        # }} ALTERNATIVE 2

        # @todo:DEBUG: is this assumption correct? Instead of average, perhaps we need to assume e.g. pn1 is "reference sequence" (ancestor), and only store pn1 rather than avg(pn1,pn2)? given an s1 codon and s2 codon, generate average potential 'N' and 'S'

    # Pickle the output (later used by ./dnds.py) rather than return value, so this needs to only be run once.
    #pdb.set_trace()

    f = open('../data/potential_changes_dict.p','w')
    pickle.dump(codonPair_to_potential,f)
    f.close()

    # @todo: check that for a single pair of codons, i.e. one key in .keys(), gives the right 's' and 'n' value, i.e. potential non-synonymous change and potential synonymous change value

    # @todo: find out what: codonPair_to_potential actually is, then 
    # add descr to func
    #pdb.set_trace()

    #return codonPair_to_potential

def observed_changes_dict(nt_to_aa):

    """ Generate a dictionary, with Sd and Nd pre-calculated for all 
    possible pairs of codons (key: pair of codons, value: (Sd,Nd).

    ARGS:
        nt_to_aa, a dict mapping codons (keys), e.g. 'TTA', to 
            amino-acid letter (values), e.g. 'L'
            
            e.g. geneticCode("standard")

    Notes:

        Sources of formulae:

        http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.html

    """

    codons      = nt_to_aa.keys()
    codonPairs  = list(permutations(codons,2))
    selfies     = [(i,i) for i in codons]
    codonPairs  = codonPairs + selfies 
    
    codonPair_to_observed = {}

    #ALTERNATIVE 1 {{


    # @DONE: Q: I have no idea what the below is trying to do... why can't we just count the differences? Am I trying to simulate the ancestral sequence? A: it's finding all possible "mutational pathways" from one codon to the other, from all possible pathways we can yeild the values we need. 
    # NOTES: commented for now, and trying out the counting of differences.
    for pair in codonPairs:
        codon1 = pair[0]
        codon2 = pair[1]
        indices_to_permute = []

        # Collect the position of the letters (1, 2, 3) where the two sequences differ... 
        for letter_i in range(0,3):
            if not codon1[letter_i] == codon2[letter_i]:
                indices_to_permute.append(letter_i)

        # We now have all the possible mutational pathways, represented as indices 
        permuted_indices = list(permutations(indices_to_permute))
        syn = []
        non = []

        for i,path in enumerate(permuted_indices):
            syn.append(int()) 
            non.append(int()) 

            codon1_path1 = list(codon1) # copies of seqs for 'mutating'
            codon2_path1 = list(codon2)

            for site in path:

                # ALTERNATIVE 1: @comparing mutants to original {{  

                # # @todo:DEBUG: I think I found it! we need to compare all mutated codons to the ORIGINAL codon, not successively mutated codons each mutation!
                # #codon1_past         = ''.join(codon1_path1)
                # codon1_path1[site]  = codon2_path1[site]        # s1 = 'TTT' , s2 = 'ATA'  ==> 'TTT' --> 'ATT' 
                # codon1_path1_str    = ''.join(codon1_path1)
                # #@comparison-step mutants to original
                # if nt_to_aa[codon1_path1_str] == nt_to_aa[codon1]:  # 'TTT --> 'ATT'
                #     syn[i] += 1.0  
                #     non[i] += 0.0
                # else:
                #     syn[i] += 0.0
                #     non[i] += 1.0

                # }} 1 ALTERNATIVE 2: comparing mutants successively (seems to give closer answer to MATLAB, esp. using dnds.dnds( ... , msCorrect='approximate', ... ) ) {{ 

                codon1_past         = ''.join(codon1_path1)
                codon1_path1[site]  = codon2_path1[site]        # s1 = 'TTT' , s2 = 'ATA'  ==> 'TTT' --> 'ATT' 
                codon1_path1        = ''.join(codon1_path1)
                #@comparison-step mutants successively
                if nt_to_aa[codon1_path1] == nt_to_aa[codon1_past]:  # 'TTT --> 'ATT'
                    syn[i] = syn[i] + 1 
                    non[i] = non[i] + 0
                else:
                    syn[i] = syn[i] + 0
                    non[i] = non[i] + 1
                codon1_path1 = list(codon1_path1)

                # }} ALTERNATIVE 2

        try:
            assert isclose(np.mean(syn)+np.mean(non),float(len(path)))
        except AssertionError:
            # pdb.set_trace()
            raise ValueError("Calculations are incorrect, mutation pathways calculation failed...")


        codonPair_to_observed[pair] = {'S':np.mean(syn),'N':np.mean(non)}

    #}} ALTERNATIVE 2 {{ 

    #
    # Count observed differences between two codons, either: synonymous difference, or non-synonymous
    # 

    # @todo:debug: maybe use explicit floats to store counts? 1.0 instead of 1

    # for codon1,codon2 in codonPairs:

    #     codonPair_to_observed[(codon1,codon2)] = {'S':int(),'N':int()}  # None instead of Int() since we want it to fail if we missed it later  


        # determine how many letters differ (if >1 then we need to determine 
        # "mutational pathways", and count the total syn and non-syn changes
        # along all pathways



        # # count no. synonymous and non-synonymous differences for this pair of codons, iterate over each letter in a codon
        # for letter_i in range(0,len(codon1)):

        #     # Do we have a difference?
        #     if codon1[letter_i]==codon2[letter_i]:
        #         pass
        #     else:
        #         # Is the difference synonymous?
        #         if nt_to_aa[codon1]==nt_to_aa[codon2]:
        #             codonPair_to_observed[(codon1,codon2)]['S'] += 1.0
        #         else:
        #             codonPair_to_observed[(codon1,codon2)]['N'] += 1.0

        # #@todo:TEST:ensure we never have more than 3 total differences
        # try:
        #     assert( (codonPair_to_observed[(codon1,codon2)]['N']+codonPair_to_observed[(codon1,codon2)]['S']) <= 3.0 )
        # except AssertionError:
        #     pdb.set_trace()

    #}} ALTERNATIVE 2

    # with open('./py/data/observed_changes_dict.p','wb') as f:
    #     pickle.dump(codonPair_to_observed,f)

    #f2 = open('./py/data/observed_changes_dict.p','wb')
    f2 = open('../data/observed_changes_dict.p','wb')
    pickle.dump(codonPair_to_observed,f2)
    f2.close()    

    #return codonPair_to_observed


if __name__ == "__main__":
    # RUN
    nt_to_aa_d = geneticCode('standard')
    potential_changes_dict(nt_to_aa_d)
    observed_changes_dict(nt_to_aa_d)
    #print('COMPLETE!')
