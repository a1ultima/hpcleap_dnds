#!/usr/bin/python

"""

Calculate dN/dS using Nei Gojobori method with Juke-Cantor's multiple-substitution correction (optional), and whole sequence or sliding window. Tested against MATLAB's dnds().

USAGE: 
    
    Run this from: 

    hpcleap_dnds/  (i.e. ../../)
 
REQUIREMENTS:

    scripts:

        - hpcleap_dnds/py/scripts/changes.py (hpcleap_dnds/py/scripts/)

    data:

        observed_changes.p, potential_changes.p: returned by changes.py
        seq1, seq2: @TODO: make cmd or file input arg

"""
############
# IMPORTS: #
############

import pickle
import math
import changes as codon_pair_data
import warnings
import pdb

# @TODO:REMOVE: \/ and \/\/
# with open('../data/observed_changes_dict.p','rb') as f_observed:
#     changes_observed = pickle.load(f_observed)


# @TODO:REMOVE: \/ and \/\/: the final script cannot use absolut paths, note: the path from which this script is executed is where current working directory is, in this case it should be a .html that is in {root} calling this script residing in {root}/scripts directory
 # changes_observed = pickle.load(open('/home/qiime/Desktop/hpcleap_wp6_compbio/hpcleap_bioinf/data/observed_changes_dict.p','rb'))
# changes_potential= pickle.load(open('/home/qiime/Desktop/hpcleap_wp6_compbio/hpcleap_bioinf/data/potential_changes_dict.p','rb'))


#############
# FUNCTIONS #
#############



def chunks(l, n):
    """ Yield successive n-sized chunks from l. """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def dnds( seq1, seq2, changes_potential, changes_observed, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):
    """ Perform dN/dS analysis, using the 'NG' algoritm, includes both whole sequence or sliding window, and either an approximate or exact multiple-substiution correction method. (@TODO: make sure it actually is exact... it could be
             something else)
            

    ARGS:
        seq1,  a DNA sequence as string of letters, AGTC. Seq1 must be equal in length 
            to, and aligned with, seq2, with gaps trimmed away. @TODO: how on earth can 
            we reliably make this work on the web service?
            
            e.g. seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
        
        seq2,  a DNA sequence similar to seq1 but with differences (substitutions), 
            representing a CDS orthologue of seq1 from a different species. Read 
            description of seq1 for other required similarities to avoid errors.
            
            e.g. seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC'

        changes_potential, a dict, with key=pair of codons tuple, e.g. ('ATG','ATG'), and value=['S':<S>,'N':<N>]. Where <S> is the number of potential synonmyous sites for each codon (averaged between the two codons), and <N> is the same but for non-synonymous sites.

            e.g. changes.potential_changes_dict(...)  (see: ./changes.py)

        changes_observed, @TODO

        msCorrect, a string to toggle between multiple-substitution correction methods:
            "approximate", "exact" (@TODO: make sure it actually is exact... it could be
             something else)
            
            e.g. msCorrect = 'approximate'

        sliding, a boolean to toggle between sliding window analysis (vector of dN/dS values at successive chunks of sequence) or whole sequence analysis (a single 
            dN/dS value for the given pair of input sequences), either: True, False
            
            e.g. sliding = False

        windowLength, an integer specifying the width of the sliding window, measured in DNA basepairs, from 1-to-length(seq1)
            
            e.g. windowLength = 50

    NOTES:

        Sources of formulae:

            http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.html

    """

    # TODO: stop codons to deal with, reject
    # TODO: ambiguous bases to deal with: gaps, Ns, Xs

    # STATS per CODON-PAIR:
    codons_seq1   = [codon for codon in chunks(seq1,3)]  #splits
    codons_seq2   = [codon for codon in chunks(seq2,3)]  #splits
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))==6] # aligned codons are paired into tuples, excess codons are truncated

    # @TODO: the next for loop is extremely innefficient, I should set the structure of the changes_potential and changes_observed dicts to how I want it to look, a priori, i.e. when it's instantiated in changes.py. 
    # OR just remove this chunk and access the data as they come in the args

    changes_all = {'observed':{'S':[],'N':[]},'potential':{'S':[],'N':[]}}

    for pair in codons_paired:
        changes_all['potential']['S'].append(changes_potential[pair]['S'])
        changes_all['potential']['N'].append(changes_potential[pair]['N'])
        changes_all['observed']['S'].append(changes_observed[pair]['S'])
        changes_all['observed']['N'].append(changes_observed[pair]['N'])

    list_S  = changes_all['potential']['S']
    list_Sd = changes_all['observed']['S']

    list_N  = changes_all['potential']['N']
    list_Nd = changes_all['observed']['N']

    if sliding:
        # STATS for each WINDOW seq
        intervals    = range(0,len(codons_paired)-windowLength+1,stepLength)
        windows      = zip(intervals,[i + windowLength - 1 for i in intervals]) 
        window_stats = {}

        for window_i,window in enumerate(windows):
            start = window[0]
            end   = window[1]+1
            window_stats[window] = {    'S':sum(list_S[start:end]),
                                        'Sd':sum(list_Sd[start:end]),
                                        'N': sum(list_N[start:end]),
                                        'Nd':sum(list_Nd[start:end])    }

            pS = window_stats[window]['Sd']/window_stats[window]['S']
            pN = window_stats[window]['Nd']/window_stats[window]['N']

            try:
                if msCorrect=='approximate':
                    #print "cow holy"
                    dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                    dS = -(3./4.)*math.log(1.-(4./3.)*pS)

                # @TODO: what is this commented code? I don't remember...
                # elif msCorrect=='exact':
                #     d=ln(1-p*4/3)/ln(1-3/(4*N))

                else: # msCorrect=='????'  # @TODO: is this the exact one? Or something else?
                    #print "holy cow"
                    dN = pN
                    dS = pS
                window_stats[window]['dNdS'] = dN/dS
            # @TODO: I'm not sure I'm treating the following exceptions in the right way...
            # technically it woud be best to exclude these from downstream analyses? 
            # e.g. missing value/datapoint on a plot of dN/dS (y-axis) vs. window interval (x-axis)
            except ZeroDivisionError:
                warnings.warn("Approximate multiple-substitutions correction cannot be achieved, for window:"+str(window_i)+", dS is zero, leading to a division error when trying dN/dS... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...") # @latest
                window_stats[window]['dNdS'] = float('Inf')
            except ValueError:
                warnings.warn("Approximate multiple-substitutions correction cannot be achieved, for window:"+str(window_i)+",  SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...") # @latest
                window_stats[window]['dNdS'] = float('nan')

        return window_stats
    else:
        # STATS for WHOLE SEQ
        S   = sum(list_S)
        Sd  = sum(list_Sd)
        pS  = Sd/S 
        N   = sum(list_N)
        Nd  = sum(list_Nd)
        pN  = Nd/N

        if msCorrect=='approximate':
            #print "cow holy (whole)"

            if (pS>=3./4.):
                pdb.set_trace()
                raise ValueError("Approximate multiple-substitutions correction cannot be achieved, SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...") 

            if (pN>=3./4.):
                pdb.set_trace()
                raise ValueError("Approximate multiple-substitutions correction cannot be achieved, NON-SYNONYMOUS changes per synonymous site, pN>=3/4, try alternative value for argument: msCorrect...") 

            # @TODO: if the proportion of syn is >3/4 it might break
            try:
                dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
                dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
                dN_dS = dN/dS
            except ValueError:
                warnings.warn("ValueError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...") 
                dN_dS = float("nan")
                return dN_dS
            except ZeroDivisionError:
                warnings.warn("ZeroDiviSionError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...") 
                dN_dS = float('Inf')
                return dN_dS

        else: # @TODO: is this the exact one? Or something else? 
            #print "holy cow (whole)"
            dS  = pS  # i.e. dS = Sd/S
            dN  = pN
            dN_dS = dN/dS

        #print "dN/dS: "+str(dN/dS)

        return dN_dS  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)


# seq1 = 'CTTTTTAACGAAAAAGAGGCAGATGA'
# seq2 = 'ATGGCCCACTTCCCAGGGTTTGGACA'


# dN/dS=0.15270463083614955 with msCorrect='approximate', dN/dS=0.16529268957638238 with 
# msCorrect='exact', MATLAB gives: dN/dS=0.0221=dnds(seq1,seq2, 'geneticCode', 1, 'Method', 'NG'). This is a huge error, an order of magnitude off. 
# seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
# seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC' 
# dnds('ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC', 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC')

# ////


# seq1 = 'ATGGAGCTCAGACTAGAGCGGGGCTGTTGACGTTTGGAGTTGAAAAAATCTATTATACCAATCGGCTTCAACGTGCTCCACGGCAGGCGCCTGACGAGGGGCCCACACCGAGGAAGTAGACTGTTGCACGTTGGGGATAGCGGTAGCTAACTAAGACGCCTGCCACAACAGCAGTATCAAACCCGTACAAAGGGAACATCCACACTTTGGTGAATCGAAGCGCGGCATCAGAATTTCCTTTTGGATACCTGATACAAAGCCCATCGTGGTCCTTAGACTTCGTACACTTACACCTGCACCGCG'
# seq2 = 'ATGCGCATGTGGAATTAGAGGCGAAGTACGATCCCTAGACCGACGTACGATGCAACTGTGTGGATGTGACGAGCTTCTTTTATATGCTTCGCCCGCCGGACCGGCCTCGGCATGGCGTAGCAGTGCACAAGCAAATGACAATTAACCACCGTGTATTCGTTATAACATCAGGCAGTTTAAGTCGGGACAATAGGAGCCGCAATACACAGTTTACCGCATCTTGACCTAACTGACATACTGCCATGGACGACTAGCCATGCCACTGGCTCTTAGATAGCCCGATACAGTGATTATGAAAGGTTT' 


# Note: msCorrect="aproximate"
#dnds_sliding = dnds( seq1,seq2, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )

# # Note: msCorrext="exact"
# dnds_sliding = dnds( seq1,seq2, msCorrect='exact', sliding=True, windowLength=50, stepLength=1 )

# dnds_whole = dnds( seq1,seq2, msCorrect='approximate', sliding=False )

#dnds_sliding

# import pprint
# #pprint.pprint(dnds_sliding)  # dict of dN/dS at each window
# pprint.pprint(dnds_whole)     # dN/dS of the whole seq, for the original test seq1,seq2 we get dN/dS=0.15270463083614955 with msCorrect='approximate', and dN/dS=0.16529268957638238 with msCorrect='exact'

#print(dnds_whole) 

# S  = sum(changes['potential']['S'])
# Sd = sum(changes['observed']['S'])
# pS = Sd/S
# N  = sum(changes['potential']['N'])
# Nd = sum(changes['observed']['N'])
# pN = Nd/N 
# #d   = -(3.0/4.0) * log_e(1 - 4/3 p)
# dN = d()
# dS = d()

if __name__ == "__main__":

    # 
    # TRY CACHED DATA:
    #   Open dictionaries that have cached computationally intensively 
    #   produced results 
    #
    try:
    # LOAD CACHED (fast)

        #
        # Unpickle codonPair-to-statistics dictionaries
        #
        f1 = open('./py/data/observed_changes_dict.p','rb')
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('./py/data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        print "LOADED!!" # @TODO:REMOVE
    except IOError:
    # CREATE NEW (slow)

        #
        # Create dictionary of codon (DNA-triplet, e.g. "ATG") -to- amino 
        #   acid residue (e.g. "M")
        #
        nt_to_aa_dict     = codon_pair_data.geneticCode("standard")

        #@TODO:Urgent:debug:2016-11-28: why does the following two lines lead to silent error? For some reason I HAVE to pickle the data to get it to work... when I just return the dicts directly I get the wrong dN/dS values. The following two lines illustrate this, when uncommented in place of the "# @2:Create" and "# @2:Unpickle" blocks of code...
        # observed_changes  = codon_pair_data.potential_changes_dict(nt_to_aa_dict)
        # potential_changes = codon_pair_data.observed_changes_dict(nt_to_aa_dict)

        #
        # @2:Create the cached codonPair-to-statistics dictionaries, then pickle
        #
        codon_pair_data.potential_changes_dict(nt_to_aa_dict)
        codon_pair_data.observed_changes_dict(nt_to_aa_dict)

        #
        # @2:Unpickle codonPair-to-statistics dictionaries
        #
        f1 = open('./py/data/observed_changes_dict.p','rb') # @TODO:coderedundancey = bad, wrap these lines into a function and call the function
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('./py/data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()
        print "CREATED!!" # @TODO:REMOVE

    #
    #  Calculate dN/dS
    #

    # @TODO: These are either taken as cmd input args or from a tmp file
    s1 = 'CGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAATG'
    s2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAA'


    # @NOTE:uncomment below to achieve dnds of 0.15.. or 0.164 if using exact method, interestingly 
    # dnds_whole = dnds( s1, s2, potential_changes, observed_changes, msCorrect='approximate', sliding=False)

    # @TODO: work with the sliding window version instead of dnds_whole

    dnds_whole = dnds( s1, s2, potential_changes, observed_changes, msCorrect='approximate', sliding=False )

    print "dN/dS: "+str(dnds_whole)
