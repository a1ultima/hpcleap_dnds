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

        - observed_changes.p, potential_changes.p: returned by changes.py
        seq1, seq2: @todo: make cmd or file input arg

"""
############
# IMPORTS: #
############

import pickle
import math
import changes as codon_pair_data
import warnings
import pdb
import time 

start_time = time.time()

# @todo: make a check to ensure the input seqs are divisible by 3 (i.e. n_aa_residues = n_dna_residues/3)

# @todo:REMOVE: \/ and \/\/
# with open('../data/observed_changes_dict.p','rb') as f_observed:
#     changes_observed = pickle.load(f_observed)


# @todo:REMOVE: \/ and \/\/: the final script cannot use absolut paths, note: the path from which this script is executed is where current working directory is, in this case it should be a .html that is in {root} calling this script residing in {root}/scripts directory
 # changes_observed = pickle.load(open('/home/qiime/Desktop/hpcleap_wp6_compbio/hpcleap_bioinf/data/observed_changes_dict.p','rb'))
# changes_potential= pickle.load(open('/home/qiime/Desktop/hpcleap_wp6_compbio/hpcleap_bioinf/data/potential_changes_dict.p','rb'))


#############
# FUNCTIONS #
#############


def dnds( seq1, seq2, changes_potential, changes_observed, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):
    """ Perform dN/dS analysis, using the 'NG' algoritm, includes both whole sequence or sliding window, and either an approximate or exact multiple-substiution correction method. (@todo: make sure it actually is exact... it could be
             something else)
            

    ARGS:
        seq1,  a DNA sequence as string of letters, AGTC. Seq1 must be equal in length 
            to, and aligned with, seq2, with gaps trimmed away. @todo: how on earth can 
            we reliably make this work on the web service?
            
            e.g. seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
        
        seq2,  a DNA sequence similar to seq1 but with differences (substitutions), 
            representing a CDS orthologue of seq1 from a different species. Read 
            description of seq1 for other required similarities to avoid errors.
            
            e.g. seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC'

        changes_potential, a dict, with key=pair of codons tuple, e.g. ('ATG','ATG'), and value=['S':<S>,'N':<N>]. Where <S> is the number of potential synonmyous sites for each codon (averaged between the two codons), and <N> is the same but for non-synonymous sites.

            e.g. changes.potential_changes_dict(...)  (see: ./changes.py)

        changes_observed, @todo

        msCorrect, a string to toggle between multiple-substitution correction methods:
            "approximate", "exact" (@todo: make sure it actually is exact... it could be
             something else)
            
            e.g. msCorrect = 'approximate'

        sliding, a boolean to toggle between sliding window analysis (vector of dN/dS values at successive chunks of sequence) or whole sequence analysis (a single 
            dN/dS value for the given pair of input sequences), either: True, False
            
            e.g. sliding = False

        windowLength, an integer specifying the width of the sliding window, measured in no. of codons in the window to measure dN/dS over, from 1-to-length(seq1)
            
            e.g. windowLength = 50

        stepLength, an integer specifying no. of codons to shift the sliding window with each iteration. If stepLength < windowLength then windows will overlap, overlapping is dealt with prior to plotting (acts to smooth values, averages along the overlaps are taken as a dN/dS value for any codon).

            e.g. stepLength = 1

    NOTES:

        Sources of formulae:

            http://www.megasoftware.net/mega4/WebHelp/part_iv___evolutionary_analysis/computing_evolutionary_distances/distance_models/synonymouse_and_nonsynonymous_substitution_models/hc_nei_gojobori_method.html

    """

    def chunks(l, n):
        """ Yield successive n-sized chunks from l. """
        for i in xrange(0, len(l), n):
            yield l[i:i+n]

    # TODO: stop codons to deal with, reject
    # TODO: ambiguous bases to deal with: gaps, Ns, Xs

    # STATS per CODON-PAIR:
    codons_seq1   = [codon for codon in chunks(seq1,3)]  #splits
    codons_seq2   = [codon for codon in chunks(seq2,3)]  #splits
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))==6] # aligned codons are paired into tuples, excess codons are truncated, @todo: in main example, we lose 5 bps of data

    # @todo: the next for loop is extremely innefficient, I should set the structure of the changes_potential and changes_observed dicts to how I want it to look, a priori, i.e. when it's instantiated in changes.py. 
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


        #pdb.set_trace()

        window_stats = {}

        window_stats_list = []

        #pdb.set_trace()  # @todo: test against matlab's sliding window, also @todo: find out what stepLength does, @todo: try to plot the sliding window version

        for window_i,window in enumerate(windows):

            start = window[0]
            end   = window[1]+1

            #print "from: "+str(start)+"  to: "+str(end)

            window_stats[window] = {    'S':sum(list_S[start:end]),
                                        'Sd':sum(list_Sd[start:end]),
                                        'N': sum(list_N[start:end]),
                                        'Nd':sum(list_Nd[start:end])    }

            pS = window_stats[window]['Sd']/window_stats[window]['S']
            pN = window_stats[window]['Nd']/window_stats[window]['N']

            try:
                if msCorrect=='approximate':
                    dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                    dS = -(3./4.)*math.log(1.-(4./3.)*pS)

                # @todo: what is this commented code? I don't remember...
                # elif msCorrect=='exact':
                #     d=ln(1-p*4/3)/ln(1-3/(4*N))

                else: # msCorrect=='????'  # @todo: is this the exact one? Or something else?
                    dN = pN
                    dS = pS
                window_stats[window]['dNdS'] = dN/dS
                window_stats_list.append(dN/dS)
            # @todo: I'm not sure I'm treating the following exceptions in the right way...
            # technically it woud be best to exclude these from downstream analyses? 
            # e.g. missing value/datapoint on a plot of dN/dS (y-axis) vs. window interval (x-axis)
            except ZeroDivisionError:
                warnings.warn("Approximate multiple-substitutions correction cannot be achieved, for window: "+str(window_i)+", dS is zero, leading to a division error when trying dN/dS... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...")
                window_stats[window]['dNdS'] = float('Inf')
                window_stats_list.append(float('Inf'))
            except ValueError:
                warnings.warn("Approximate multiple-substitutions correction cannot be achieved, for window: "+str(window_i)+",  SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...") # @latest
                window_stats[window]['dNdS'] = float('nan')
                window_stats_list.append(float('nan'))

        return window_stats_list,window_stats  # list of dnds per window interval // dict of dnds, key=(<from #base pair>,<to #base pair>), value=<dN/dS of the window specified in the key>
    else:
        # STATS for WHOLE SEQ
        S   = sum(list_S)
        Sd  = sum(list_Sd)
        pS  = Sd/S 
        N   = sum(list_N)
        Nd  = sum(list_Nd)
        pN  = Nd/N

        try:
            if msCorrect=='approximate':

                if (pS>=3./4.):
                    raise ValueError("Approximate multiple-substitutions correction cannot be achieved, SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 

                if (pN>=3./4.):
                    raise ValueError("Approximate multiple-substitutions correction cannot be achieved, NON-SYNONYMOUS changes per synonymous site, pN>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 

                dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
                dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
                dN_dS = dN/dS

            else: # @todo: is this the exact one? Or something else? 
                
                # @DONE: one day the following three lines of code will error, giving a ZeroDivisionError, this needs to be handled with try
                dS  = pS  # i.e. dS = Sd/S
                dN  = pN
                dN_dS = dN/dS
        except ValueError:
            warnings.warn("ValueError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...") 
            dN_dS = float("nan")
            return dN_dS
        except ZeroDivisionError:
            warnings.warn("ZeroDiviSionError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...") 
            dN_dS = float('Inf')
            return dN_dS

        return dN_dS  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)


def plot_dnds_sliding(dnds_slide_dict):

    """ Plots sliding dN/dS values (y-axis) along the input sequence's aligned codons (x-axis). If the sliding windows overlap (see below), then plot_dnds_sliding will also average the overlapping dN/dS values for each codon. 

        ----     window 1
         ----    window 2
          ----   ...
        ^^^^^^^  take average along each column (codon)

    ARGS:
        dnds_slide_dict,    the output of dnds() if the 'sliding' optional argument is set to True

            e.g. dnds_slide_dict = dnds( s1, s2, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )

    RETURNS:
        None,   a plot is generated and written to: py/data/dnds_sliding_test.png

    """
    import matplotlib.pyplot as plt
    import numpy as np
    #import pickle

    # #
    # # dnds_slide_dict has overlapping windows of dnds calculated, windows overlap by stepLength bps, and are uniformly windowLength wide
    # #
    # with open("py/data/dnds_slide_dict.p","r") as fi:
    #     dnds_slide_dict=pickle.load(fi)

    #window_intervals = sorted(dnds_slide_dict.keys()) # @todo: I dont think sorting the windows will make a difference to final result, but it will make it slower, @todo: test this just in case
    window_intervals = dnds_slide_dict.keys() # @todo: I dont think sorting the windows will make a difference to final result, but it will make it slower, @todo: test this just in case

    #
    # vectorize and find max of the window_intervals
    #
    max_window_position = np.amax(window_intervals) # e.g. 243 @done: a whole order of magnitude faster than: max_window_interval = max(window_intervals, key=lambda x: x[1])
    # @todo: ^ doing equivalent operations on np.array() (instead of list) is much faster, so maybe this means we should  np.arrays() in the rest of the code in dnds.py and changes.py
    # @DONE: Q: would finding the max() be faster than sorting then indexing? A: yes

    #
    # Initialise matrix with NaN, each row is the values for a specific window, these will cascade 
    #       and overlap to various degrees depending on stepLength and windowLength
    #

    # ----     window 1
    #  ----    window 2
    #   ----   ...
    # ^^^^^^^  take average along each column
    # 
    overlap_matrix      = np.empty((len(window_intervals),max_window_position+1))  # @todo: are you sure it's +1? initialize empty matrix, note: entries are not actually NaN, just near-zero
    overlap_matrix[:]   = np.NAN # initiate empty np array with NaN, so later we can mask

    for window_i,window in enumerate(window_intervals):

        start = window[0] # e.g. 0
        end   = window[1] # e.g. 49 

        # in the i-th row, fill all elements from the window[0]-th to window[1]-th with the dN/dS value for this window
        overlap_matrix[window_i,start:end+1] = dnds_slide_dict[window]['dNdS'] # @todo: are you sure it's +1? test, keep in mind for these indices it does -1 for the "to" part

    #
    # Mask all non-finite values, to allow proper plotting
    #
    nan_masker              = ~np.isfinite(overlap_matrix) # boolean matrix, True if element is finite, False if element is Inf or NaN
    overlap_matrix_masked   = np.ma.masked_array(overlap_matrix,mask=nan_masker)

    #
    # Columwise averages (i.e. avg all values in each column, leading to a 1D vector of averages)
    #
    overlap_matrix_avg      = overlap_matrix_masked.mean(axis=0)
    #overlap_matrix_avg.mean() # sanity: 0.15584966052233765 (whole seq average: 0.152307100775) not bad!

    #
    # Plot dN/dS along the input sequences
    #
    plt.plot(overlap_matrix_avg)
    # plt.show()
    plt.savefig("py/data/dnds_sliding_test.png")
    #avg_matrix = overlap_matrix.mean(axis=1)

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
        print "LOADED!!" # @todo:REMOVE
    except IOError:
    # CREATE NEW (slow)

        #
        # Create dictionary of codon (DNA-triplet, e.g. "ATG") -to- amino 
        #   acid residue (e.g. "M")
        #
        nt_to_aa_dict     = codon_pair_data.geneticCode("standard")


        # alternative 1 {{ 

        # #@todo:Urgent:debug:2016-11-28: why does the following two lines (commented) lead to silent error? For some reason I HAVE to pickle the data to get it to work... when I just return the dicts directly I get the wrong dN/dS values. The following two lines illustrate this, when uncommented in place of the "# @2:Create" and "# @2:Unpickle" blocks of code...
        # observed_changes  = codon_pair_data.potential_changes_dict(nt_to_aa_dict)
        # potential_changes = codon_pair_data.observed_changes_dict(nt_to_aa_dict)

        # }} 1 alternative 2 {{

        #
        # @2:Create the cached codonPair-to-statistics dictionaries, then pickle
        #
        codon_pair_data.potential_changes_dict(nt_to_aa_dict)
        codon_pair_data.observed_changes_dict(nt_to_aa_dict)

        #
        # @2:Unpickle codonPair-to-statistics dictionaries
        #
        f1 = open('./py/data/observed_changes_dict.p','rb') # @todo:coderedundancey = bad, wrap these lines into a function and call the function
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('./py/data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()

        # }} alternative 2

        print "CREATED!!" # @todo:REMOVE

    #
    #  Calculate dN/dS
    #

    # @todo: These are either taken as cmd input args or from a tmp file
    s1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAA'
    s2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAA'

    # # @NOTE:uncomment below to achieve dnds of 0.15.. or 0.164 if using exact method, interestingly 
    # dnds_whole = dnds( s1, s2, potential_changes, observed_changes, msCorrect='exact', sliding=False)
    # print "dN/dS: "+str(dnds_whole)

    # @todo: work with the sliding window version instead of dnds_whole

    dnds_slide_list, dnds_slide_dict = dnds( s1, s2, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )

    # with open("./py/data/dnds_slide_dict.p","w") as fo:
    #     pickle.dump(dnds_slide_dict, fo)
    # #print dnds_slide_list 

    #
    # Plot the sliding window values
    #
    plot_dnds_sliding(dnds_slide_dict)

    # @todo: show the user also the whole dnds value somewhere in the webpage

print "Elapsed time: "+str(time.time() - start_time)    


