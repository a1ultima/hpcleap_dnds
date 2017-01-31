#!/usr/bin/python

"""

Calculate dN/dS using Nei Gojobori method with Juke-Cantor's multiple-substitution correction (optional), and whole sequence or sliding window. Tested against MATLAB's dnds().

USAGE: 
    
    Run this from command-line, from path (../../): 

    hpcleap_dnds/py/scripts/dnds.py <query sequence> <reference sequence>

        For more information on query and reference sequence command-line arguments see ** in DETAILS section below.

 
REQUIREMENTS:

    scripts:

        - hpcleap_dnds/py/scripts/changes.py (hpcleap_dnds/py/scripts/)

    data:

        - observed_changes.p, potential_changes.p: returned by changes.py
        seq1, seq2: @todo: make cmd or file input arg

DETAILS:



** Query/Reference Input sequences:  both sequences are DNA CDS-transcripts 
of interest (exons-only). To give an example, say the "query" sequence 
corresponds to the Anopheles gambiae TEP1 protein's transcript id 
(VectorBase: AGAP010815), and the other "reference" sequence corresponds 
to the orthologous Aedes aegypti TEP1 protein's trancript id (VectorBase: 
AAEL001802). The orthologous reference sequence is required, along with 
the query sequence of interest, in order to estimate dN/dS values (sliding
 or whole), which requires simulations/approximations of ancestral DNA 
substitution mutations (synonymous, non-synonymous changes).

In the main webpage: vg-genes.html, the query sequence is what the user 
will request aggregate informations from (after clicking "GO!"), whereas 
the orthologous reference sequence is only used to calculate dnds.py, 
without the aggregate information requested.


"""
############
# IMPORTS: #
############

# default modules
import sys
import pickle
import math
import warnings

#import matplotlib.pyplot as plt 
import numpy as np 
# andy-developed modules: imported from py/scripts/ (i.e. "." relative to where this file is)
import changes as codon_pair_data  # /hpcleap_dnds/py/scripts/changes.py
import align as align_then_trim    # /hpcleap_dnds/py/scripts/align.py

# just for testing, @todo: remove testing imports
import pdb
import time 
start_time = time.time()  # @time

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

    # todo: stop codons to deal with, reject
    # todo: ambiguous bases to deal with: 
        # gaps,  @done 
        # Ns,    @todo
        # Xs     @todo

    warning_count = 0

    # STATS per CODON-PAIR:
    codons_seq1   = [codon for codon in chunks(seq1,3)]  #splits
    codons_seq2   = [codon for codon in chunks(seq2,3)]  
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))==6] # aligned codons are paired into tuples, excess codons are truncated, @todo: in main example, we lose 5 bps of data

    # @done: the next for loop is extremely innefficient, I should set the structure of the changes_potential and changes_observed dicts to how I want it to look, a priori, i.e. when it's instantiated in changes.py. 
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

        #window_stats_list = []

        # @done: test against matlab's sliding window, also @todo: find out what stepLength does, @todo: try to plot the sliding window version

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
                    dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                    dS = -(3./4.)*math.log(1.-(4./3.)*pS)

                # @todo: what is this commented code? I don't remember...
                # elif msCorrect=='exact':
                #     d=ln(1-p*4/3)/ln(1-3/(4*N))

                else: # msCorrect=='????'  # @todo: is this the exact one? Or something else?
                    dN = pN
                    dS = pS
                window_stats[window]['dNdS'] = dN/dS
            # @todo: I'm not sure I'm treating the following exceptions in the right way...
            # technically it woud be best to exclude these from downstream analyses? 
            # e.g. missing value/datapoint on a plot of dN/dS (y-axis) vs. window interval (x-axis)
            except ZeroDivisionError:
                warning_count += 1
                #warn_msg = "Query and Reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, for window: "+str(window_i)+", dS is zero, leading to a division error when trying dN/dS... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...\n"   # # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                #warnings.warn(warn_msg)  # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                window_stats[window]['dNdS'] = float('Inf')
            except ValueError:
                warning_count += 1
                #warn_msg="Query and Reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, for window: "+str(window_i)+",  SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined... try alternative value for argument: msCorrect (e.g. 'exact') OR alternative value for argument: windowLength (e.g. "+str(windowLength+20)+") ...\n"  # # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                #warnings.warn(warn_msg)  # @TODO: uncomment for verbose warning message prints // @ANDY-2017-01-30
                window_stats[window]['dNdS'] = float('nan')

        return window_stats,warning_count  # list of dnds per window interval // dict of dnds, key=(<from #base pair>,<to #base pair>), value=<dN/dS of the window specified in the key>
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
                    raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, SYNONYMOUS changes per synonymous site, pS>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 

                if (pN>=3./4.):
                    raise ValueError("Query and reference sequences are too divergent. Approximate multiple-substitutions correction cannot be achieved, NON-SYNONYMOUS changes per synonymous site, pN>=3/4, log() operation will yeild return undefined. Try alternative value for argument: msCorrect (e.g. 'exact')...") 

                dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
                dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
                dN_dS = dN/dS

            else: # @todo: is this the exact one? Or something else? 
                
                # @DONE: one day the following three lines of code will error, giving a ZeroDivisionError, this needs to be handled with try
                dS = pS  # i.e. dS = Sd/S
                dN = pN
                dN_dS = dN/dS
        except ValueError:
            warning_count += 1
            warnings.warn("Query and reference sequencea are too divergent. ValueError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n") 
            dN_dS = float("nan")
        except ZeroDivisionError:
            warning_count += 1
            warnings.warn("Query and reference sequences are too divergent. ZeroDiviSionError: Approximate multiple-substitutions correction cannot be achieved: UNKNOWN reason, probably due to illegal numbers in a log() function...\n") 
            dN_dS = float('Inf')

        return dN_dS, warning_count  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)


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
    # plt.plot(overlap_matrix_avg)
    # plt.show()
    # plt.savefig("py/data/dnds_sliding_test.png")
    # avg_matrix = overlap_matrix.mean(axis=1)
    # print overlap_matrix_avg
    return list(overlap_matrix_avg), overlap_matrix_avg.mean()


# @ANDY:code redundancy, can uncomment if we need to load this whole pipeline as a function

def dnds_pipeline(qry_seq_in, ref_seq_in):

    """ Runs the whole dnds pipeline for sliding windows, returns the smoothed vector
    of dN/dS values along the codons of the query sequence.

    Essentially the same code as: if __name__ == "__main__": ...


    """

    # 
    # TRY CACHED DATA:
    #   Open dictionaries that have cached computationally intensively 
    #   produced results, these are all possible codon pairs with cached statistics for all possible pairs regardless of user input seqs 

    print("\t============================================================================") # @todo:REMOVE
    print("\tdN/dS sliding analysis: pre-cached statistics for all possible codon pairs...")
    print("\t============================================================================") # @todo:REMOVE

    #
    try:       
    #if (os.path.exists("./py/data/observed_changes_dict.p") and os.path.exists("./py/data/potential_changes_dict.p")):
    # LOAD CACHED (fast)

        #
        # Unpickle codonPair-to-statistics dictionaries
        #
        f1 = open('../data/observed_changes_dict.p','rb')
        observed_changes  = pickle.load(f1)
        f1.close()

        f2 = open('../data/potential_changes_dict.p','rb')
        potential_changes = pickle.load(f2)
        f2.close()

        print("\t\tLOADED!") # @todo:REMOVE
        
    except IOError:
    #else:
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

        #}} 1 alternative 2 {{

        #
        # @2:Create the cached codonPair-to-statistics dictionaries, then pickle
        #
        codon_pair_data.potential_changes_dict(nt_to_aa_dict)
        codon_pair_data.observed_changes_dict(nt_to_aa_dict)

        #
        # @2:Unpickle codonPair-to-statistics dictionaries
        #
        #f1 = open('./py/data/observed_changes_dict.p','rb') # @todo:coderedundancey = bad, wrap these lines into a function and call the function
        f1 = open('../data/observed_changes_dict.p','r') # @todo:coderedundancey = bad, wrap these lines into a function and call the function
        observed_changes  = pickle.load(f1)
        f1.close()

        #f2 = open('./py/data/potential_changes_dict.p','rb')
        f2 = open('../data/potential_changes_dict.p','r')
        potential_changes = pickle.load(f2)
        f2.close()

        # }} alternative 2

        print("\t\tCREATED!") # @todo:REMOVE


    #
    #  Calculate dN/dS
    #

    # @todo: 
    #   Q:These are either taken as cmd input args or from a tmp file?
    #   A: No, we take them from the output of align.align_query_vs_reference() and align.trim_gaps_from_aligned_seqs() (align. imported as align_then_trim.)


    # alternative 1 {{

    # Aligned and Trimmed HTLV-1 vs. STLV-1 proteins (MATLAB example to benhmark dnds calculations)
    # s1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAA'
    # s2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAA'


    # }} 1 alternative 2  {{  # @todo: debug:ever since the .p loading bug, 
    #                               which is still un-resolved (see: "@todo:Urgent:debug:2016-11-28")
    #                               I think we should test a load( .p ) version just in case and check

    #
    # AGAP010815 (gambiae)
    #     https://www.vectorbase.org/Anopheles_gambiae/Gene/Sequence?db=core;g=AGAP010815;r=3L:11202091-11206882;t=AGAP010815-RA
    qry_seq_raw = qry_seq_in
    # s1_AGAP010815_RA = 'ATGTGGCAGTTCATAAGGTCACGAATATTAACGGTGATAATCTTCATAGGTGCTGCTCATGGGCTACTGGTTGTGGGTCCGAAATTTATACGGGCCAACCAGGAATACACTCTGGTGATCAGCAACTTTAACTCACAGCTAAGCAAAGTGGACCTGCTGTTAAAACTGGAAGGCGAAACTGATAATGGTTTAAGCGTTCTGAACGTTACCAAGATGGTTGACGTGCGACGTAATATGAACCGAATGATCAACTTCAATATGCCTGAGGATCTGACGGCTGGAAACTACAAAATAACTATCGATGGACAGCGTGGCTTCAGCTTTCACAAGGAGGCAGAGCTGGTGTATCTCAGCAAATCGATATCGGGGCTAATACAGGTCGATAAGCCCGTATTTAAACCTGGGGATACGGTGAACTTCCGTGTGATCGTGCTGGACACGGAGCTGAAACCGCCGGCGAGGGTCAAGTCGGTTTATGTAACTATACGAGATCCTCAGCGCAATGTGATTCGCAAATGGTCCACGGCAAAACTGTATGCCGGTGTGTTCGAGAGCGATCTACAGATAGCGCCTACTCCAATGCTCGGGGTCTGGAATATCTCGGTGGAGGTGGAAGGAGAAGAGCTTGTGTCAAAGACGTTTGAGGTGAAGGAGTACGTGTTGTCAACGTTCGACGTGCAGGTCATGCCATCGGTGATTCCACTGGAAGAGCATCAAGCTGTGAATCTTACAATCGAAGCGAACTATCACTTTGGTAAGCCAGTGCAAGGAGTGGCCAAGGTGGAGCTGTACCTAGACGACGATAAGCTAAAACTGAAAAAAGAGCTGACTGTGTACGGAAAGGGCCAGGTAGAGTTGCGCTTTGACAATTTTGCAATGGATGCGGATCAGCAGGATGTACCAGTGAAGGTGTCGTTCGTCGAGCAGTACACAAATCGTACGGTGGTCAAACAGTCACAAATCACGGTATATAGGTATGCGTACCGAGTAGAGTTGATAAAAGAGAGTCCACAGTTTCGTCCGGGACTCCCGTTCAAATGTGCGCTTCAGTTTACACACCATGATGGAACACCGGCTAAAGGCATTAGCGGTAAGGTAGAGGTATCCGATGTACGATTCGAAACGACAACAACGAGTGATAACGATGGATTGATTAAGCTCGAGCTGCAACCAAGTGAGGGTACTGAACAACTCAGTATTCACTTCAATGCTGTTGATGGATTCTTTTTTTATGAAGATGTGAATAAGGTAGAAACGGTTACGGATGCGTATATTAAACTGGAGCTGAAATCACCGCATCAAACGGAACAAATTGATGCGTTTCATGGTGACGTGCACGGAGCGCATGACATTCTTCGTGTACTATGTCATGTCAAAGGGCAATATCATCGATGCAGGATTCATGCGACCCAACAAGCAACCGAAGTACCTGTTGCAGCTGAACGCAACAGAAAAGATGATTCCGAGGGCGAAAATTCTCATCGCTACCGTAGCGGGCCGCACGGTGGTGTACGACTTCGCAGACCTCGATTTCCAAGAGCTTCGCAATAATTTTGATTTAAGCATTGACGAGCAAGAGATCAAGCCGGGACGACAAATCGAGCTGAGCATGTCTGGACGCCCAGGAGCGTACGTTGGGCTGGCCGCGTATGACAAAGCCTTGCTGCTTTTCAACAAGAACCACGACCTGTTCTGGGAGGACATTGGGCAGGTGTTTGATGGGTTCCATGCAATCAATGAGAACGAGTTTGACATATTCCACAGCTTGGGTCTGTTCGCCAGGACATTGGACGATATCTTGTTCGACAGTGCAAATGAAAAGACGGGGCGTAATGCACTGCAGTCAGGCAAGCCGATCGGCAAGCTGGTGTCGTATCGGACGAACTTCCAGGAATCGTGGTTGTGGAAAAATGTTTCCATCGGACGATCGGGAAGTCGCAAGTTGATCGAGGTAGTACCGGACACGACCACCTCCTGGTATCTGACGGGCTTCTCGATCGATCCCGTGTACGGGTTGGGTATCATCAAGAAGCCAATCCAGTTCACAACAGTCCAGCCGTTCTACATCGTAGAGAACTTACCATATTCAATCAAACGAGGCGAAGCGGTTGTGTTGCAGTTTACGCTGTTCAACAACCTTGGAGCGGAGTATATAGCCGATGTGACGCTGTACAATGTGGCCAACCAGACCGAGTTCGTCGGACGTCCAAATACGGATCTCAGCTACACCAAATCCGTGAGCGTTCCTCCAAAAGTTGGTGTGCCAATCTCGTTCCTCATCAAGGCCCGCAAGCTCGGCGAGATGGCGGTTCGTGTAAAGGCTTCGATAATGCTGGGACACGAAACGGACGCCCTGGAAAAGGTAATACGGGTGATGCCTGAAAGTTTGGTGCAGCCGAGAATGGATACACGCTTTTTCTGCTTCGACGATCACAAAAATCAAACGTTTCCGATCAACTTGGACATCAACAAGAAGGCCGACAGTGGATCGACAAAGATTGAGTTTCGACTAAATCCCAATTTGTTGACCACGGTCATCAAGAACCTGGACCATCTTCTCGGCGTTCCGACGGGATGTGGTGAGCAGAATATGGTCAAATTTGTTCCCAACATTTTGGTACTGGATTATTTGCATGCCATCGGGTCGAAAGAACAGCATCTAATCGACAAAGCTACGAATTTGTTGCGTCAAGGATATCAAAACCAGATGCGCTACCGTCAGACGGATGGTTCATTTGGTTTGTGGGAGACTACTAATGGTAGCGTGTTTCTCACCGCGTTCGTTGGCACATCGATGCAAACTGCAGTAAAATACATAAGCGATATTGATGCAGCAATGGTGGAGAAGGCATTGGATTGGTTAGCCTCGAAGCAGCATTTCTCGGGACGGTTTGACAAGGCCGGTGCAGAGTATCACAAAGAAATGCAAGGAGGGTTGCGCAATGGTGTGGCCCTCACATCATATGTGTTGATGGCATTGCTGGAGAATGACATTGCCAAAGCAAAGCACGCAGAGGTGATTCAAAAAGGAATGACCTATCTGAGCAATCAGTTTGGATCCATCAACAATGCATACGACCTATCGATAGCAACCTACGCGATGATGTTGAACGGACACACCATGAAGGAGGAGGCACTCAATAAGCTGATTGATATGTCTTTCATTGATGCTGATAAAAACGAACGGTTCTGGAACACAACGAATCCAATAGAAACCACCGCATATGCTCTGCTGTCGTTTGTGATGGCCGAGAAGTACACAGACGGTATACCGGTCATGAATTGGTTGGTGAATCAACGTTACGTTACCGGTAGCTTTCCGAGCACGCAAGACACGTTTGTGGGGCTGAAAGCGCTGACCAAAATGGCGGAAAAGATATCTCCGTCCCGAAACGACTACACCGTTCAACTGAAGTACAAGAAGAGTGCAAAATACTTCAAAATAAACTCGGAGCAAATTGATGTGGAAAACTTCGTGGATATACCGGAGGACACAAAAAAGCTCGAGATCAATGTGGGGGGCATTGGATTTGGGTTGTTAGAGGTGGTTTATCAATTTAATTTGAATCTCGTCAACTTTGAGAATAGATTCCAACTAGACCTGGAGAAACAGAACACAGGCTCTGACTACGAGCTGAGGCTGAAGGTCTGTGCCAGCTACATACCCCAGCTGACCGACAGACGATCGAACATGGCACTGATTGAGGTAACCTTACCGAGCGGTTACGTGGTTGATCGCAATCCGATCAGCGAGCAGACGAAGGTGAATCCGATTCAGAAAACTGAAATCCGTTACGGTGGCACTTCAGTCGTTTTATACTACGACAATATGGGCAGCGAGCGTAACTGTTTCACCCTGACCGCGTACAGACGCTTTAAGGTCGCATTGAAGCGTCCAGCGTATGTGGTTGTGTATGATTATTATAATACAAATCTGAACGCCATCAAAGTGTACGAAGTGGACAAGCAGAATTTGTGCGAAATCTGTGACGAAGAAGACTGTCCTGCAGAGTGCAAAAAATAG'
    #
    # AAEL001802 (aegypti)
    #     https://www.vectorbase.org/Aedes_aegypti/Gene/Sequence?db=core;g=AAEL001802;r=supercont1.43:685886-717122;t=AAEL001802-RA
    ref_seq_raw = ref_seq_in
    #s2_AAEL001802_RA = 'ATGTCGGTATTCATACAAACGGACAAACCGGTGTATACCCCGGGAGATCTGATACGTTTTCGGGTAATCGTGGTGGATGCTGACACTAGACCTGTGACTAGTATTAAAACGGTAAATATAGCGATCGACGATTCTGCAAAAAATTCCATTCGAAAGTGGCCTTATGCCAAGTTGTTAAACGGCATCTTTGAGTCACAAGTGCAATTAGCTTCTTCGCCTGTTCTTGGCACCTGGATTATCAACGTAACAGCTTCCGACGACATCATTGTCACCAAACAGATAGAAGTTAAGGAATATGTGTTGCCAAAATTTTTCGTGAAAGTTTACCCTTCGGAGGTTCTATTGGGGAAAAATAAGAAGGTTTCTCTTACCTTAGATGCCTATTACACGTTCAAAGAACCCGTCGACGGCAATTACAAAGTTGAGTTATTTTTGGACCATACCAAGAGAAAGCCTGACTTCATAAAAAGTGATCGAATCACCGGTAAAACATCACTTGAGTTTCAATTGAAAAATGAAGTAGACATTGATGGCGACGAGCAGTACACTGATGTCACGGTTGAAGTTGAAGTTGTCGAGGCATTTTCTAATCGCACAGTTAGTATAACTGAGAATATTCCGATTTATCGTCAGCCTTATACCGTGACCCTTCTTCCATCTGCACCATCATTTCGACCAGGAGTTCCATTCAATGTACAAATAGTTGTGAAAGATCAGCTTGGACACCCTCCTGCCGAAGAAAAAGCGGCATCAATTGACCTTACTGTAGAGTTCCATTTGCCCATTGACAGTGACACCAAATCTATCACTGTAGATCTGGACGAGAAAGGAACAGGTCAGCTCACATTAGAGCCCCGCCCAGACGCCCAAGAACTGAAAGTGAACGCTACATATGACTCTCAACAATACGATGTAATTCACGATCCGATACATGGTTTCAGTTCGCAAAGTAAGCAGTACATCACAGTAACTCTGAATCCAAAATACTATAACAACATTAAAGTCGATAAGGACATCGTACTGGACATCTCCTGCACTGAAACAATGACGCACTTCTCGTACATCGTTGTCACCAGAGGAAACATAGTGGAAGCATCGAACGTTCCTGTCAGGATAAAAAAGAAACATTCTCTGAGATTGAAAATGACTTCAAAAATGTCTCCGGAGTCGAGGCTTCTAGTGTACTATACAAACAGGGAGTATCTCATCTTTGATGATATTGAGCTGAAGTTCGATTCGTTCAACAACGACTTCAAATTCGATTTGAACGATGATGAGTATTTTCCAGGGCAATCAGTTTATATCGATGTATACGCTTCAAAGGATTCATACGTTGCGTTCAGTGGAATCGATGAAAGTGTACTCCTGGTAGGCAAAGAGCGCCATGACTTCAACAAAGGAGATGTGCTCAAGGAACTCGCTCTTTACGGAGCAACAAATGATGCCGAGTTTGACTTGTTCCACGTAAGTTTCATGTCAAATGGTTTGATTATTCCAGTTAATGTATCTGTAACTCGCTCACAGAATGCACGATTTGGTACTCTACTAGGAAGGACTAGGCAGCAAGCGATTGAAATTCGAACTCAATTCCTAGAATCCTGGTTATGGAAATCCTTTTCCATGGATGGTCGAAACAACTTCAAAGCAATAGAAGACTCGGTTCCGGATACTATTACAACGTATCACGTGTCAGGATTTGCTTTAAGTCCAACACTAGGTCTTGGAGTAATCCAACAACCAGTGAGTTTCACCGTTCGTAAAAAATTCTACTTGGTTGCAAATTTGCCTTACTCGATCAAACGGGGTGAAGTGGCGTTGATTCAGGTTACCGTCTTCAACTTCCTAGGAAGCAGCATAACAACCGATGTGACGCTGTTCAATAAACGCGATGAAATTGAGTTTGTCGAGAATGCATCCACTAATAATACACATCGAACAAAGGCGGTAATTGTCCCGAATAACAATGGAAAATCTGTATCATTTATGGTGAAAGCAAAGAAATTAGGACAGATTGCGATCAAATTCCAGGCGGTAAACCTGCTGGAAACGGATGCATTGGAGCACATGTTACGAGTAACCCCAGAGAGCCATCGCTATGAGAAAAATGTAGCTCGATTCGTTGAGCTACCAAAGTTTGAGACGCAAACTTTCGATGTGAAGCTGGACATTCCCAAAAATATCGACGAGGGTTCTGCTCAAATCAAATTCACGTTAGACCCGGACATTTTGGGAACAGCCATCAGCAACCTAGACGGGTTGATCCGGAAACCCTTTGGATGTGGCGAACAAAATATGCTCCATTTTGTGCCAAATATAGTCGTTTTGGATTATCTTAACGAAACCAACACAGCGGCAGAAGATGTGAGGACCAAAGCGATAAATTTTCTTAGCAGCGGATATCAAAACCAGCTACGCTACAAACGTTCGGATGGGGCCTTCAGTGTCTGGGGACAATCGTATGCTGGCAGTACATTTTTGACGGCCTTTGTGGCGAAATCATTCAAAATAGCAGCCAAATACATTCAGGTGGATAAGTCTATAGTAGACGCGGCATTCGACTGGTTAGTGAAACAACAACAATCAGATGGGCGGTTCCCAGAAGTGGGGCAAGTATTCCAAGCAGATATGCAGGGTGGGCTTCGTAATAACGGTTTTGCGCTTACCGCGTATGTTCTGATCGCTTTTGCTGAAAATAAGGAAGTATACAGAAAATACCAATCACAACTGAACAAAACTACTAACTTCATAGCAGATAGACTTGCTAATATGGAGAATCCATACGACCTCTCGCTGTCCACTTATGCGTTGATGCTAACAAATCATGGCAAGCGCACCGAGTTTCTTCACAAATTAGTCGAAAAGTCGATATTTGACCGCAATCAAACTGAGAGATATTGGGACAGCAAACCAGTTGATATTGAAGTTGCTGGATATGCTCTATTGTCATACGTAGCTGCCGGTAAATTATTGGATGCAACGCCTATCATGCGGTGGCTCAACAAGCAGCGTTATGGTCTCGGAGGCTATCCTGGAACTCAGGAAACATTCGTTGGATTGAAAGCATTGGCAACGTTCGCTGCAAATGTAACTAGTAGGAGAAACGAATATACTGTAAGGATATTCTACGAACCAAATGGTCGACGAACATTCGACGTACACATGCACAATTCGTTTAATATTCAAGAGCTTGACATTCCTAATAACATCAGAAAAATGAAGGTGGAAGTTGAAGGCATCGGCAGAGGCTTCTTCCAAGTGGCATATCAGTACTATCAAAATATGCAGGTGGCTAAGCCCAGTTTCAGCATTACAATTAATCAGCTTAACACCACGACGGAACACATGCAGCAATTGGACGTGTGTGTGAAATACATACCAAAAGAGGCTTATCAAAAATCGAATATGGCTTTGGTGGAAATATTCTTGCCTAGTGGGCTTGTAGCAGACTCAGATGCCATTACGGACAAGACTGGAGGAATTCGAAGAATTGAAAGACGTTTTTCGGACACCTCAGTAGTTATATATTATGATAATTTGGACCCCGAAGACAAGTGCTTCCGAGTGACTGCTTATCGTCGGTATAAAATTGCATTGCATTTGCCATCATATATTATAGTTTATGATTATTATAATTTTGAGCGCTTTGCCATTCAAAAGTACGAAGGAAAGGTGCTGCAGCTCTGCGATATTTGTGAAGACGAGGACTGCGAAACTTTATCATGTCAAAATAGCTCGAAATTGGCAATAATGTAA'
    

    ###########################
    # Benchmarking vs. MATLAB #
    ###########################

    print("\t=========================================")
    print("\tProcessing heuristic codon alignments....")
    print("\t=========================================")

    # @done:@testing:benchmark vs. matlab: @todo: place in README these details
    #   using s1_AGAP010815_RA and s2_AAEL001802_RA (after aligning/trimming)
    # MATLAB whole dnds vs. andy-wenping whole dnds (msCorrect='approximate'): dN/dS: 
    #   1.0958 vs. 1.04939841193 (Elapsed time: 1.952491045, Total warnings: 0), respectively
    # MATLAB mean(sliding dnds) vs. andy-wenping mean(sliding dnds): Mean dN/dS over all 50-codon long windows: 
    #   1.7254581782566112 (smoothing?) vs. 1.49715496621 (incl. smoothing, stepLength=1) (Elapsed time: 2.29854488373, Total warnings: 173) (when taking the mean of all 50 window intervals in MATLAB's sliding dnds()
    # MATLAB sliding dnds plot vs. andy-wenping sliding dnds plot: similar sliding dnds plots: 
    #   see:  D:\Dropbox\00_HPC-LEAP\WP5 Thematic Cyprus Computational Biology\giannis_database

    # @benchmark:matlab's dnds()
    # dN/dS=0.15270463083614955 with msCorrect='approximate', dN/dS=0.16529268957638238 with 
    # msCorrect='exact', MATLAB gives: dN/dS=0.0221=dnds(seq1,seq2, 'geneticCode', 1, 'Method', 'NG'). This is a huge error, an order of magnitude off. 
    # seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
    # seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC' 
    # dnds('ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC', 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC')

    #
    # ALIGN REF vs. QUERY seq, then TRIM GAPS: 
    #     Next we want to remove all characters (trim) in each of the sequences that is aligned to a gap position**** in the top alignment
    #

    #  Alignment:
    #     ATGTGC----TAA  qry_seq
    #     ATG--A--TTTGA  ref_seq
    #        ^^ ^^^^     gap positions****
    #  Trimming:
    #     ATGCTAA        qry_seq (after trimming)
    #     ATGATGA        ref_seq (after trimming)

    # alignment scores are calculated from these two parameters, best scoring alignments have few gaps and gaps are small
    aln_gap_open = -10      
    aln_gap_extend = -0.5 # @todo: make cmd args later?

    qry_seq_aln, ref_seq_aln                          = align_then_trim.align_query_vs_reference(qry_seq_raw, ref_seq_raw, aln_gap_open, aln_gap_extend)
    qry_seq_trimmed, ref_seq_trimmed, qry_seq_indices = align_then_trim.trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )

    # @Note: benhamrking: it is qry_seq_trimmed and ref_seq_trimmed that should be used as input to MATLAB's dnds() function, for benchmarking

    print("\t\tCOMPLETE!")

    #}} alternative 2

    # ///

    # alternative 1 {{

    print("\t===========================")
    print("\tSliding window analysis....")
    print("\t===========================")

    # # @NOTE:uncomment below to achieve dnds of 0.15.. or 0.164 if using exact method, interestingly 
    # dnds_whole, warning_count = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=False)
    # print "dN/dS: "+str(dnds_whole)

    # }} 1 alternative 2 {{

    # @DONE: work with the sliding window version instead of dnds_whole
    dnds_slide_dict, warning_count = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=25, stepLength=1 )
    #
    # Plot the sliding window values
    #
    dnds_sliding_vec, dnds_sliding_mean = plot_dnds_sliding(dnds_slide_dict)

    print("\t\tCOMPLETE!")


    # 
    # Summary statistics
    #
    
    print("\t===============================================")
    print("\tElapsed time: "+str(time.time() - start_time))   # @time
    print("\tTotal warnings (missing values): "+str(warning_count))
    print "\tAvg. dN/dS over all windows: "+str(dnds_sliding_mean)
    print("\t===============================================")
    # @todo: show the user also the whole dnds value somewhere in the webpage    
    print("\t\tdN/dS Sliding Analysis: All Jobs Complete!")

    return dnds_sliding_vec, qry_seq_indices



########
# MAIN # (still testing, cmd-args-mode not available yet)
########
if __name__ == "__main__":

    # 
    # TRY CACHED DATA:
    #   Open dictionaries that have cached computationally intensively 
    #   produced results 
    #
    try:       
    #if (os.path.exists("./py/data/observed_changes_dict.p") and os.path.exists("./py/data/potential_changes_dict.p")):
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
    #else:
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

        #}} 1 alternative 2 {{

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

    # @todo: 
    #   Q:These are either taken as cmd input args or from a tmp file?
    #   A: No, we take them from the output of align.align_query_vs_reference() and align.trim_gaps_from_aligned_seqs() (align. imported as align_then_trim.)


    # alternative 1 {{

    # Aligned and Trimmed HTLV-1 vs. STLV-1 proteins (MATLAB example to benhmark dnds calculations)
    # s1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAA'
    # s2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAA'


    # }} 1 alternative 2  {{  # @todo: debug:ever since the .p loading bug, 
    #                               which is still un-resolved (see: "@todo:Urgent:debug:2016-11-28")
    #                               I think we should test a load( .p ) version just in case and check

    #
    # AGAP010815 (gambiae)
    #     https://www.vectorbase.org/Anopheles_gambiae/Gene/Sequence?db=core;g=AGAP010815;r=3L:11202091-11206882;t=AGAP010815-RA
    qry_seq_raw = sys.argv[1]
    # OBP7 gambiae: v
    # qry_seq_raw = "ATGTGTGAATATTCGAATACGCGCAACAAGATGAGCAACCTGGTCGTCGTCCTCGTCCTGCTGACGATGTACATTGTGCTTTCGGCCCCATTCGAAATACCGGACCGGTACAAAAAGCCGGCTAAAATGTTGCACGAAATTTGTATCGCCGAGTCGGGCGCCTCGGAGGAGCAGCTGCGCACCTGTCTCGATGGAACCGTACCGACAGCTCCGGCCGCCAAGTGCTACATCCACTGCCTGTTCGACAAGATCGACGTGGTGGACGAGGCGACTGGGCGCATCCTGCTCGACCGACTGCTTTACATCATCCCGGACGACGTGAAGGCAGCGGTGGACCATTTAACGCGCGAATGTAGCCACATCGTAACGCCGGATAAGTGCGAAACCGCCTACGAGACGGTCAAATGTTATTTCAATGCGCACGACGAGGTGATCAAATTCTGCCACCTACTAGTGCTGGAGTGA"
    #s1_AGAP010815_RA = 'ATGTGGCAGTTCATAAGGTCACGAATATTAACGGTGATAATCTTCATAGGTGCTGCTCATGGGCTACTGGTTGTGGGTCCGAAATTTATACGGGCCAACCAGGAATACACTCTGGTGATCAGCAACTTTAACTCACAGCTAAGCAAAGTGGACCTGCTGTTAAAACTGGAAGGCGAAACTGATAATGGTTTAAGCGTTCTGAACGTTACCAAGATGGTTGACGTGCGACGTAATATGAACCGAATGATCAACTTCAATATGCCTGAGGATCTGACGGCTGGAAACTACAAAATAACTATCGATGGACAGCGTGGCTTCAGCTTTCACAAGGAGGCAGAGCTGGTGTATCTCAGCAAATCGATATCGGGGCTAATACAGGTCGATAAGCCCGTATTTAAACCTGGGGATACGGTGAACTTCCGTGTGATCGTGCTGGACACGGAGCTGAAACCGCCGGCGAGGGTCAAGTCGGTTTATGTAACTATACGAGATCCTCAGCGCAATGTGATTCGCAAATGGTCCACGGCAAAACTGTATGCCGGTGTGTTCGAGAGCGATCTACAGATAGCGCCTACTCCAATGCTCGGGGTCTGGAATATCTCGGTGGAGGTGGAAGGAGAAGAGCTTGTGTCAAAGACGTTTGAGGTGAAGGAGTACGTGTTGTCAACGTTCGACGTGCAGGTCATGCCATCGGTGATTCCACTGGAAGAGCATCAAGCTGTGAATCTTACAATCGAAGCGAACTATCACTTTGGTAAGCCAGTGCAAGGAGTGGCCAAGGTGGAGCTGTACCTAGACGACGATAAGCTAAAACTGAAAAAAGAGCTGACTGTGTACGGAAAGGGCCAGGTAGAGTTGCGCTTTGACAATTTTGCAATGGATGCGGATCAGCAGGATGTACCAGTGAAGGTGTCGTTCGTCGAGCAGTACACAAATCGTACGGTGGTCAAACAGTCACAAATCACGGTATATAGGTATGCGTACCGAGTAGAGTTGATAAAAGAGAGTCCACAGTTTCGTCCGGGACTCCCGTTCAAATGTGCGCTTCAGTTTACACACCATGATGGAACACCGGCTAAAGGCATTAGCGGTAAGGTAGAGGTATCCGATGTACGATTCGAAACGACAACAACGAGTGATAACGATGGATTGATTAAGCTCGAGCTGCAACCAAGTGAGGGTACTGAACAACTCAGTATTCACTTCAATGCTGTTGATGGATTCTTTTTTTATGAAGATGTGAATAAGGTAGAAACGGTTACGGATGCGTATATTAAACTGGAGCTGAAATCACCGCATCAAACGGAACAAATTGATGCGTTTCATGGTGACGTGCACGGAGCGCATGACATTCTTCGTGTACTATGTCATGTCAAAGGGCAATATCATCGATGCAGGATTCATGCGACCCAACAAGCAACCGAAGTACCTGTTGCAGCTGAACGCAACAGAAAAGATGATTCCGAGGGCGAAAATTCTCATCGCTACCGTAGCGGGCCGCACGGTGGTGTACGACTTCGCAGACCTCGATTTCCAAGAGCTTCGCAATAATTTTGATTTAAGCATTGACGAGCAAGAGATCAAGCCGGGACGACAAATCGAGCTGAGCATGTCTGGACGCCCAGGAGCGTACGTTGGGCTGGCCGCGTATGACAAAGCCTTGCTGCTTTTCAACAAGAACCACGACCTGTTCTGGGAGGACATTGGGCAGGTGTTTGATGGGTTCCATGCAATCAATGAGAACGAGTTTGACATATTCCACAGCTTGGGTCTGTTCGCCAGGACATTGGACGATATCTTGTTCGACAGTGCAAATGAAAAGACGGGGCGTAATGCACTGCAGTCAGGCAAGCCGATCGGCAAGCTGGTGTCGTATCGGACGAACTTCCAGGAATCGTGGTTGTGGAAAAATGTTTCCATCGGACGATCGGGAAGTCGCAAGTTGATCGAGGTAGTACCGGACACGACCACCTCCTGGTATCTGACGGGCTTCTCGATCGATCCCGTGTACGGGTTGGGTATCATCAAGAAGCCAATCCAGTTCACAACAGTCCAGCCGTTCTACATCGTAGAGAACTTACCATATTCAATCAAACGAGGCGAAGCGGTTGTGTTGCAGTTTACGCTGTTCAACAACCTTGGAGCGGAGTATATAGCCGATGTGACGCTGTACAATGTGGCCAACCAGACCGAGTTCGTCGGACGTCCAAATACGGATCTCAGCTACACCAAATCCGTGAGCGTTCCTCCAAAAGTTGGTGTGCCAATCTCGTTCCTCATCAAGGCCCGCAAGCTCGGCGAGATGGCGGTTCGTGTAAAGGCTTCGATAATGCTGGGACACGAAACGGACGCCCTGGAAAAGGTAATACGGGTGATGCCTGAAAGTTTGGTGCAGCCGAGAATGGATACACGCTTTTTCTGCTTCGACGATCACAAAAATCAAACGTTTCCGATCAACTTGGACATCAACAAGAAGGCCGACAGTGGATCGACAAAGATTGAGTTTCGACTAAATCCCAATTTGTTGACCACGGTCATCAAGAACCTGGACCATCTTCTCGGCGTTCCGACGGGATGTGGTGAGCAGAATATGGTCAAATTTGTTCCCAACATTTTGGTACTGGATTATTTGCATGCCATCGGGTCGAAAGAACAGCATCTAATCGACAAAGCTACGAATTTGTTGCGTCAAGGATATCAAAACCAGATGCGCTACCGTCAGACGGATGGTTCATTTGGTTTGTGGGAGACTACTAATGGTAGCGTGTTTCTCACCGCGTTCGTTGGCACATCGATGCAAACTGCAGTAAAATACATAAGCGATATTGATGCAGCAATGGTGGAGAAGGCATTGGATTGGTTAGCCTCGAAGCAGCATTTCTCGGGACGGTTTGACAAGGCCGGTGCAGAGTATCACAAAGAAATGCAAGGAGGGTTGCGCAATGGTGTGGCCCTCACATCATATGTGTTGATGGCATTGCTGGAGAATGACATTGCCAAAGCAAAGCACGCAGAGGTGATTCAAAAAGGAATGACCTATCTGAGCAATCAGTTTGGATCCATCAACAATGCATACGACCTATCGATAGCAACCTACGCGATGATGTTGAACGGACACACCATGAAGGAGGAGGCACTCAATAAGCTGATTGATATGTCTTTCATTGATGCTGATAAAAACGAACGGTTCTGGAACACAACGAATCCAATAGAAACCACCGCATATGCTCTGCTGTCGTTTGTGATGGCCGAGAAGTACACAGACGGTATACCGGTCATGAATTGGTTGGTGAATCAACGTTACGTTACCGGTAGCTTTCCGAGCACGCAAGACACGTTTGTGGGGCTGAAAGCGCTGACCAAAATGGCGGAAAAGATATCTCCGTCCCGAAACGACTACACCGTTCAACTGAAGTACAAGAAGAGTGCAAAATACTTCAAAATAAACTCGGAGCAAATTGATGTGGAAAACTTCGTGGATATACCGGAGGACACAAAAAAGCTCGAGATCAATGTGGGGGGCATTGGATTTGGGTTGTTAGAGGTGGTTTATCAATTTAATTTGAATCTCGTCAACTTTGAGAATAGATTCCAACTAGACCTGGAGAAACAGAACACAGGCTCTGACTACGAGCTGAGGCTGAAGGTCTGTGCCAGCTACATACCCCAGCTGACCGACAGACGATCGAACATGGCACTGATTGAGGTAACCTTACCGAGCGGTTACGTGGTTGATCGCAATCCGATCAGCGAGCAGACGAAGGTGAATCCGATTCAGAAAACTGAAATCCGTTACGGTGGCACTTCAGTCGTTTTATACTACGACAATATGGGCAGCGAGCGTAACTGTTTCACCCTGACCGCGTACAGACGCTTTAAGGTCGCATTGAAGCGTCCAGCGTATGTGGTTGTGTATGATTATTATAATACAAATCTGAACGCCATCAAAGTGTACGAAGTGGACAAGCAGAATTTGTGCGAAATCTGTGACGAAGAAGACTGTCCTGCAGAGTGCAAAAAATAG'
    #qry_seq_raw = s1_AGAP010815_RA
    #
    # AAEL001802 (aegypti)
    #     https://www.vectorbase.org/Aedes_aegypti/Gene/Sequence?db=core;g=AAEL001802;r=supercont1.43:685886-717122;t=AAEL001802-RA
    ref_seq_raw = sys.argv[2]
    # OBP7 aegypti: v
    #ref_seq_raw = "ATGATGGAACAGCTTATGCTGGCAGTTTTGCTGGCGGTTTTTCTCGGGCTCGTAGCAGATGTTACGATGGCCGCTCAAATCAAGGACAATTTGGAGCTACCCGAATATTACAAACGTCCGGCCAAAATTCTGCACAACATCTGTCTGGCAGAATCCGGTGCCATGGAGAGCAAACTAAAGCAGTGCATGGACGGAGTGCTTCATGACGACCGGGAAGTCAAGTGCTACATCCATTGTCTATTCGACAAGGTGGACGTAATCGACGAAGCAACCGGGCAGATCCTGTTGGACCGATTGGCACCACTGGCACCGGACAACGATGTGAAGGATGTGTTCAATCATTTGACCAAAGAGTGTGGTCATATCAAACTACAAGATTCCTGCGATACGGCGTACGAAGTGGCCAAATGTTACTTCGCGGCACACGATCAGGTCGTCAAATTCTGTCACCTGTTGATGGCTGATGTTACCAGCTAG"
    #s2_AAEL001802_RA = 'ATGTCGGTATTCATACAAACGGACAAACCGGTGTATACCCCGGGAGATCTGATACGTTTTCGGGTAATCGTGGTGGATGCTGACACTAGACCTGTGACTAGTATTAAAACGGTAAATATAGCGATCGACGATTCTGCAAAAAATTCCATTCGAAAGTGGCCTTATGCCAAGTTGTTAAACGGCATCTTTGAGTCACAAGTGCAATTAGCTTCTTCGCCTGTTCTTGGCACCTGGATTATCAACGTAACAGCTTCCGACGACATCATTGTCACCAAACAGATAGAAGTTAAGGAATATGTGTTGCCAAAATTTTTCGTGAAAGTTTACCCTTCGGAGGTTCTATTGGGGAAAAATAAGAAGGTTTCTCTTACCTTAGATGCCTATTACACGTTCAAAGAACCCGTCGACGGCAATTACAAAGTTGAGTTATTTTTGGACCATACCAAGAGAAAGCCTGACTTCATAAAAAGTGATCGAATCACCGGTAAAACATCACTTGAGTTTCAATTGAAAAATGAAGTAGACATTGATGGCGACGAGCAGTACACTGATGTCACGGTTGAAGTTGAAGTTGTCGAGGCATTTTCTAATCGCACAGTTAGTATAACTGAGAATATTCCGATTTATCGTCAGCCTTATACCGTGACCCTTCTTCCATCTGCACCATCATTTCGACCAGGAGTTCCATTCAATGTACAAATAGTTGTGAAAGATCAGCTTGGACACCCTCCTGCCGAAGAAAAAGCGGCATCAATTGACCTTACTGTAGAGTTCCATTTGCCCATTGACAGTGACACCAAATCTATCACTGTAGATCTGGACGAGAAAGGAACAGGTCAGCTCACATTAGAGCCCCGCCCAGACGCCCAAGAACTGAAAGTGAACGCTACATATGACTCTCAACAATACGATGTAATTCACGATCCGATACATGGTTTCAGTTCGCAAAGTAAGCAGTACATCACAGTAACTCTGAATCCAAAATACTATAACAACATTAAAGTCGATAAGGACATCGTACTGGACATCTCCTGCACTGAAACAATGACGCACTTCTCGTACATCGTTGTCACCAGAGGAAACATAGTGGAAGCATCGAACGTTCCTGTCAGGATAAAAAAGAAACATTCTCTGAGATTGAAAATGACTTCAAAAATGTCTCCGGAGTCGAGGCTTCTAGTGTACTATACAAACAGGGAGTATCTCATCTTTGATGATATTGAGCTGAAGTTCGATTCGTTCAACAACGACTTCAAATTCGATTTGAACGATGATGAGTATTTTCCAGGGCAATCAGTTTATATCGATGTATACGCTTCAAAGGATTCATACGTTGCGTTCAGTGGAATCGATGAAAGTGTACTCCTGGTAGGCAAAGAGCGCCATGACTTCAACAAAGGAGATGTGCTCAAGGAACTCGCTCTTTACGGAGCAACAAATGATGCCGAGTTTGACTTGTTCCACGTAAGTTTCATGTCAAATGGTTTGATTATTCCAGTTAATGTATCTGTAACTCGCTCACAGAATGCACGATTTGGTACTCTACTAGGAAGGACTAGGCAGCAAGCGATTGAAATTCGAACTCAATTCCTAGAATCCTGGTTATGGAAATCCTTTTCCATGGATGGTCGAAACAACTTCAAAGCAATAGAAGACTCGGTTCCGGATACTATTACAACGTATCACGTGTCAGGATTTGCTTTAAGTCCAACACTAGGTCTTGGAGTAATCCAACAACCAGTGAGTTTCACCGTTCGTAAAAAATTCTACTTGGTTGCAAATTTGCCTTACTCGATCAAACGGGGTGAAGTGGCGTTGATTCAGGTTACCGTCTTCAACTTCCTAGGAAGCAGCATAACAACCGATGTGACGCTGTTCAATAAACGCGATGAAATTGAGTTTGTCGAGAATGCATCCACTAATAATACACATCGAACAAAGGCGGTAATTGTCCCGAATAACAATGGAAAATCTGTATCATTTATGGTGAAAGCAAAGAAATTAGGACAGATTGCGATCAAATTCCAGGCGGTAAACCTGCTGGAAACGGATGCATTGGAGCACATGTTACGAGTAACCCCAGAGAGCCATCGCTATGAGAAAAATGTAGCTCGATTCGTTGAGCTACCAAAGTTTGAGACGCAAACTTTCGATGTGAAGCTGGACATTCCCAAAAATATCGACGAGGGTTCTGCTCAAATCAAATTCACGTTAGACCCGGACATTTTGGGAACAGCCATCAGCAACCTAGACGGGTTGATCCGGAAACCCTTTGGATGTGGCGAACAAAATATGCTCCATTTTGTGCCAAATATAGTCGTTTTGGATTATCTTAACGAAACCAACACAGCGGCAGAAGATGTGAGGACCAAAGCGATAAATTTTCTTAGCAGCGGATATCAAAACCAGCTACGCTACAAACGTTCGGATGGGGCCTTCAGTGTCTGGGGACAATCGTATGCTGGCAGTACATTTTTGACGGCCTTTGTGGCGAAATCATTCAAAATAGCAGCCAAATACATTCAGGTGGATAAGTCTATAGTAGACGCGGCATTCGACTGGTTAGTGAAACAACAACAATCAGATGGGCGGTTCCCAGAAGTGGGGCAAGTATTCCAAGCAGATATGCAGGGTGGGCTTCGTAATAACGGTTTTGCGCTTACCGCGTATGTTCTGATCGCTTTTGCTGAAAATAAGGAAGTATACAGAAAATACCAATCACAACTGAACAAAACTACTAACTTCATAGCAGATAGACTTGCTAATATGGAGAATCCATACGACCTCTCGCTGTCCACTTATGCGTTGATGCTAACAAATCATGGCAAGCGCACCGAGTTTCTTCACAAATTAGTCGAAAAGTCGATATTTGACCGCAATCAAACTGAGAGATATTGGGACAGCAAACCAGTTGATATTGAAGTTGCTGGATATGCTCTATTGTCATACGTAGCTGCCGGTAAATTATTGGATGCAACGCCTATCATGCGGTGGCTCAACAAGCAGCGTTATGGTCTCGGAGGCTATCCTGGAACTCAGGAAACATTCGTTGGATTGAAAGCATTGGCAACGTTCGCTGCAAATGTAACTAGTAGGAGAAACGAATATACTGTAAGGATATTCTACGAACCAAATGGTCGACGAACATTCGACGTACACATGCACAATTCGTTTAATATTCAAGAGCTTGACATTCCTAATAACATCAGAAAAATGAAGGTGGAAGTTGAAGGCATCGGCAGAGGCTTCTTCCAAGTGGCATATCAGTACTATCAAAATATGCAGGTGGCTAAGCCCAGTTTCAGCATTACAATTAATCAGCTTAACACCACGACGGAACACATGCAGCAATTGGACGTGTGTGTGAAATACATACCAAAAGAGGCTTATCAAAAATCGAATATGGCTTTGGTGGAAATATTCTTGCCTAGTGGGCTTGTAGCAGACTCAGATGCCATTACGGACAAGACTGGAGGAATTCGAAGAATTGAAAGACGTTTTTCGGACACCTCAGTAGTTATATATTATGATAATTTGGACCCCGAAGACAAGTGCTTCCGAGTGACTGCTTATCGTCGGTATAAAATTGCATTGCATTTGCCATCATATATTATAGTTTATGATTATTATAATTTTGAGCGCTTTGCCATTCAAAAGTACGAAGGAAAGGTGCTGCAGCTCTGCGATATTTGTGAAGACGAGGACTGCGAAACTTTATCATGTCAAAATAGCTCGAAATTGGCAATAATGTAA'
    #ref_seq_raw = s2_AAEL001802_RA

    ###########################
    # Benchmarking vs. MATLAB #
    ###########################

    # @done:@testing:benchmark vs. matlab: @todo: place in README these details
    #   using s1_AGAP010815_RA and s2_AAEL001802_RA (after aligning/trimming)
    # MATLAB whole dnds vs. andy-wenping whole dnds (msCorrect='approximate'): dN/dS: 
    #   1.0958 vs. 1.04939841193 (Elapsed time: 1.952491045, Total warnings: 0), respectively
    # MATLAB mean(sliding dnds) vs. andy-wenping mean(sliding dnds): Mean dN/dS over all 50-codon long windows: 
    #   1.7254581782566112 (smoothing?) vs. 1.49715496621 (incl. smoothing, stepLength=1) (Elapsed time: 2.29854488373, Total warnings: 173) (when taking the mean of all 50 window intervals in MATLAB's sliding dnds()
    # MATLAB sliding dnds plot vs. andy-wenping sliding dnds plot: similar sliding dnds plots: 
    #   see:  D:\Dropbox\00_HPC-LEAP\WP5 Thematic Cyprus Computational Biology\giannis_database

    # @benchmark:matlab's dnds()
    # dN/dS=0.15270463083614955 with msCorrect='approximate', dN/dS=0.16529268957638238 with 
    # msCorrect='exact', MATLAB gives: dN/dS=0.0221=dnds(seq1,seq2, 'geneticCode', 1, 'Method', 'NG'). This is a huge error, an order of magnitude off. 
    # seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
    # seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC' 
    # dnds('ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC', 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC')

    #
    # ALIGN REF vs. QUERY seq, then TRIM GAPS: 
    #     Next we want to remove all characters (trim) in each of the sequences that is aligned to a gap position**** in the top alignment
    #

    #  Raw input:
    #     ATGTGCTAA      qry
    #     ATGATTTGA      ref
    #
    # Alignment:
    #     ATGTGC----TAA  qry_seq
    #     ATG--A--TTTGA  ref_seq
    #        ^^ ^^^^     gap positions****
    #  Trimming:
    #     ATGCTAA        qry_seq (after trimming)
    #     ATGATGA        ref_seq (after trimming)

    # alignment scores are calculated from these two parameters, best scoring alignments have few gaps and gaps are small
    aln_gap_open = -10      
    aln_gap_extend = -0.5 # @todo: make cmd args later?

    qry_seq_aln, ref_seq_aln         = align_then_trim.align_query_vs_reference(qry_seq_raw, ref_seq_raw, aln_gap_open, aln_gap_extend)
    qry_seq_trimmed, ref_seq_trimmed = align_then_trim.trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )

    # @Note: benhamrking: it is qry_seq_trimmed and ref_seq_trimmed that should be used as input to MATLAB's dnds() function, for benchmarking

    #}} alternative 2

    # ///

    # alternative 1 {{

    # # @NOTE:uncomment below to achieve dnds of 0.15.. or 0.164 if using exact method, interestingly 
    # dnds_whole, warning_count = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=False)
    # print "dN/dS: "+str(dnds_whole)

    # }} 1 alternative 2 {{

    # @DONE: work with the sliding window version instead of dnds_whole
    dnds_slide_dict, warning_count = dnds( qry_seq_trimmed, ref_seq_trimmed, potential_changes, observed_changes, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )
    #
    # Plot the sliding window values
    #
    dnds_sliding_vec, dnds_sliding_mean = plot_dnds_sliding(dnds_slide_dict)
    
    print "Mean dN/dS over all windows: "+str(dnds_sliding_mean)
    # @todo: show the user also the whole dnds value somewhere in the webpage

    # }} alternative 2  
