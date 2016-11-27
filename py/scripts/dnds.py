############
# IMPORTS: #
############

import pickle
import math
import pdb

# @TODO:REMOVE: \/ and \/\/
# with open('../data/observed_changes_dict.p','rb') as f_observed:
#     changes_observed = pickle.load(f_observed)


# @TODO:REMOVE: \/ and \/\/: the final script cannot use absolut paths, note: the path from which this script is executed is where current working directory is, in this case it should be a .html that is in {root} calling this script residing in {root}/scripts directory
 # changes_observed = pickle.load(open('/home/qiime/Desktop/hpcleap_wp6_compbio/hpcleap_bioinf/data/observed_changes_dict.p','rb'))
# changes_potential= pickle.load(open('/home/qiime/Desktop/hpcleap_wp6_compbio/hpcleap_bioinf/data/potential_changes_dict.p','rb'))

with open('./py/data/observed_changes_dict.p','rb') as f:
    changes_observed = pickle.load(f)

with open('./py/data/potential_changes_dict.p','rb') as f:
    changes_potential= pickle.load(f)

#############
# FUNCTIONS #
#############



def chunks(l, n):
    """ Yield successive n-sized chunks from l. """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def dnds( seq1, seq2, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):
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

        msCorrect, a string to toggle between multiple-substitution correction methods:
            "approximate", "exact" (@TODO: make sure it actually is exact... it could be
             something else)
            
            e.g. msCorrect = 'approximate'

        sliding, a boolean to toggle between sliding window analysis (vector of dN/dS values at successive chunks of sequence) or whole sequence analysis (a single 
            dN/dS value for the given pair of input sequences), either: True, False
            
            e.g. sliding = False

        windowLength, an integer specifying the width of the sliding window, measured in DNA basepairs, from 1-to-length(seq1)
            e.g. 

    """

    # TODO: stop codons to deal with, reject
    # TODO: ambiguous bases to deal with: gaps, Ns, Xs

    # STATS per CODON-PAIR:
    codons_seq1   = [codon for codon in chunks(seq1,3)]  #splits
    codons_seq2   = [codon for codon in chunks(seq2,3)]  #splits
    codons_paired = [pair for pair in zip(codons_seq1,codons_seq2) if (len(pair[0])+len(pair[1]))==6] # aligned codons are paired into tuples, excess codons are truncated

    changes = {'observed':{'S':[],'N':[]},'potential':{'S':[],'N':[]}}

    for pair in codons_paired:
        changes['potential']['S'].append(changes_potential[pair]['S'])
        changes['potential']['N'].append(changes_potential[pair]['N'])
        changes['observed']['S'].append(changes_observed[pair]['S'])
        changes['observed']['N'].append(changes_observed[pair]['N'])

    list_S  = changes['potential']['S']
    list_Sd = changes['observed']['S']

    list_N  = changes['potential']['N']
    list_Nd = changes['observed']['N']

    #pdb.set_trace()

    if sliding:
        # STATS per WINDOW
        intervals    = range(0,len(codons_paired)-windowLength+1,stepLength)
        windows      = zip(intervals,[i + windowLength - 1 for i in intervals]) 
        window_stats = {}

        for window in windows:
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
            # technically it woud be best so exclude these from downstream analyses? 
            # e.g. missing value/datapoint on a plot of dN/dS (y-axis) vs. window interval (x-axis)
            except ValueError:
                window_stats[window]['dNdS'] = None
            except ZeroDivisionError:
                window_stats[window]['dNdS'] = float('Inf')
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
                raise ValueError("Approximate multiple-substitutions correction cannot be achieved, pS>=3/4, try alternative value for argument: msCorrect...") 

            if (pN>=3./4.):
                raise ValueError("Approximate multiple-substitutions correction cannot be achieved, pN>=3/4, try alternative value for argument: msCorrect...") 

            # @TODO: if the proportion of syn is >3/4 it might break
            try:
                dS  = -(3./4.)*math.log(1.-((4./3.)*pS))
                dN  = -(3./4.)*math.log(1.-((4./3.)*pN))
            except ValueError:
                pdb.set_trace()

        else: # @TODO: is this the exact one? Or something else? 
            #print "holy cow (whole)"
            dS  = pS  # i.e. dS = Sd/S
            dN  = pN

        print "dN/dS: "+str(dN/dS)

        return (dN/dS)  # i.e. omega = dN/dS = (Nd/N)/(Sd/S)


# seq1 = 'CTTTTTAACGAAAAAGAGGCAGATGA'
# seq2 = 'ATGGCCCACTTCCCAGGGTTTGGACA'


# dN/dS=0.15270463083614955 with msCorrect='approximate', dN/dS=0.16529268957638238 with 
# msCorrect='exact', MATLAB gives: dN/dS=0.0221=dnds(seq1,seq2, 'geneticCode', 1, 'Method', 'NG'). This is a huge error, an order of magnitude off. 
seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC' 

# seq1 = 'CGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
# seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAA' 

# dnds('ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC', 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC')

# ////


# seq1 = 'ATGGAGCTCAGACTAGAGCGGGGCTGTTGACGTTTGGAGTTGAAAAAATCTATTATACCAATCGGCTTCAACGTGCTCCACGGCAGGCGCCTGACGAGGGGCCCACACCGAGGAAGTAGACTGTTGCACGTTGGGGATAGCGGTAGCTAACTAAGACGCCTGCCACAACAGCAGTATCAAACCCGTACAAAGGGAACATCCACACTTTGGTGAATCGAAGCGCGGCATCAGAATTTCCTTTTGGATACCTGATACAAAGCCCATCGTGGTCCTTAGACTTCGTACACTTACACCTGCACCGCG'
# seq2 = 'ATGCGCATGTGGAATTAGAGGCGAAGTACGATCCCTAGACCGACGTACGATGCAACTGTGTGGATGTGACGAGCTTCTTTTATATGCTTCGCCCGCCGGACCGGCCTCGGCATGGCGTAGCAGTGCACAAGCAAATGACAATTAACCACCGTGTATTCGTTATAACATCAGGCAGTTTAAGTCGGGACAATAGGAGCCGCAATACACAGTTTACCGCATCTTGACCTAACTGACATACTGCCATGGACGACTAGCCATGCCACTGGCTCTTAGATAGCCCGATACAGTGATTATGAAAGGTTT' 







# Note: msCorrect="aproximate"
#dnds_sliding = dnds( seq1,seq2, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )

# # Note: msCorrext="exact"
# dnds_sliding = dnds( seq1,seq2, msCorrect='exact', sliding=True, windowLength=50, stepLength=1 )

dnds_whole = dnds( seq1,seq2, msCorrect='approximate', sliding=False )

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

