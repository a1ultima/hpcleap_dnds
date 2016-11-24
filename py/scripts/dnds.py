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

with open('./data/observed_changes_dict.p','rb') as f:
    changes_observed = pickle.load(f)

with open('./data/potential_changes_dict.p','rb') as f:
    changes_potential= pickle.load(f)

#############
# FUNCTIONS #
#############

pdb.set_trace()

def chunks(l, n):
    """ Yield successive n-sized chunks from l. """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def dnds( seq1, seq2, msCorrect='approximate', sliding=False, windowLength=3, stepLength=1):

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
                    dN = -(3./4.)*math.log(1.-(4./3.)*pN)
                    dS = -(3./4.)*math.log(1.-(4./3.)*pS)

                # elif msCorrect=='exact':
                #     d=ln(1-p*4/3)/ln(1-3/(4*N))

                else:
                    dN = pN
                    dS = pS
                window_stats[window]['dNdS'] = dN/dS
            except ValueError:
                window_stats[window]['dNdS'] = None
            except ZeroDivisionError:
                window_stats[window]['dNdS'] = float('Inf')
        return window_stats
    else:
        S   = sum(list_S)
        Sd  = sum(list_Sd)
        pS  = Sd/S 
        N   = sum(list_N)
        Nd  = sum(list_Nd)
        pN  = Nd/N

        if msCorrect:
            dS  = -(3./4.)*math.log(1.-(4./3.)*pS)
            dN  = -(3./4.)*math.log(1.-(4./3.)*pN)
        else: 
            dS  = pS
            dN  = pN

        return dN/dS


# seq1 = 'CTTTTTAACGAAAAAGAGGCAGATGA'
# seq2 = 'ATGGCCCACTTCCCAGGGTTTGGACA'

seq1 = 'ATGCGCAAATACTCCCCCTTCCGAAATGGATACATGGAACCCACCCTTGGGCAGCACCTCCCAACCCTGTCTTTTCCAGACCCCGGACTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGGCTCCGTTGTCTGCATGTACCTCTACCAGCTTTCCCCCCCCATCACCTGGCCCCTCCTGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAACGAATAGAAAAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTGCCCACCACCCTTTTCCAGCCTGCTAGGGCACCCGTCACGCTGACAGCCTGGCAAAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGATTTCCGGGCCCTGCCCTAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCCTTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCCTCATTTCTACTCTCACACGGCCTCATACAGTACTCTTCCTTTCATAATTTGCATCTCCTATTTGAAGAATACACCAACATCCCCATTTCTCTACTTTTTAACGAAAAAGAGGCAGATGACAATGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCTCAGTGAAAAACATTTCCGTGAAACAGAAGTC'
seq2 = 'ATGCGCAAGTACTCCCCCTTCCGAAACGGATACATGGAACCCACCCTTGGGCAACACCTCCCAACCCTGTCTTTTCCAGACCCCGGCCTCCGGCCCCAAAACCTGTACACCCTCTGGGGAGACTCTGTTGTCTGCCTGTACCTCTACCAGCTCTCCCCCCCCATCACCTGGCCCCTCCCGCCCCATGTGATTTTTTGCCACCCCGGCCAGCTCGGGGCCTTCCTCACCAATGTTCCCTACAAGCGTATGGAAGAACTCCTCTATAAAATTTCCCTTACCACAGGGGCCCTAATAATTCTACCCGAGGACTGTTTACCAACCACCCTTTTCCAGCCTGCTAGGGCCCCCGTCACGTTGACCGCCTGGCAGAACGGCCTCCTTCCGTTCCACTCAACCCTCACCACTCCAGGCCTTATTTGGACATTTACCGATGGCACGCCTATGGTTTCCGGACCCTGCCCCAAAGATGGCCAGCCATCTTTAGTACTACAGTCCTCCTCATTTATATTTCACAAATTTCAAACCAAGGCCTACCACCCTTCATTTCTACTCTCACACGGCCTCATACAGTACTCCTCCTTTCACAATTTACATCTCCTTTTTGAAGAATACACCAACATCCCCGTTTCTCTACTTTTTAACGAAAAAGAGGCAAATGACACTGACCATGAGCCCCAAATATCCCCCGGGGGCTTAGAGCCTCCCGCTGAAAAACATTTCCGCGAAACAGAAGTC'

# Note: msCorrect="aproximate"
dnds_sliding = dnds( seq1,seq2, msCorrect='approximate', sliding=True, windowLength=50, stepLength=1 )

# # Note: msCorrext="exact"
# dnds_sliding = dnds( seq1,seq2, msCorrect='exact', sliding=True, windowLength=50, stepLength=1 )

dnds_whole = dnds( seq1,seq2, msCorrect=True, sliding=False )

#dnds_sliding

import pprint
pprint.pprint(dnds_sliding)
pprint.pprint(dnds_whole)

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

