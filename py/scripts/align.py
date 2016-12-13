
""" Performs Needleman-Wunsch global alignment using DNA substitution matrix, 
on a pair of DNA sequences(see below)**, one representing the "query" sequence, 
and the other a "reference" sequence. Followed by a gap trimming process,
as outlined below*** to format each sequence prior to running the dN/dS 
calculations later on (e.g. dnds.py). 

   See Bio.pairwise2.alignxs for more details: 
   http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html

** Input sequences: both sequences represent a CDS-transcript of 
interest (exons-only). To give an example, say the "query" sequence 
corresponds to the Anopheles gambiae TEP1 protein's transcript id 
(VectorBase: AGAP010815), and the other "reference" sequence corresponds 
to the orthologous Aedes aegypti TEP1 protein's trancript id (VectorBase: 
AAEL001802). The orthologous reference sequence is required, along with 
the query sequence of interest, in order to estimate dN/dS values (sliding
 or whole), which requires simulations/approximations of ancestral DNA 
substitution mutations (synonymous, non-synonymous changes).  

***Pipeline outline:

   1. Once the two orthologous transcript IDs are supplied by the user, 
      use REST queries, to VectorBase, to retrieve the DNA sequence data, 
      one FASTA format DNA sequence per user-input CDS transcript id
      (one for the "query" and one for the "reference").

   2. Globally align the two FASTA sequences (above) using Needleman-Wunsch 
      algorithm. 

   3. Trim out all basepairs (DNA letters), from the pair of aligned 
      sequences, corresponding to gap positions****. Mutations leading 
      to such gaps are not considered in dN/dS calculation, and so have 
      undefined values, hence trimming is necessary for downstream 
      analyses (e.g. dnds.py). 

****Search: "gap (genetics)" online, for more information.

"""

###########
# IMPORTS #
###########

import numpy as np
import pdb

#############
# FUNCTIONS #
#############


def align_query_vs_reference( qry_seq, ref_seq, aln_gap_open = -10, aln_gap_extend = -0.5 ):

   """ Globally align two input sequences using NW-algorithm, return highest-scoring alignment.

   ARGS:
      
      qry_seq, ref_seq,    two different orthologous CDS (exon-only) DNA 
                           sequences to be aligned (see: "** Input sequences").

         e.g. qry_seq = s1_AGAP010815_RA = "ATGTGGCAGTTCATAAGGTCACGAATATTAACGGTGATAATCTTCATAGGTGCTGCTCATGGGCTACTGGTTGTGGGTCCGAAATTTATACGGGCCAACCAGGAATACACTCTGGTGATCAGCAACTTTAACTCACAGCTAAGCAAAGTGGACCTGCTGTTAAAACTGGAAGGCGAAACTGATAATGGTTTAAGCGTTCTGAACGTTACCAAGATGGTTGACGTGCGACGTAATATGAACCGAATGATCAACTTCAATATGCCTGAGGATCTGACGGCTGGAAACTACAAAATAACTATCGATGGACAGCGTGGCTTCAGCTTTCACAAGGAGGCAGAGCTGGTGTATCTCAGCAAATCGATATCGGGGCTAATACAGGTCGATAAGCCCGTATTTAAACCTGGGGATACGGTGAACTTCCGTGTGATCGTGCTGGACACGGAGCTGAAACCGCCGGCGAGGGTCAAGTCGGTTTATGTAACTATACGAGATCCTCAGCGCAATGTGATTCGCAAATGGTCCACGGCAAAACTGTATGCCGGTGTGTTCGAGAGCGATCTACAGATAGCGCCTACTCCAATGCTCGGGGTCTGGAATATCTCGGTGGAGGTGGAAGGAGAAGAGCTTGTGTCAAAGACGTTTGAGGTGAAGGAGTACGTGTTGTCAACGTTCGACGTGCAGGTCATGCCATCGGTGATTCCACTGGAAGAGCATCAAGCTGTGAATCTTACAATCGAAGCGAACTATCACTTTGGTAAGCCAGTGCAAGGAGTGGCCAAGGTGGAGCTGTACCTAGACGACGATAAGCTAAAACTGAAAAAAGAGCTGACTGTGTACGGAAAGGGCCAGGTAGAGTTGCGCTTTGACAATTTTGCAATGGATGCGGATCAGCAGGATGTACCAGTGAAGGTGTCGTTCGTCGAGCAGTACACAAATCGTACGGTGGTCAAACAGTCACAAATCACGGTATATAGGTATGCGTACCGAGTAGAGTTGATAAAAGAGAGTCCACAGTTTCGTCCGGGACTCCCGTTCAAATGTGCGCTTCAGTTTACACACCATGATGGAACACCGGCTAAAGGCATTAGCGGTAAGGTAGAGGTATCCGATGTACGATTCGAAACGACAACAACGAGTGATAACGATGGATTGATTAAGCTCGAGCTGCAACCAAGTGAGGGTACTGAACAACTCAGTATTCACTTCAATGCTGTTGATGGATTCTTTTTTTATGAAGATGTGAATAAGGTAGAAACGGTTACGGATGCGTATATTAAACTGGAGCTGAAATCACCGCATCAAACGGAACAAATTGATGCGTTTCATGGTGACGTGCACGGAGCGCATGACATTCTTCGTGTACTATGTCATGTCAAAGGGCAATATCATCGATGCAGGATTCATGCGACCCAACAAGCAACCGAAGTACCTGTTGCAGCTGAACGCAACAGAAAAGATGATTCCGAGGGCGAAAATTCTCATCGCTACCGTAGCGGGCCGCACGGTGGTGTACGACTTCGCAGACCTCGATTTCCAAGAGCTTCGCAATAATTTTGATTTAAGCATTGACGAGCAAGAGATCAAGCCGGGACGACAAATCGAGCTGAGCATGTCTGGACGCCCAGGAGCGTACGTTGGGCTGGCCGCGTATGACAAAGCCTTGCTGCTTTTCAACAAGAACCACGACCTGTTCTGGGAGGACATTGGGCAGGTGTTTGATGGGTTCCATGCAATCAATGAGAACGAGTTTGACATATTCCACAGCTTGGGTCTGTTCGCCAGGACATTGGACGATATCTTGTTCGACAGTGCAAATGAAAAGACGGGGCGTAATGCACTGCAGTCAGGCAAGCCGATCGGCAAGCTGGTGTCGTATCGGACGAACTTCCAGGAATCGTGGTTGTGGAAAAATGTTTCCATCGGACGATCGGGAAGTCGCAAGTTGATCGAGGTAGTACCGGACACGACCACCTCCTGGTATCTGACGGGCTTCTCGATCGATCCCGTGTACGGGTTGGGTATCATCAAGAAGCCAATCCAGTTCACAACAGTCCAGCCGTTCTACATCGTAGAGAACTTACCATATTCAATCAAACGAGGCGAAGCGGTTGTGTTGCAGTTTACGCTGTTCAACAACCTTGGAGCGGAGTATATAGCCGATGTGACGCTGTACAATGTGGCCAACCAGACCGAGTTCGTCGGACGTCCAAATACGGATCTCAGCTACACCAAATCCGTGAGCGTTCCTCCAAAAGTTGGTGTGCCAATCTCGTTCCTCATCAAGGCCCGCAAGCTCGGCGAGATGGCGGTTCGTGTAAAGGCTTCGATAATGCTGGGACACGAAACGGACGCCCTGGAAAAGGTAATACGGGTGATGCCTGAAAGTTTGGTGCAGCCGAGAATGGATACACGCTTTTTCTGCTTCGACGATCACAAAAATCAAACGTTTCCGATCAACTTGGACATCAACAAGAAGGCCGACAGTGGATCGACAAAGATTGAGTTTCGACTAAATCCCAATTTGTTGACCACGGTCATCAAGAACCTGGACCATCTTCTCGGCGTTCCGACGGGATGTGGTGAGCAGAATATGGTCAAATTTGTTCCCAACATTTTGGTACTGGATTATTTGCATGCCATCGGGTCGAAAGAACAGCATCTAATCGACAAAGCTACGAATTTGTTGCGTCAAGGATATCAAAACCAGATGCGCTACCGTCAGACGGATGGTTCATTTGGTTTGTGGGAGACTACTAATGGTAGCGTGTTTCTCACCGCGTTCGTTGGCACATCGATGCAAACTGCAGTAAAATACATAAGCGATATTGATGCAGCAATGGTGGAGAAGGCATTGGATTGGTTAGCCTCGAAGCAGCATTTCTCGGGACGGTTTGACAAGGCCGGTGCAGAGTATCACAAAGAAATGCAAGGAGGGTTGCGCAATGGTGTGGCCCTCACATCATATGTGTTGATGGCATTGCTGGAGAATGACATTGCCAAAGCAAAGCACGCAGAGGTGATTCAAAAAGGAATGACCTATCTGAGCAATCAGTTTGGATCCATCAACAATGCATACGACCTATCGATAGCAACCTACGCGATGATGTTGAACGGACACACCATGAAGGAGGAGGCACTCAATAAGCTGATTGATATGTCTTTCATTGATGCTGATAAAAACGAACGGTTCTGGAACACAACGAATCCAATAGAAACCACCGCATATGCTCTGCTGTCGTTTGTGATGGCCGAGAAGTACACAGACGGTATACCGGTCATGAATTGGTTGGTGAATCAACGTTACGTTACCGGTAGCTTTCCGAGCACGCAAGACACGTTTGTGGGGCTGAAAGCGCTGACCAAAATGGCGGAAAAGATATCTCCGTCCCGAAACGACTACACCGTTCAACTGAAGTACAAGAAGAGTGCAAAATACTTCAAAATAAACTCGGAGCAAATTGATGTGGAAAACTTCGTGGATATACCGGAGGACACAAAAAAGCTCGAGATCAATGTGGGGGGCATTGGATTTGGGTTGTTAGAGGTGGTTTATCAATTTAATTTGAATCTCGTCAACTTTGAGAATAGATTCCAACTAGACCTGGAGAAACAGAACACAGGCTCTGACTACGAGCTGAGGCTGAAGGTCTGTGCCAGCTACATACCCCAGCTGACCGACAGACGATCGAACATGGCACTGATTGAGGTAACCTTACCGAGCGGTTACGTGGTTGATCGCAATCCGATCAGCGAGCAGACGAAGGTGAATCCGATTCAGAAAACTGAAATCCGTTACGGTGGCACTTCAGTCGTTTTATACTACGACAATATGGGCAGCGAGCGTAACTGTTTCACCCTGACCGCGTACAGACGCTTTAAGGTCGCATTGAAGCGTCCAGCGTATGTGGTTGTGTATGATTATTATAATACAAATCTGAACGCCATCAAAGTGTACGAAGTGGACAAGCAGAATTTGTGCGAAATCTGTGACGAAGAAGACTGTCCTGCAGAGTGCAAAAAATAG"
         e.g. ref_seq = s2_AAEL001802_RA = "ATGTCGGTATTCATACAAACGGACAAACCGGTGTATACCCCGGGAGATCTGATACGTTTTCGGGTAATCGTGGTGGATGCTGACACTAGACCTGTGACTAGTATTAAAACGGTAAATATAGCGATCGACGATTCTGCAAAAAATTCCATTCGAAAGTGGCCTTATGCCAAGTTGTTAAACGGCATCTTTGAGTCACAAGTGCAATTAGCTTCTTCGCCTGTTCTTGGCACCTGGATTATCAACGTAACAGCTTCCGACGACATCATTGTCACCAAACAGATAGAAGTTAAGGAATATGTGTTGCCAAAATTTTTCGTGAAAGTTTACCCTTCGGAGGTTCTATTGGGGAAAAATAAGAAGGTTTCTCTTACCTTAGATGCCTATTACACGTTCAAAGAACCCGTCGACGGCAATTACAAAGTTGAGTTATTTTTGGACCATACCAAGAGAAAGCCTGACTTCATAAAAAGTGATCGAATCACCGGTAAAACATCACTTGAGTTTCAATTGAAAAATGAAGTAGACATTGATGGCGACGAGCAGTACACTGATGTCACGGTTGAAGTTGAAGTTGTCGAGGCATTTTCTAATCGCACAGTTAGTATAACTGAGAATATTCCGATTTATCGTCAGCCTTATACCGTGACCCTTCTTCCATCTGCACCATCATTTCGACCAGGAGTTCCATTCAATGTACAAATAGTTGTGAAAGATCAGCTTGGACACCCTCCTGCCGAAGAAAAAGCGGCATCAATTGACCTTACTGTAGAGTTCCATTTGCCCATTGACAGTGACACCAAATCTATCACTGTAGATCTGGACGAGAAAGGAACAGGTCAGCTCACATTAGAGCCCCGCCCAGACGCCCAAGAACTGAAAGTGAACGCTACATATGACTCTCAACAATACGATGTAATTCACGATCCGATACATGGTTTCAGTTCGCAAAGTAAGCAGTACATCACAGTAACTCTGAATCCAAAATACTATAACAACATTAAAGTCGATAAGGACATCGTACTGGACATCTCCTGCACTGAAACAATGACGCACTTCTCGTACATCGTTGTCACCAGAGGAAACATAGTGGAAGCATCGAACGTTCCTGTCAGGATAAAAAAGAAACATTCTCTGAGATTGAAAATGACTTCAAAAATGTCTCCGGAGTCGAGGCTTCTAGTGTACTATACAAACAGGGAGTATCTCATCTTTGATGATATTGAGCTGAAGTTCGATTCGTTCAACAACGACTTCAAATTCGATTTGAACGATGATGAGTATTTTCCAGGGCAATCAGTTTATATCGATGTATACGCTTCAAAGGATTCATACGTTGCGTTCAGTGGAATCGATGAAAGTGTACTCCTGGTAGGCAAAGAGCGCCATGACTTCAACAAAGGAGATGTGCTCAAGGAACTCGCTCTTTACGGAGCAACAAATGATGCCGAGTTTGACTTGTTCCACGTAAGTTTCATGTCAAATGGTTTGATTATTCCAGTTAATGTATCTGTAACTCGCTCACAGAATGCACGATTTGGTACTCTACTAGGAAGGACTAGGCAGCAAGCGATTGAAATTCGAACTCAATTCCTAGAATCCTGGTTATGGAAATCCTTTTCCATGGATGGTCGAAACAACTTCAAAGCAATAGAAGACTCGGTTCCGGATACTATTACAACGTATCACGTGTCAGGATTTGCTTTAAGTCCAACACTAGGTCTTGGAGTAATCCAACAACCAGTGAGTTTCACCGTTCGTAAAAAATTCTACTTGGTTGCAAATTTGCCTTACTCGATCAAACGGGGTGAAGTGGCGTTGATTCAGGTTACCGTCTTCAACTTCCTAGGAAGCAGCATAACAACCGATGTGACGCTGTTCAATAAACGCGATGAAATTGAGTTTGTCGAGAATGCATCCACTAATAATACACATCGAACAAAGGCGGTAATTGTCCCGAATAACAATGGAAAATCTGTATCATTTATGGTGAAAGCAAAGAAATTAGGACAGATTGCGATCAAATTCCAGGCGGTAAACCTGCTGGAAACGGATGCATTGGAGCACATGTTACGAGTAACCCCAGAGAGCCATCGCTATGAGAAAAATGTAGCTCGATTCGTTGAGCTACCAAAGTTTGAGACGCAAACTTTCGATGTGAAGCTGGACATTCCCAAAAATATCGACGAGGGTTCTGCTCAAATCAAATTCACGTTAGACCCGGACATTTTGGGAACAGCCATCAGCAACCTAGACGGGTTGATCCGGAAACCCTTTGGATGTGGCGAACAAAATATGCTCCATTTTGTGCCAAATATAGTCGTTTTGGATTATCTTAACGAAACCAACACAGCGGCAGAAGATGTGAGGACCAAAGCGATAAATTTTCTTAGCAGCGGATATCAAAACCAGCTACGCTACAAACGTTCGGATGGGGCCTTCAGTGTCTGGGGACAATCGTATGCTGGCAGTACATTTTTGACGGCCTTTGTGGCGAAATCATTCAAAATAGCAGCCAAATACATTCAGGTGGATAAGTCTATAGTAGACGCGGCATTCGACTGGTTAGTGAAACAACAACAATCAGATGGGCGGTTCCCAGAAGTGGGGCAAGTATTCCAAGCAGATATGCAGGGTGGGCTTCGTAATAACGGTTTTGCGCTTACCGCGTATGTTCTGATCGCTTTTGCTGAAAATAAGGAAGTATACAGAAAATACCAATCACAACTGAACAAAACTACTAACTTCATAGCAGATAGACTTGCTAATATGGAGAATCCATACGACCTCTCGCTGTCCACTTATGCGTTGATGCTAACAAATCATGGCAAGCGCACCGAGTTTCTTCACAAATTAGTCGAAAAGTCGATATTTGACCGCAATCAAACTGAGAGATATTGGGACAGCAAACCAGTTGATATTGAAGTTGCTGGATATGCTCTATTGTCATACGTAGCTGCCGGTAAATTATTGGATGCAACGCCTATCATGCGGTGGCTCAACAAGCAGCGTTATGGTCTCGGAGGCTATCCTGGAACTCAGGAAACATTCGTTGGATTGAAAGCATTGGCAACGTTCGCTGCAAATGTAACTAGTAGGAGAAACGAATATACTGTAAGGATATTCTACGAACCAAATGGTCGACGAACATTCGACGTACACATGCACAATTCGTTTAATATTCAAGAGCTTGACATTCCTAATAACATCAGAAAAATGAAGGTGGAAGTTGAAGGCATCGGCAGAGGCTTCTTCCAAGTGGCATATCAGTACTATCAAAATATGCAGGTGGCTAAGCCCAGTTTCAGCATTACAATTAATCAGCTTAACACCACGACGGAACACATGCAGCAATTGGACGTGTGTGTGAAATACATACCAAAAGAGGCTTATCAAAAATCGAATATGGCTTTGGTGGAAATATTCTTGCCTAGTGGGCTTGTAGCAGACTCAGATGCCATTACGGACAAGACTGGAGGAATTCGAAGAATTGAAAGACGTTTTTCGGACACCTCAGTAGTTATATATTATGATAATTTGGACCCCGAAGACAAGTGCTTCCGAGTGACTGCTTATCGTCGGTATAAAATTGCATTGCATTTGCCATCATATATTATAGTTTATGATTATTATAATTTTGAGCGCTTTGCCATTCAAAAGTACGAAGGAAAGGTGCTGCAGCTCTGCGATATTTGTGAAGACGAGGACTGCGAAACTTTATCATGTCAAAATAGCTCGAAATTGGCAATAATGTAA"

      aln_gap_open, aln_gap_extend -   float/integer specifying how to score alignments, i.e. fewest gaps and smallest gaps get higher score
      
         e.g. aln_gap_open = -10      
         e.g. aln_gap_extend = -0.5


  
   RETURNS:

      @todo

   """

   from Bio import pairwise2

   #from Bio.SubsMat import MatrixInfo as matlist
   #matrix = matlist.blosum62

   #alns = pairwise2.align.globalds(p53_human, p53_mouse, matrix, gap_open, gap_extend)
   # 
   # @todo: Q: are you sure about: we chose alignxs, which does not penaliste substitutions and instead only penalises gap open and gap extensions, i.e. is alignxs just penalise gaps and not subs?
   alns = pairwise2.align.globalxs(qry_seq, ref_seq, aln_gap_open, aln_gap_extend)  # @todo: make sure open and extend are in right order

   top_aln = alns[0] # sorte list of alternative alignments, highest-scoring alignent is alns[0]

   # @todo:make sure the 0 is actually query seq and 1 is actually reference seq

   qry_seq_aligned = top_aln[0]
   ref_seq_aligned = top_aln[1]

   # #
   # # ASCII diagram of alignment (s1 on top, s2 on bot and differences in between)
   # #
   # from Bio.pairwise2 import format_alignment  
   # print(format_alignment(*top_aln))

   return qry_seq_aligned, ref_seq_aligned


def trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln ):

   """ Trim the gaps away from aligned sequences generated by align_query_vs_reference. 
   Aligned then gap-trimmed sequenes are thus in a format ready for input to dnds.py.

      e.g.  Aligned qry and ref seqs:
               ATGTGC----TAA  qry_seq
               ATG--A--TTTGA  ref_seq
                  ^^ ^^^^     gap positions****
            After trimming:
               ATGCTAA        qry_seq (after trimming)
               ATGATGA        ref_seq (after trimming)

   ARGS:

      qry_seq_aln, ref_seq_aln -    @todo

         e.g. qry_seq_aln = 'ATGTGC----TAA'
         e.g. ref_seq_aln = 'ATG--A--TTTGA'

   RETURNS:

      qry_seq_trimmed, ref_seq_trimmed -    @todo

         e.g. qry_seq_trimmed = 'ATGCTAA'
         e.g. ref_seq_trimmed = 'ATGATGA'

   """

   # alternative 1 {{  # progressiely delete gaps and re-normalise index to delete with

   # qry_gap_positions = [i for i, letter in enumerate(qry_seq_aln) if letter.upper() not in ['A','G','T','C'] ]

   # # temp variables start with a copy of a and b as lists 
   # #  (to allow index-based deletions of elements)
   # qry_after_trim_qry_gaps = list(qry_seq_aln)
   # ref_after_trim_qry_gaps = list(ref_seq_aln)

   # re_normalise_i = 0  # since the sequences get progressively shorter with each element deletion, we must correct the index used to delete with, for every deletion

   # # corresponding to gaps in a: delete letters in a and b 
   # for i in qry_gap_positions:
   #    del qry_after_trim_qry_gaps[i-re_normalise_i]
   #    del ref_after_trim_qry_gaps[i-re_normalise_i]
   #    re_normalise_i += 1

   # ref_gap_positions_after_trim_qry_gaps = [i for i, letter in enumerate("".join(ref_after_trim_qry_gaps)) if '-' in letter]

   # # temp variables start with a copy of a and b as lists 
   # #  (to allow index-based deletions of elements)
   # qry_after_trim_qry_and_ref_gaps = qry_after_trim_qry_gaps
   # ref_after_trim_qry_and_ref_gaps = ref_after_trim_qry_gaps

   # re_normalise_i = 0  # since the sequences get progressively shorter with each element deletion, we must correct the index used to delete with, for every deletion

   # # corresponding to gaps in b: delete letters in a and b 
   # for i in ref_gap_positions_after_trim_qry_gaps:
   #    del qry_after_trim_qry_and_ref_gaps[i-re_normalise_i]
   #    del ref_after_trim_qry_and_ref_gaps[i-re_normalise_i]
   #    re_normalise_i += 1

   # # format back into strings
   # qry_trimmed = "".join(qry_after_trim_qry_and_ref_gaps)
   # ref_trimmed = "".join(ref_after_trim_qry_and_ref_gaps)

   # # # @TEST:passed:simulated the following seqs, a and b, and predicted exprected 
   # # #     correct output, and indeed output of trim_gaps_from_aligned_seqs() met expected
   # # qry_seq_aln = 'ATGTGC----TAA'
   # # ref_seq_aln = 'ATG--A--TTTGA'
   # # qry_expected_after_trim = 'ATGCTAA'
   # # ref_expected_after_trim = 'ATGATGA'

   # }} 1 alternative 2  {{ # wenping style: print only letters that both lists have in a nesteed list

   #qry_and_ref_arr = [qry_seq_aln, ref_seq_aln]

   qry_and_ref_arr = [(j,ref_seq_aln[i],i) for i,j in enumerate(qry_seq_aln) if (("-" != ref_seq_aln[i]) and ("-" != j))]

   qry_and_ref_arr = np.array(qry_and_ref_arr)

   qry_trimmed = "".join(list(qry_and_ref_arr[:,0]))
   ref_trimmed = "".join(list(qry_and_ref_arr[:,1]))
   qry_indices = list(qry_and_ref_arr[:,2])

   # }} alternative 2 

   #pdb.set_trace()
 
   return qry_trimmed, ref_trimmed, qry_indices


##########
## MAIN ## (still testing, cmd-args-mode not available yet)
##########
if __name__ == "__main__":

   #
   # PAIRWISE ALIGN
   #

   # @todo: make the s1_AGAP010815_RA and s2_AAEL001802_RA cmd args via the html

   #
   # AGAP010815 (gambiae)
   #     https://www.vectorbase.org/Anopheles_gambiae/Gene/Sequence?db=core;g=AGAP010815;r=3L:11202091-11206882;t=AGAP010815-RA
   s1_AGAP010815_RA = 'ATGTGGCAGTTCATAAGGTCACGAATATTAACGGTGATAATCTTCATAGGTGCTGCTCATGGGCTACTGGTTGTGGGTCCGAAATTTATACGGGCCAACCAGGAATACACTCTGGTGATCAGCAACTTTAACTCACAGCTAAGCAAAGTGGACCTGCTGTTAAAACTGGAAGGCGAAACTGATAATGGTTTAAGCGTTCTGAACGTTACCAAGATGGTTGACGTGCGACGTAATATGAACCGAATGATCAACTTCAATATGCCTGAGGATCTGACGGCTGGAAACTACAAAATAACTATCGATGGACAGCGTGGCTTCAGCTTTCACAAGGAGGCAGAGCTGGTGTATCTCAGCAAATCGATATCGGGGCTAATACAGGTCGATAAGCCCGTATTTAAACCTGGGGATACGGTGAACTTCCGTGTGATCGTGCTGGACACGGAGCTGAAACCGCCGGCGAGGGTCAAGTCGGTTTATGTAACTATACGAGATCCTCAGCGCAATGTGATTCGCAAATGGTCCACGGCAAAACTGTATGCCGGTGTGTTCGAGAGCGATCTACAGATAGCGCCTACTCCAATGCTCGGGGTCTGGAATATCTCGGTGGAGGTGGAAGGAGAAGAGCTTGTGTCAAAGACGTTTGAGGTGAAGGAGTACGTGTTGTCAACGTTCGACGTGCAGGTCATGCCATCGGTGATTCCACTGGAAGAGCATCAAGCTGTGAATCTTACAATCGAAGCGAACTATCACTTTGGTAAGCCAGTGCAAGGAGTGGCCAAGGTGGAGCTGTACCTAGACGACGATAAGCTAAAACTGAAAAAAGAGCTGACTGTGTACGGAAAGGGCCAGGTAGAGTTGCGCTTTGACAATTTTGCAATGGATGCGGATCAGCAGGATGTACCAGTGAAGGTGTCGTTCGTCGAGCAGTACACAAATCGTACGGTGGTCAAACAGTCACAAATCACGGTATATAGGTATGCGTACCGAGTAGAGTTGATAAAAGAGAGTCCACAGTTTCGTCCGGGACTCCCGTTCAAATGTGCGCTTCAGTTTACACACCATGATGGAACACCGGCTAAAGGCATTAGCGGTAAGGTAGAGGTATCCGATGTACGATTCGAAACGACAACAACGAGTGATAACGATGGATTGATTAAGCTCGAGCTGCAACCAAGTGAGGGTACTGAACAACTCAGTATTCACTTCAATGCTGTTGATGGATTCTTTTTTTATGAAGATGTGAATAAGGTAGAAACGGTTACGGATGCGTATATTAAACTGGAGCTGAAATCACCGCATCAAACGGAACAAATTGATGCGTTTCATGGTGACGTGCACGGAGCGCATGACATTCTTCGTGTACTATGTCATGTCAAAGGGCAATATCATCGATGCAGGATTCATGCGACCCAACAAGCAACCGAAGTACCTGTTGCAGCTGAACGCAACAGAAAAGATGATTCCGAGGGCGAAAATTCTCATCGCTACCGTAGCGGGCCGCACGGTGGTGTACGACTTCGCAGACCTCGATTTCCAAGAGCTTCGCAATAATTTTGATTTAAGCATTGACGAGCAAGAGATCAAGCCGGGACGACAAATCGAGCTGAGCATGTCTGGACGCCCAGGAGCGTACGTTGGGCTGGCCGCGTATGACAAAGCCTTGCTGCTTTTCAACAAGAACCACGACCTGTTCTGGGAGGACATTGGGCAGGTGTTTGATGGGTTCCATGCAATCAATGAGAACGAGTTTGACATATTCCACAGCTTGGGTCTGTTCGCCAGGACATTGGACGATATCTTGTTCGACAGTGCAAATGAAAAGACGGGGCGTAATGCACTGCAGTCAGGCAAGCCGATCGGCAAGCTGGTGTCGTATCGGACGAACTTCCAGGAATCGTGGTTGTGGAAAAATGTTTCCATCGGACGATCGGGAAGTCGCAAGTTGATCGAGGTAGTACCGGACACGACCACCTCCTGGTATCTGACGGGCTTCTCGATCGATCCCGTGTACGGGTTGGGTATCATCAAGAAGCCAATCCAGTTCACAACAGTCCAGCCGTTCTACATCGTAGAGAACTTACCATATTCAATCAAACGAGGCGAAGCGGTTGTGTTGCAGTTTACGCTGTTCAACAACCTTGGAGCGGAGTATATAGCCGATGTGACGCTGTACAATGTGGCCAACCAGACCGAGTTCGTCGGACGTCCAAATACGGATCTCAGCTACACCAAATCCGTGAGCGTTCCTCCAAAAGTTGGTGTGCCAATCTCGTTCCTCATCAAGGCCCGCAAGCTCGGCGAGATGGCGGTTCGTGTAAAGGCTTCGATAATGCTGGGACACGAAACGGACGCCCTGGAAAAGGTAATACGGGTGATGCCTGAAAGTTTGGTGCAGCCGAGAATGGATACACGCTTTTTCTGCTTCGACGATCACAAAAATCAAACGTTTCCGATCAACTTGGACATCAACAAGAAGGCCGACAGTGGATCGACAAAGATTGAGTTTCGACTAAATCCCAATTTGTTGACCACGGTCATCAAGAACCTGGACCATCTTCTCGGCGTTCCGACGGGATGTGGTGAGCAGAATATGGTCAAATTTGTTCCCAACATTTTGGTACTGGATTATTTGCATGCCATCGGGTCGAAAGAACAGCATCTAATCGACAAAGCTACGAATTTGTTGCGTCAAGGATATCAAAACCAGATGCGCTACCGTCAGACGGATGGTTCATTTGGTTTGTGGGAGACTACTAATGGTAGCGTGTTTCTCACCGCGTTCGTTGGCACATCGATGCAAACTGCAGTAAAATACATAAGCGATATTGATGCAGCAATGGTGGAGAAGGCATTGGATTGGTTAGCCTCGAAGCAGCATTTCTCGGGACGGTTTGACAAGGCCGGTGCAGAGTATCACAAAGAAATGCAAGGAGGGTTGCGCAATGGTGTGGCCCTCACATCATATGTGTTGATGGCATTGCTGGAGAATGACATTGCCAAAGCAAAGCACGCAGAGGTGATTCAAAAAGGAATGACCTATCTGAGCAATCAGTTTGGATCCATCAACAATGCATACGACCTATCGATAGCAACCTACGCGATGATGTTGAACGGACACACCATGAAGGAGGAGGCACTCAATAAGCTGATTGATATGTCTTTCATTGATGCTGATAAAAACGAACGGTTCTGGAACACAACGAATCCAATAGAAACCACCGCATATGCTCTGCTGTCGTTTGTGATGGCCGAGAAGTACACAGACGGTATACCGGTCATGAATTGGTTGGTGAATCAACGTTACGTTACCGGTAGCTTTCCGAGCACGCAAGACACGTTTGTGGGGCTGAAAGCGCTGACCAAAATGGCGGAAAAGATATCTCCGTCCCGAAACGACTACACCGTTCAACTGAAGTACAAGAAGAGTGCAAAATACTTCAAAATAAACTCGGAGCAAATTGATGTGGAAAACTTCGTGGATATACCGGAGGACACAAAAAAGCTCGAGATCAATGTGGGGGGCATTGGATTTGGGTTGTTAGAGGTGGTTTATCAATTTAATTTGAATCTCGTCAACTTTGAGAATAGATTCCAACTAGACCTGGAGAAACAGAACACAGGCTCTGACTACGAGCTGAGGCTGAAGGTCTGTGCCAGCTACATACCCCAGCTGACCGACAGACGATCGAACATGGCACTGATTGAGGTAACCTTACCGAGCGGTTACGTGGTTGATCGCAATCCGATCAGCGAGCAGACGAAGGTGAATCCGATTCAGAAAACTGAAATCCGTTACGGTGGCACTTCAGTCGTTTTATACTACGACAATATGGGCAGCGAGCGTAACTGTTTCACCCTGACCGCGTACAGACGCTTTAAGGTCGCATTGAAGCGTCCAGCGTATGTGGTTGTGTATGATTATTATAATACAAATCTGAACGCCATCAAAGTGTACGAAGTGGACAAGCAGAATTTGTGCGAAATCTGTGACGAAGAAGACTGTCCTGCAGAGTGCAAAAAATAG'
   #
   # AAEL001802 (aegypti)
   #     https://www.vectorbase.org/Aedes_aegypti/Gene/Sequence?db=core;g=AAEL001802;r=supercont1.43:685886-717122;t=AAEL001802-RA
   s2_AAEL001802_RA = 'ATGTCGGTATTCATACAAACGGACAAACCGGTGTATACCCCGGGAGATCTGATACGTTTTCGGGTAATCGTGGTGGATGCTGACACTAGACCTGTGACTAGTATTAAAACGGTAAATATAGCGATCGACGATTCTGCAAAAAATTCCATTCGAAAGTGGCCTTATGCCAAGTTGTTAAACGGCATCTTTGAGTCACAAGTGCAATTAGCTTCTTCGCCTGTTCTTGGCACCTGGATTATCAACGTAACAGCTTCCGACGACATCATTGTCACCAAACAGATAGAAGTTAAGGAATATGTGTTGCCAAAATTTTTCGTGAAAGTTTACCCTTCGGAGGTTCTATTGGGGAAAAATAAGAAGGTTTCTCTTACCTTAGATGCCTATTACACGTTCAAAGAACCCGTCGACGGCAATTACAAAGTTGAGTTATTTTTGGACCATACCAAGAGAAAGCCTGACTTCATAAAAAGTGATCGAATCACCGGTAAAACATCACTTGAGTTTCAATTGAAAAATGAAGTAGACATTGATGGCGACGAGCAGTACACTGATGTCACGGTTGAAGTTGAAGTTGTCGAGGCATTTTCTAATCGCACAGTTAGTATAACTGAGAATATTCCGATTTATCGTCAGCCTTATACCGTGACCCTTCTTCCATCTGCACCATCATTTCGACCAGGAGTTCCATTCAATGTACAAATAGTTGTGAAAGATCAGCTTGGACACCCTCCTGCCGAAGAAAAAGCGGCATCAATTGACCTTACTGTAGAGTTCCATTTGCCCATTGACAGTGACACCAAATCTATCACTGTAGATCTGGACGAGAAAGGAACAGGTCAGCTCACATTAGAGCCCCGCCCAGACGCCCAAGAACTGAAAGTGAACGCTACATATGACTCTCAACAATACGATGTAATTCACGATCCGATACATGGTTTCAGTTCGCAAAGTAAGCAGTACATCACAGTAACTCTGAATCCAAAATACTATAACAACATTAAAGTCGATAAGGACATCGTACTGGACATCTCCTGCACTGAAACAATGACGCACTTCTCGTACATCGTTGTCACCAGAGGAAACATAGTGGAAGCATCGAACGTTCCTGTCAGGATAAAAAAGAAACATTCTCTGAGATTGAAAATGACTTCAAAAATGTCTCCGGAGTCGAGGCTTCTAGTGTACTATACAAACAGGGAGTATCTCATCTTTGATGATATTGAGCTGAAGTTCGATTCGTTCAACAACGACTTCAAATTCGATTTGAACGATGATGAGTATTTTCCAGGGCAATCAGTTTATATCGATGTATACGCTTCAAAGGATTCATACGTTGCGTTCAGTGGAATCGATGAAAGTGTACTCCTGGTAGGCAAAGAGCGCCATGACTTCAACAAAGGAGATGTGCTCAAGGAACTCGCTCTTTACGGAGCAACAAATGATGCCGAGTTTGACTTGTTCCACGTAAGTTTCATGTCAAATGGTTTGATTATTCCAGTTAATGTATCTGTAACTCGCTCACAGAATGCACGATTTGGTACTCTACTAGGAAGGACTAGGCAGCAAGCGATTGAAATTCGAACTCAATTCCTAGAATCCTGGTTATGGAAATCCTTTTCCATGGATGGTCGAAACAACTTCAAAGCAATAGAAGACTCGGTTCCGGATACTATTACAACGTATCACGTGTCAGGATTTGCTTTAAGTCCAACACTAGGTCTTGGAGTAATCCAACAACCAGTGAGTTTCACCGTTCGTAAAAAATTCTACTTGGTTGCAAATTTGCCTTACTCGATCAAACGGGGTGAAGTGGCGTTGATTCAGGTTACCGTCTTCAACTTCCTAGGAAGCAGCATAACAACCGATGTGACGCTGTTCAATAAACGCGATGAAATTGAGTTTGTCGAGAATGCATCCACTAATAATACACATCGAACAAAGGCGGTAATTGTCCCGAATAACAATGGAAAATCTGTATCATTTATGGTGAAAGCAAAGAAATTAGGACAGATTGCGATCAAATTCCAGGCGGTAAACCTGCTGGAAACGGATGCATTGGAGCACATGTTACGAGTAACCCCAGAGAGCCATCGCTATGAGAAAAATGTAGCTCGATTCGTTGAGCTACCAAAGTTTGAGACGCAAACTTTCGATGTGAAGCTGGACATTCCCAAAAATATCGACGAGGGTTCTGCTCAAATCAAATTCACGTTAGACCCGGACATTTTGGGAACAGCCATCAGCAACCTAGACGGGTTGATCCGGAAACCCTTTGGATGTGGCGAACAAAATATGCTCCATTTTGTGCCAAATATAGTCGTTTTGGATTATCTTAACGAAACCAACACAGCGGCAGAAGATGTGAGGACCAAAGCGATAAATTTTCTTAGCAGCGGATATCAAAACCAGCTACGCTACAAACGTTCGGATGGGGCCTTCAGTGTCTGGGGACAATCGTATGCTGGCAGTACATTTTTGACGGCCTTTGTGGCGAAATCATTCAAAATAGCAGCCAAATACATTCAGGTGGATAAGTCTATAGTAGACGCGGCATTCGACTGGTTAGTGAAACAACAACAATCAGATGGGCGGTTCCCAGAAGTGGGGCAAGTATTCCAAGCAGATATGCAGGGTGGGCTTCGTAATAACGGTTTTGCGCTTACCGCGTATGTTCTGATCGCTTTTGCTGAAAATAAGGAAGTATACAGAAAATACCAATCACAACTGAACAAAACTACTAACTTCATAGCAGATAGACTTGCTAATATGGAGAATCCATACGACCTCTCGCTGTCCACTTATGCGTTGATGCTAACAAATCATGGCAAGCGCACCGAGTTTCTTCACAAATTAGTCGAAAAGTCGATATTTGACCGCAATCAAACTGAGAGATATTGGGACAGCAAACCAGTTGATATTGAAGTTGCTGGATATGCTCTATTGTCATACGTAGCTGCCGGTAAATTATTGGATGCAACGCCTATCATGCGGTGGCTCAACAAGCAGCGTTATGGTCTCGGAGGCTATCCTGGAACTCAGGAAACATTCGTTGGATTGAAAGCATTGGCAACGTTCGCTGCAAATGTAACTAGTAGGAGAAACGAATATACTGTAAGGATATTCTACGAACCAAATGGTCGACGAACATTCGACGTACACATGCACAATTCGTTTAATATTCAAGAGCTTGACATTCCTAATAACATCAGAAAAATGAAGGTGGAAGTTGAAGGCATCGGCAGAGGCTTCTTCCAAGTGGCATATCAGTACTATCAAAATATGCAGGTGGCTAAGCCCAGTTTCAGCATTACAATTAATCAGCTTAACACCACGACGGAACACATGCAGCAATTGGACGTGTGTGTGAAATACATACCAAAAGAGGCTTATCAAAAATCGAATATGGCTTTGGTGGAAATATTCTTGCCTAGTGGGCTTGTAGCAGACTCAGATGCCATTACGGACAAGACTGGAGGAATTCGAAGAATTGAAAGACGTTTTTCGGACACCTCAGTAGTTATATATTATGATAATTTGGACCCCGAAGACAAGTGCTTCCGAGTGACTGCTTATCGTCGGTATAAAATTGCATTGCATTTGCCATCATATATTATAGTTTATGATTATTATAATTTTGAGCGCTTTGCCATTCAAAAGTACGAAGGAAAGGTGCTGCAGCTCTGCGATATTTGTGAAGACGAGGACTGCGAAACTTTATCATGTCAAAATAGCTCGAAATTGGCAATAATGTAA'

   qry_seq_aln, ref_seq_aln = align_query_vs_reference(s1_AGAP010815_RA,s2_AAEL001802_RA)


   #
   # TRIM GAPS: 
   #     Next we want to remove all characters (trim) in each of the sequences that is aligned to a gap position**** in the top alignment
   #

   #  E.g. alignment:
   #     ATGTGC----TAA  qry_seq
   #     ATG--A--TTTGA  ref_seq
   #        ^^ ^^^^     gap positions****
   #  Perform trimming:
   #     ATGCTAA        qry_seq (after trimming)
   #     ATGATGA        ref_seq (after trimming)

   qry_seq_trimmed, ref_seq_trimmed, qry_seq_indices = trim_gaps_from_aligned_seqs( qry_seq_aln, ref_seq_aln )



   #
   # @test:if the alignment and trimmer functions did their job correctly there should be exactly equal len(seq) for qry and ref
   #
   try:
      assert len(ref_seq_trimmed) == len(qry_seq_trimmed)
   except AssertionError:
      # pdb.set_trace()
      raise AttributeError("qry_seq_trimmed and ref_seq_trimmed are not of equal length, unexpectedly (i.e. the qry and ref seqs after alignment + gap trimming), downstream analyses such as dnds.py will be erroneous...")


