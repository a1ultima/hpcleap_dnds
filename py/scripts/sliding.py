
import matplotlib.pyplot as plt
import numpy as np
import pickle


# dnds_slide_dict has overlapping windows of dnds calculated, windows overlap by stepLength bps, and are uniformly windowLength wide
with open("py/data/dnds_slide_dict.p","r") as fi:
    dnds_slide_dict=pickle.load(fi)

window_intervals = dnds_slide_dict.keys() # @todo: I dont think sorting the windows will make a difference to final result, but it will make it slower, @todo: test this just in case
#window_intervals = sorted(dnds_slide_dict.keys())

# @DONE: Q: would finding the max() be faster than sorting then indexing? A: yes

#
# vectorize and find max of the window_intervals
#
window_intervals_np = np.asarray(window_intervals)
max_window_position = np.amax(window_intervals) # e.g. 243 @done: a whole order of magnitude faster than: max_window_interval = max(window_intervals, key=lambda x: x[1])

# @todo: maybe this means we should  np.arrays() in the rest of the code in dnds.py and changes.py

# @todo: deal with non-finite values, e.g. NaN and Inf, what can we do about these?

#
# Initialise matrix with NaN, each row is the values for a specific window, these will cascade 
#       and overlap to various degrees depending on stepLength and windowLength

# ----     window 1
#  ----    window 2
#   ----   ...
# ^^^^^^^  take average along each column
# 
overlap_matrix = np.empty((len(window_intervals),max_window_position+1))  # @todo: are you sure it's +1? initialize empty matrix, note: entries are not actually NaN, just near-zero
overlap_matrix[:] = np.NAN
#overlap_matrix[:] = 0.0

for window_i,window in enumerate(window_intervals):

    start = window[0] # e.g. 0
    end   = window[1] # e.g. 49 

    # in the i-th row, fill all elements from the window[0]-th to window[1]-th with the dN/dS value for this window
    overlap_matrix[window_i,start:end+1] = dnds_slide_dict[window]['dNdS'] # @todo: are you sure it's +1? test, keep in mind for these indices it does -1 for the "to" part


# columwise averages (e.g. avg all values in each column)
nan_masker = ~np.isfinite(overlap_matrix) # boolean matrix, True if element is finite, False if element is Inf or NaN

overlap_matrix_masked = np.ma.masked_array(overlap_matrix,mask=nan_masker)

overlap_matrix_avg = overlap_matrix_masked.mean(axis=1)

overlap_matrix_avg.mean() # sanity: 0.15584966052233765 (whole seq average: 0.152307100775) not bad!

plt.plot(overlap_matrix_avg)

plt.show()
#avg_matrix = overlap_matrix.mean(axis=1)

#   window 

# mask the non-finite values and take means

# for i in 
