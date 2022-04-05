#!/usr/bin/env python

import sys

coord = open(sys.argv[1])
threshold = sys.argv[2]
outputfilename=sys.argv[1].replace(".coord","{}uniqfilter.txt".format(threshold))
output=open(outputfilename,'w')

# Store all alignments as ranges from start to end position in a dictionnary
position_dict = {}  # I initiate the dictionary
i = 0
for line1 in coord.readlines()[1:]:
    i = i + 1
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line only at first occurring tab
    start = int(splitted_line1[1])
    end = int(splitted_line1[1])+int(splitted_line1[2])
    position_dict[i] = range(start, end)

# Now iterate over the each line of the coord file again and find out if given alignment has overlap with any other alignment
coord = open(sys.argv[1])
i = 0 # keeps track of which line in the coord file we're at
k = 0 # number of non-unique alignments
l = 0 # number of unique alignments
for line2 in coord.readlines()[1:]:
    i = i + 1
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t')
    start = int(splitted_line2[1])
    end = int(splitted_line2[1])+int(splitted_line2[2])
    alignment = range(start,end)
    j = 0
    # for each line in the coord file, we iterate over all alignments in the dictionnary 
    for line, pos in position_dict.items():
        # if current line is not the same alignment in dictionnary
        if i is not line: 
            # calculate the overlap and the proportion of overlap
            overlap = [n for n in alignment if n in pos]
            prop = float(len(overlap)) /float(len(alignment))
            # if the proportion of overlap exceeds the threshold, don't print the line
            if (prop >= threshold):
                k = k+1
            else:
                j = j + 1
    # If current alignment passed all iterations without encountering another alignment with overlap > threshold, print the line in the outfile
    if j == (len(position_dict.keys()) - 1):
        l = l+1
        line_to_print = line2+"\n"
        output.write(line_to_print)
# Final message reporting a filtering summary
print("I'm done filtering! {} alignments had > {} overlap with another alignment. After filtering, {} alignments were kept.".format(k, threshold, l))

