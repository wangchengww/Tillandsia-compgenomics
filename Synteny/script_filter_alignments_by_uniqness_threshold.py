#!/usr/bin/env python

import sys

coord = open(sys.argv[1])
threshold = float(sys.argv[2]) # alignments overlapping more than the threshold will be removed

# Define a function to see if there's an overlap between two ranges
def isthereoverlap(interval1, interval2):
    """Does the range (start1, end1) overlap with (start2, end2)?"""
    start1 = interval1[0]
    end1 = interval1[1]
    start2 = interval2[0]
    end2 = interval2[1]
    return end1 >= start2 and end2 >= start1

# Define a function to calculate percentage of overlap
def perc_overlap(interval1, interval2):
    """
    Given [0, 4] and [1, 10] returns [1, 4]
    """
    if interval2[0] <= interval1[0] <= interval2[1]:
        start = interval1[0]
    elif interval1[0] <= interval2[0] <= interval1[1]:
        start = interval2[0]
    if interval2[0] <= interval1[1] <= interval2[1]:
        end = interval1[1]
    elif interval1[0] <= interval2[1] <= interval1[1]:
        end = interval2[1]
    overlap = [start,end]
    return (float(overlap[1] - overlap[0]) / (interval1[1] - interval1[0]))

# Store start and end positions of each alignment in a dictionnary
position_dict = {}  # I initiate the dictionary
i = 0
for line1 in coord.readlines()[1:]:
    i = i + 1
    line1 = line1.replace('\n','') # rm the return carriage
    splitted_line1 = line1.split('\t') # split the line only at first occurring tab
    chrom_ref = splitted_line1[0]
    chrom_aln = splitted_line1[4]
    start = int(splitted_line1[1])
    end = int(splitted_line1[1])+int(splitted_line1[2])
    position_dict[i] = [chrom_ref, chrom_aln, start, end]

# Now iterate over the each line of the coord file again and find out if given alignment has overlap with any other alignment from a different chromosome
coord = open(sys.argv[1])
outputfilename=sys.argv[1].replace(".coord",".{}uniqfilter.coord".format(threshold))
output=open(outputfilename,'w')
firstline = "Qname\tQpos\tQlength\tQstrand\tTname\tTpos\tTlength\tTstrand\tIdentity\tScore\n"
output.write(firstline)
i = 0 # keeps track of which line in the coord file we're at
k = 0 # number of non-unique alignments
l = 0 # number of unique alignments
total = len(position_dict)
for line2 in coord.readlines()[1:]:
    i = i + 1
    print("Doing alignment {} of {}".format(i,total))
    line2 = line2.replace('\n','') # rm the return carriage
    splitted_line2 = line2.split('\t')
    chrom_ref = splitted_line2[0]
    chrom_aln = splitted_line2[4]
    start = int(splitted_line2[1])
    end = int(splitted_line2[1])+int(splitted_line2[2])
    alignment = [start,end]
    j = 0 # number of times an alignment passed
    count = 0 # number of comparisons between different chromosomes
    v = 0 # to not filter out alignments with no overlap
    # for each line in the coord file, we iterate over all alignments in the dictionnary
    for line, pos in position_dict.items():
        # if current line is not the same alignment in dictionnary
        if i is not line:
            # if current alignment is from a different chromosome than the alignment from dictionnary
            if chrom_ref == pos[0]:
                if chrom_aln != pos[1]:
                    count = count + 1
                    # calculate the overlap and the proportion of overlap
                    isoverlapping = isthereoverlap(alignment, [pos[2], pos[3]])
                    # if there's no overlap, we move to the next
                    if isoverlapping == False:
                        j = j + 1
                    elif isoverlapping == True:
                        prop = perc_overlap(alignment,[pos[2], pos[3]])
                        # if the proportion of overlap exceeds the threshold, don't print the line
                        if (prop >= threshold):
                            k = k+1
                            break
                        else:
                            j = j + 1
                else:
                    count = count + 1
                    j = j + 1

    # If current alignment passed all iterations without encountering another alignment with overlap > threshold, print the line in the outfile
    if j == count and j != 0 and count != 0:
        l = l+1
        line_to_print = line2+"\n"
        output.write(line_to_print)
# Final message reporting a filtering summary
print("I'm done filtering! {} alignments had > {} overlap with another alignment. After filtering, {} alignments were kept.".format(k, threshold, l))
