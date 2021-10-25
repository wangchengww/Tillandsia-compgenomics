import sys
from Bio import SeqIO

fasta_sequences = SeqIO.parse(open(sys.argv[1]),'fasta')
window_size = float(sys.argv[2])
outputfilename="TE_content_per{}window_py.txt".format(window_size)
output=open(outputfilename,'w')

def chunks(seq, win, step):
    seqlen = len(seq)
    for i in range(0,seqlen,step):
        j = seqlen if i+win>seqlen else i+win
        yield seq[i:j]
        if j==seqlen: break

firstline = "chrom\tstart_window\tend_window\tperc_te\n"
output.write(firstline)

for fasta in fasta_sequences:
    name = fasta.id
    sequence = str(fasta.seq)
    start = 0
    end = window_size
    for subseq in chunks(sequence, window_size, window_size):
        masked_a = subseq.count("a")
        masked_c = subseq.count("c")
        masked_g = subseq.count("g")
        masked_t = subseq.count("t")
        masked_totalcount = masked_a + masked_c + masked_g + masked_t
        masked_perc = (masked_totalcount/window_size)*100
        line_to_write = name+"\t"+str(start)+"\t"+str(end)+"\t"+str(masked_perc)+"\n"
        print(line_to_write)
        start = start + window_size
        end = end + window_size
