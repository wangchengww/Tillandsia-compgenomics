#! /usr/env python

import HTSeq
import sys

fp = open(sys.argv[1])
gff_file = HTSeq.GFF_Reader(fp, end_included=True )

#gff_file = HTSeq.GFF_Reader( "subset_with_introns.gff", end_included=True )

transcripts={}

### Original part to get intron lengths

for feature in gff_file:
   if feature.type == "intron":
      transcript_id = feature.attr['Parent']
      if transcript_id not in transcripts:
         transcripts[ transcript_id ] = list()
         transcripts[ transcript_id ].append( feature )

for transcript_id in sorted( transcripts ):
   transcript_length = 0
   for intron in transcripts[ transcript_id ]:
      transcript_length += intron.iv.length
   print transcript_id, transcript_length, len( transcripts[ transcript_id ] )

### Extra part to add genes without introns

genes=[]

for feature in gff_file:
    if feature.type == "mRNA":
        gene_id = feature.attr['ID']
        if gene_id not in genes:
            genes.append(gene_id)

introns = []

for feature in gff_file:
    if feature.type == "intron":
        intron_id = feature.attr['Parent']
        if intron_id not in introns:
            introns.append(intron_id)

for gene_id in genes:
    if gene_id not in introns:
        print gene_id, "0", "0"

gff_file.close()
