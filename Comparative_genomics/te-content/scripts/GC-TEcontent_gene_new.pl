#!/usr/bin/perl

# file: gc_content.pl
# Calculate the GC content in a sliding window

use Getopt::Long;


while ($line = <>) {
  chomp;
  if ( $line =~ /^>/ ) {
     $id = $line ;
     $id =~ s/\n//g                         # remove the return carriage
#     print "$id";                          # extra line for neatness
  }else {
     $sequence = $line; 
#     print($sequence);
     print gc_content($sequence);       # gc_content is a local function
}
}

sub gc_content {
  my $seq = shift;                        # sequence
  my $win = shift;                        # window
  my $n_count = $seq =~ tr/Nn/Nn/;
  my $atgc_count = $seq =~ tr/ATGCatgc/ATGCatgc/;  # trick alert -- see manual entry for tr////
  my $atgcTE_count = $seq =~ tr/atgc/atgc/;  # trick alert -- see manual entry for tr////
  my $gc_count = $seq =~ tr/GCgc/GCgc/;  # trick alert -- see manual entry for tr////
  my $gcTE_count = $seq =~ tr/gc/gc/;  # trick alert -- see manual entry for tr////
#    my $others_count = $segment =~ tr/MRWSYKmrwsyk/MRWSYKmrwsyk/;  # trick alert -- see manual entry for tr////
  my $n_content = $n_count/length($seq);
  my $te_content = $atgcTE_count/$atgc_count;
  my $gc_content = $gc_count/$atgc_count;
  if ($atgcTE_count>=100) {
    my $gcTE_content = $gcTE_count/$atgcTE_count;
    print $id,"\t",$i+1,"\t",$gc_count,"\t",$gcTE_count,"\t",$atgc_count,"\t",$atgcTE_count,"\t",$n_count,"\t",$gc_content,"\t",$n_content,"\t",$te_content,"\t",$gcTE_content,"\n"}
  if ($atgcTE_count<100) {
    print $id,"\t",$i+1,"\t",$gc_count,"\t",$gcTE_count,"\t",$atgc_count,"\t",$atgcTE_count,"\t",$n_count,"\t",$gc_content,"\t",$n_content,"\t",$te_content,"\t","NA","\n"}
}
