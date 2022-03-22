te_gff="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/EDTA/leiboldiana/26_scaffold_full_anno/Tillandsia_leiboldiana_26_scaffolds.fasta.mod.EDTA.TEanno.no-unknown.gff3"
gff="/gpfs/data/fs71400/grootcrego/REFERENCES_TILLANDSIA/Tlei_assembly/assembly_26_scaffolds/Tillandsia_leiboldiana_v1.2.edited_allfeatures.26chrom.curated-orthologs.genes.gff"

bedtools intersect -a $gff -b $te_gff -wa -c > Tillandsia_leiboldiana_GENE-TE-intersection.counts.txt
bedtools intersect -a $gff -b $te_gff -wa -wb > Tillandsia_leiboldiana_GENE-TE-intersection.txt
