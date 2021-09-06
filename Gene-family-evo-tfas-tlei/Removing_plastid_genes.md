# Filtering out genes outside of the main genome

It is not so easy to separate plastid and mitochondrial genes from nuclear genes, as sometimes proteins active in the organelles are stored in the nucleus. I tried first by simply removing all genes with an annotation containing the words "chloroplast", "mitochondrial" or "ribosomal", but this seemed a leaky method. Genes from the photosystem, for example, are located in the chloroplast but this is not mentioned in the annotation. Therefore, I was not actively removing them from my set of genes.

It is key to properly filter out non-nuclear genes, since they will otherwise come out as the main large gene-families in our analysis. So I decided to try another approach - use uniprot annotations to properly remove genes located in the organelles.

To do this, I had to return to the annotation file from Blast2Go to obtain gene IDs from other assemblies that blasted to my gene models. I can then search these genes against a list of typical chloroplast genes.

I extracted the gene names of all annotated A. comosus genes located on the chloroplast. To do this, I downloaded the GFF annotation file of the F153 assembly from NCBI (here)[https://www.ncbi.nlm.nih.gov/genome/?term=txid4615].

I obtained the names of the plastid genes with the one-liner:
`grep "NC_026220" GCF_001540865.1_ASM154086v1_genomic.gff | awk '$3 == "gene" {print $0}' | cut -f 9 | sed 's/;/\t/g' | cut -f 3 | sed 's/=/\t/g' | cut -f 2 | sort -u > plastid_genes_Acomosus_F153.txt`

Exactly the same was done for mitochondrial genes from Oryza sativa (found (here)[https://www.ncbi.nlm.nih.gov/genome/10])

I also extracted the genes we are examining from the original Blast2Go gff file for both species:

Using these files, I accumulated the genes from mitochondrial, plastid and ribosomal origin in the following way:

    grep -w -f plastid_genes_Acomosus_F153.txt Tfas_genes_blast2go_annotations.txt | cut -f 1 > Tfas_mito_plastid_ribo_genes.txt
    grep -w -f mitochondrial_genes_oryza_IRSGP1.txt Tfas_genes_blast2go_annotations.txt | cut -f 1 >> Tfas_mito_plastid_ribo_genes.txt
    grep "ribosomal" Tfas_genes_blast2go_annotations.txt | cut -f 1>> Tfas_mito_plastid_ribo_genes.txt

This resulted in 353 genes in Tfas and 418 genes in Tlei.
I then combined both files:

    `cat Tfas_mito_plastid_ribo_genes.txt Tlei_mito_plastid_ribo_genes.txt | sort -u > mito_plastid_ribo_genes_to_remove.txt`

Note: I realized that this approach still left plastid genes out. So additionally, I blasted all the T.fasciculata genes to the fasta file of A.comosus plastid genes. To obtain the A.comosus fasta sequences of the plastid genes, I used getfasta from bedtools with the gff file. This is because the original fastafile on NCBI doesn't contain all genes.

    `grep "NC_011033.1" GCF_001433935.1_IRGSP-1.0_genomic.gff | awk '$3 == "gene" {print $0}' > mitochondrial_genes_oryza_IRSGP1.gff
    gff2bed < mitochondrial_genes_oryza_IRSGP1.gff > mitochondrial_genes_oryza_IRSGP1.bed
    bedtools getfasta -fi motichondrion_full_seq_Oryza_IRGSP1.fasta -bed mitochondrial_genes_oryza_IRSGP1.bed -name > mitochondrial_genes_oryza_IRSGP1.fasta`

Blast was run as following:

    `/apps/ncbiblastplus/2.10.0/bin/blastn -query plastid_genes_Acomosus_F153.fasta -subject   Tillandsia_fasciculata_v1.transcripts.in_orthogroups.fasta -outfmt "6 qseqid sseqid pident length qcovs qlen slen gapopen mismatch evalue" -evalue 1e-10 -out blastn_Acomosus_chloroplast_Tfas.out`

I then simply combined the identified genes from the blast run and the name search, as the lists only partially overlap:

    `grep -w -f plastid_genes_Acomosus_F153.txt Tfas_genes_blast2go_annotations.txt | cut -f 1 | sort -u > namesearch_genes
    cut -f 2 blastn_Acomosus_chloroplast_Tfas.out | sort -u > blast_genes
    cat blast_genes namesearch_genes | sort -u > plastid_genes_Acomosus_F153_blast_and_namesearch.to_remove.txt`

Then, I combined all mitochondrial and chloroplast genes (of both species) into one file:

    `cat Tfas_plastid_genes_Acomosus_F153_blast_and_namesearch.to_remove.txt Tlei_plastid_genes_Acomosus_F153_blast_and_namesearch.to_remove.txt Tfas_mito_genes_Oryza_IRSGP1_blast_and_namesearch.to_remove.txt Tlei_mito_genes_Oryza_IRSGP1_blast_and_namesearch.to_remove.txt | sort -u > all_mito_plastid_genes_to_remove.txt`

Lastly, I added ribosomal genes with the grep function as showed above. The total set of mitochondrial, chloroplast and ribosomal genes to remove from both species contained 852 entries.

I don't just want to remove these genes, but the entire orthogroup they belong to. So, to get hold of that information, I searched for the genes in the original orthogroups file, extracted the orthogroups:

    `grep -f mito_plastid_ribo_genes_to_remove.txt orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.txt | cut -f 7 | sort -u > mito_plastid_ribo_OGs_to_remove.txt`

This lead to 404 orthogroups to remove, which is a much more select group than by plainly using the search function.

    `grep -v -f mito_plastid_ribo_OGs_to_remove.txt ../orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.txt > orthogroups_Tfas_Tlei_Acom.per_gene.with_functional_info.no_TEs.size_corrections.no_plastid-mito-ribo.blastandsearch.txt`

The final list contained 69,229 genes and 18,697 orthogroups.
