# Extraction of contigs removed for 'hitting multiple probes'


Requirements:   

1. Duplicates file generated using the `--keep-duplicates` flag in the `phyluce_assembly_match_contigs_to_probes` script step. 
1. dataset1.incomplete, dataset1.conf, and dataset1.fasta from `phyluce_assembly_get_match_counts` & `phyluce_assembly_get_fastas_from_match_counts`
2. MAFFT (should be installed with Phyluce).
2. BioPython (should be installed with Phyluce). 

Call:

`phyluce_dupe_contig_parser.py --output <OUTPUT_DIR> --contigs <TRINITY_ASSEMBLIES_DIR> --dupes <PATH/TO/DUPLICATES_FILE> --fasta <PATH/TO/dataset1.fasta> --conf <PATH/TO/dataset1.conf> --incomplete <PATH/TO/dataset1.incomplete>`

How it works: 

1.	Find all the contigs that hit locus2.	Align contigs against one another3.	Generate consensus sequence of alignment with mafft4.	Add consensus sequence to the large fasta file phyluce_assembly_get_fastas_from_match_counts5.	Modify the dataset1.incomplete and dataset1.conf files generated after `phyluce_assembly_get_match_counts` and `phyluce_assembly_get_fastas_from_match_counts` to correct 'missing' loci/taxa and add any loci that are new. 