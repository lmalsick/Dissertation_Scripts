This is where the scripts for the RNAseq pipeline live

You must convert the GTF file with this line:
sed 's/\(gene_id "[^"]*\).*/\1"/' GCF_021234435.1_CSHL_Jam_final_genomic.gtf > Fixed_AJ6file.gtf

Use the fixed GTF for subread/featurecounts
There is a space in the 9th column that causes it to fail, the rest of the scripts are updated accordingly
