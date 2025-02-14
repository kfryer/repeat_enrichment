# Double sided CASK ChAR-seq analysis steps
Follow these steps to count and calculate enrichment scores for repetitve RNAs contacting repetitive DNA sides annotated by CASK (double sided analysis). 


## Step 1: Count reads
Counts the number of reads in each input fastq and saves to a file. 


Inputs:

1. Path to input fastq (same as input #3 from Step 3 in cask_docs.md). Can handle gzip files. 

Steps:

1. Counts the number of reads in the fastq file 

Outputs:

1. {input_fastq_filename}.counts.txt

Run as:
```
cask_countreads_kf.sh {/path/to/input_fastq}
```
Example:
```
cask_countreads_kf.sh /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep1/dna.fastq.gz
```

## Step 2: Count all RNAs

Count the number of RNA reads for each ENSG/ENST. This will be used to normalize contact density in enrichment score calculations.

Inputs:

1. Full kmatch directory (Output directory for CASK read annotation steps)
2. Pairs dna.bed.gz file from charseq output (/sample_name/data/pairs/gencondeV29_hg38/all/dna.bed.gz)

Steps:

1. Pulls ENSGs from pairs bedfile for each read. Processes exon derived RNAs and intron derived RNAs independently. Reads with ENSTs are considered exon derived and those with just ENSGs are intron derived. Counts the number of reads in the file for each ENSG

Outputs:

1. {/path/to/kmatch_dir}/dnaside.bed.exons.countallRNAs.txt
2. {/path/to/kmatch_dir}/dnaside.bed.introns.countallRNAs.txt

 - Output fields: Count`\t`ENSG


Run as:

```
cask/scripts/count_allRNAs_kf.sh {kmatch_dir} {paired_dna_bedfile.gz}
```

Example:

```
cask/scripts/count_allRNAs_kf.sh /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep1/cen_kmatch/k25.minusBKG /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun/data/TFONM2_ES/pairs/gencondeV29_hg38/all/dna.bed.gz

```

## Step 3: Generate genome bin file

Generate genomic bins of various sizes based on chromosome coordinates.

Inputs:

1. Chromosome sizes file 
    Columns:
    -  chromosome
    - start_coordinate (should be 0 for each chromosome)
    - stop_coordinate

2. Output file prefix 

Steps:

1. Creates a bed file of genomic coordinates of equal sized bins for each chromosome.


Outputs:

1. Genome bin bed file for each bin size (currently includes 1kb, 10kb, 100kb). Modify bin sizes by changing the bin_size_dict keys and values in generate_gbins.py

    Output file is saved as:

    - "output_file_prefix" + binsize + "kb_genome_bins.bed"

Run as:
```
python generate_gbins.py my_chrom_sizes_file output_file_prefix
```

Example:

```
#Run:
python generate_gbins.py /oak/stanford/groups/astraigh/T2T/CURRENT/T2T-CHM13v2.0-hs1_jan2022/hs1.chromosomes.bed /home/groups/astraigh/kelsey/chm13_genome/v2/hs1_

#Outputs:
/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_1kb_genome_bins.bed
/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_10kb_genome_bins.bed
/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_100kb_genome_bins.bed

OR including only chromosomes 1-22,X, and M

#Run:
python generate_gbins.py /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/chrNameLength_1-22XM.sorted.txt /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_

#Outputs:
/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_1kb_genome_bins.bed
/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_10kb_genome_bins.bed
/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_100kb_genome_bins.bed

```

## Step 4: Convert DNAside BAM file to BED
Use bedtools to convert the dna.bam file output from charseq pipeline to a bed file. This bed file will then be intersected with genome bin bedfile to sort each read into a genomic bin.

Inputs:

1. dna.bam - output from charseq pipeline that contains alignment information for all reads (not just those that have an aligned RNA pair)
    ```/path/to/charseq_pipeline_output/data/sample_name/alignments/dna/bowtie_hg38/dna.bam```

Steps:

1. Use bedtools bamtobed function to convert dna.bam to bed file 

Outputs:

1. dna.all.bed - bed file containing genomic coordinates of all dna reads that aligned to the genome

Run as:

```
bedtools bamtobed -i dna.bam > dna.all.bed
```

Example:

```
bedtools bamtobed -i /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun/data/TFONM2_ES/alignments/dna/bowtie_hg38/dna.bam > /scratch/groups/astraigh/kelsey/charseq/diffchar/data/ESrep1/alignments/dna/hg38/dna.all.bed
```


## Step 5: Intersect dna bed file with genome bins bed file 

Use bedtools to assign each read to a genomic bin based on DNA side alignment. This is to include reads that are not classified as repetitive so that they can be included in background contact rate calculations.

Inputs:

1. DNA Bed file from Step 4 in this doc (created from dna.bam file output from charseq pipeline)

    ```/scratch/groups/astraigh/kelsey/charseq/diffchar/data/ESrep1/alignments/dna/hg38/dna.all.bed```

2. Genome bins bedfile from Step 3 in this doc

Steps:

Annotates each read with the genomic bin overlapping at least 50% of the read

Outputs:

1. Genome bin annotated dna bed file

Run as:

```
bedtools intersect -a dna.all.bed -b genome_bins.bed -wo > allDNA_gbins.bed
```

Example:
```
bedtools intersect -a /scratch/groups/astraigh/kelsey/charseq/diffchar/data/ESrep1/alignments/dna/hg38/dna.all.bed -b /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_{2}kb_genome_bins.bed -f 0.5 -wo > /scratch/groups/astraigh/kelsey/charseq/diffchar/data/ESrep1/alignments/dna/hg38/all_dna_{2}kb1-22XMgbins.bed


OR parallelize 

parallel -j4 'bedtools intersect -a /scratch/groups/astraigh/kelsey/charseq/diffchar/data/{1}/alignments/dna/hg38/dna.all.bed -b /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_{2}kb_genome_bins.bed -f 0.5 -wo > /scratch/groups/astraigh/kelsey/charseq/diffchar/data/{1}/alignments/dna/hg38/all_dna_{2}kb1-22XMgbins.bed' ::: ESrep1 ESrep2 DErep1 DErep2 ::: 10 100
```

## Step 6: Count DpnII sites in each repeat domain and in each genomic bin
Use bedtools to get fasta records for each repeat domain or genomic bin from coordinates provided in a bedfile. Use bbduk to annotate DpnII sites each each repeat domain or genomic bin. Use custom python script to summarize counts.

Inputs:

1. Repeat fasta (output from CASK Step 2 in cask_docs.md)
2. Genome fasta
3. Genomic bin bedfile (output from Step 4 in this doc)

Steps: 

1. Create dpnII.fa

    ```
    >dpnII
    GATC
    ```

2. Get fasta for domains

    Repeats (use repeats.fa generated by Step 1 in cask_docs.md)
    ``` 
    /path/to/CASK/output/fa/repeats.fa
    ```

    Genomic bins (use bedtools)
    ```
    #Run as:
    bedtools getfasta -nameOnly -fo gbins.fasta -fi genome.fasta -bed genome_bins.bed

    #Example:
    bedtools getfasta -nameOnly -fo /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_100kb_genome_bins.fa -fi /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -bed /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_100kb_genome_bins.bed

    ```
3. Annotate DpnII sites in each fasta record using bbduk.sh

    For repeats

    ```
    #Run as:
    bbduk.sh -Xmx200g threads=12 ref=dpnII.fa in=repeats.fa int=f k=4 overwrite=t outm=dpnII_annotated_repeats.fa maskmiddle=f rename=t

    #Example:
    bbduk.sh -Xmx200g threads=12 ref="/scratch/groups/astraigh/kelsey/centrochar/dpnii/dpnII.fasta" in="/scratch/groups/astraigh/kelsey/centrochar/cen_only/backup/cen_only/fa/repeats.fa" int=f k=4 overwrite=t outm="/scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_dpnII.fa" maskmiddle=f rename=t

    ```
    For genomic bins:

    ```
    #Run as:
    bbduk.sh -Xmx200g threads=12 ref=dpnII.fasta in=gbins.fasta int=f k=4 overwrite=t outm=gbins_dpnII.fa maskmiddle=f rename=t

    #Example
    bbduk.sh -Xmx200g threads=12 ref="/scratch/groups/astraigh/kelsey/centrochar/dpnii/dpnII.fasta" in="/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_100kb_genome_bins.fa" int=f k=4 overwrite=t outm="/scratch/groups/astraigh/kelsey/centrochar/dpnii/hg38_1-22XM_100kb_genome_bins_dpnII.fa" maskmiddle=f rename=t

    ```
4. Clean up bbduk output to get just the read ID and number of dpnII sites found
    
    Run for both repeats and gbins
    ```
    #Run as:
    grep '>' dpnII_annotated_output_bbduk.fa | awk -F ">|\t|=" '{OFS="\t"; print $2,$4}' > dpnII_info_cut.txt

    #Example:
    grep '>' /scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_dpnII.fa | awk -F ">|\t|=" '{OFS="\t"; print $2,$4}'> /scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_dpnII_cut.txt


    grep '>' /scratch/groups/astraigh/kelsey/centrochar/dpnii/hg38_1-22XM_100kb_genome_bins_dpnII.fa | awk -F ">|\t|=" '{OFS="\t"; print $2,$4}' > /scratch/groups/astraigh/kelsey/centrochar/dpnii/hg38_1-22XM_100kb_genome_bins_dpnII_cut.txt 

    ```

5. Map reptype dpnII counts to class (for repeats only)
    ```
    #Run as:
    join <(sort -k1,1 dpnII_info_cut.txt) <(sort -k1,1 /path/to/CASK/output/repeats.withClass.txt) > output_withClass_dpnII_cut.txt

    #Example
    join <(sort -k1,1 /scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_dpnII_cut.txt) <(sort -k1,1 /scratch/groups/astraigh/kelsey/centrochar/cen_only/backup/cen_only/repeats.withClass.txt) > /scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_withClass_dpnII_cut.txt
    ```
    
Outputs:

1. dpnII.fa - single entry fasta containing dpnII sequence
2. gbins.fa - fasta file where each read is a genomic bin and the its corresponding sequence (output from bedtools getfasta) 
3. gbins_dpnII.fa - gbins.fa with dpnII site annotation (output from bbduk)
4. repeats_dpnII.fa - repeats.fa file with dpnII site annotation (output from bbduk)
5. gbins_dpnII_cut.txt - File containing number of dpnII sites in each genomic bin(output from Step 4)
6. repeats_dpnII_cut.txt - File containing number of dpnII sites in each repeat domain (output from Step 4)
7. repeat_swithClass_dpnII_cut.txt - File containing number of dpnII sites in each repeat class (output from Step 5)



## Step 7: Calculate enrichment scores for each contact


Command Line Inputs:

1. sample name (denoted as i in run file)

Run File Inputs:

1. cask_parent_dir_base: base directory of annotated reads

2. cask_subdir: subdirectory of annotated reads (this is also the output dir)

3. chrom_sizes_file = Chromosome sizes file used to generate genomic bins in step 4

4. mybedfile = Repeat annotation bed file used to classify reads using CASK

5. genome_bin_size = size of genomic bins in bp (ex. ```10000``` for 10kb)

6. g = Name of genomic bin size (this was included to differentiate between 10kb bins including all chromosomes and 10kb bins including just chromosomes 1-22,X, and M) (ex ```10kb``` OR ```10kb1-22XM```)

7. gbin_parent_dir_base = base directory of binned paired reads (output from Step 4). This is the file path prefix prior to sample name

8. gbinpath_suffix = subdirectory and file name of binned paired reads. This is the filename and path suffix after the sample name 

9. pairsfile_dict = dictionary with sample_name as key and full path to pairs file from charseq pipeline

10. gbin_dpnII_f = file containing the number of dpnII sites found in each genomic bin (output from Step 6)

11. rep_dpnII_f = file containing the number of dpnII sites found in each repeat (output from Step 6)

9. RNA_types = a list of RNA types (ex. ```["exons", "introns"]```)

10. levels = a list of levels (ex. ```["class", "reptype"]```)

11. max_ag_chrs = The maximum number of chromosomes a signle amibivalence group can span. This is to exclude reads classified as ambivalence groups encompassing repeats on every chromosome which slows down enrichment calculations and complicates background modeling. 

Internal Inputs:

1. RNA and DNA ag files from CASK in 
    
    ```
    dna_cask_file:
    /cask_parent_dir_base/{sample_name}/cask_subdir/dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt

    rna_cask_file:
    /cask_parent_dir_base/{sample_name}/cask_subdir/rna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt
    ```

2. Genomic bin information for DNA reads (output from Step 5)
    ```
    gbin_file = gbin_parent_dir_base + gbinpath_suffix
    ```

3. ci2_pickle and unq_pickle (Generated by CASK from Step 4: Process Reads)
    ```
    /cask_parent_dir_base/{sample_name}/cask_subdir/{ag_side}.fastq.gz.ANNOTATED.ci2.aglut.pickle
    /cask_parent_dir_base/{sample_name}/cask_subdir/{ag_side}.fastq.gz.ANNOTATED.unq.aglut.pickle
    ```

4. Total RNA counts file (Generated by Step 2 in this doc: Count All RNAs)
    ``` 
    /cask_parent_dir_base/{sample_name}/cask_subdir/dnaside.bed.{RNA_type}.countallRNAs.txt
    ```

Steps:

1. Loads ambivalence group assignments from ag files (internal inputs #1)

2. Maps ambivalence groups to repeat or class types based on pickle file output from CASK. CASK output file still has repeat types annotated as their unique IDs. This step looks up the ambiv group ID and fills in the repeat types that are possible for that ambiv group.

3. Resolves reptypes by comparing possible repeat types according to unq k-mer and ci2 k-mers. This step takes assignments from both ci2 and unq k-mer databased (see CASK docs for explanation of ci2 and unq k-mers) to more specifically assign reads.

4. Counts the total number of RNA side reads annotated for each ambivalence group. This will be used in enrichment score calculation to account for RNA expression levels (highly abundant RNAs will make more DNA contacts just because they are more likely to be in proximity to chromatin).

5. Merges CASK annotations for RNA and DNA sides.

6. Combines annotations to get RNA and DNA side information for every read.
	
    
    -  **Both sides have repeat annotations:** Both RNA and DNA sides have repeat annotations (that are not NA). RNA side is either a multimapper or does not have an alignment annotation.

    - **Only RNA side has repeat annotation:** RNA side has repeat annotation and DNA side did not, DNA side assigned to genomic bin based on gbin file

    - **Only DNA side has repeat annotation:** DNA side has repeat annotation and RNA side did not, RNA side assigned gene name (based on charseq alignment and paired file output from charseq pipeline)

    - **Neither side has repeat annotation:** RNA assigned gene name & DNA side assigned genomic bin

    - **RNA side has multiple annotations:** RNA side has both repeat and gene (alignment) assignment. Multimapping reads (map qual < 255) are assigned the repeat annotation. Unique alignments (map qual >= 255) are assigned alignment annotation. DNA side reads with multiple annotations were just assigned the repeat annotation since the alignment annotation is just a genomic bin.

7. Counts number of reads for each RNA-DNA contact 
    
    Example:

    - 10 reads have DNA = [repeat_1,repeat_2] and RNA = 'SNORD3A'

    - 25 reads have DNA = [repeat_1,repeat_2] and RNA=[repeat_1]

8. Calculates an enrichment score (E score) for each contact. See calc_enrich_doublesided_v2.py for details

9. Models background contacts for each chromosome or group of chromosomes using all protein coding RNAs derived from all other chromosomes

10. Uses model of background contacts to predict Escore for each contact based on the N_tot of the RNA side and the domain length of the DNA side. Subtracts predicted_Escore from Escore to get residual score (R)



Outputs:

1. Resolved and repeat type/class annotated reads parquet for each level (class and reptype) and each ag_side (rna and dna).

    ```
    #Saved in: 
    cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File Names:
    dna_class.resolved_agmapped_reads.parquet
    dna_reptype.resolved_agmapped_reads.parquet
    rna_class.resolved_agmapped_reads.parquet
    rna_reptype.resolved_agmapped_reads.parquet
    ```
    
    Columns (_class_):
    
    - index
    - readID
    - class_aG_ci2: class level ambivalence group ID (based on ci2 k-mers)
    - ci2kcount_inread: number of ci2 k-mers in read
    - class_aG_unq: class level ambivalence group ID (based on unq k-mers)
    - unqkcount_inread: number of unq k-mers in read
    - ag_reptypes_ci2: repeat classes corresponding to ci2 ag ID (col 3)
    - ag_reptypes_unq: repeat classes corresponding to unq ag ID (col 5)
    - ag_reptypes_str: resolved repeat classes resulting from comparison between unq and ci2 k-mers
    
    Columns (_reptype_):
    
    - index
    - readID
    - reptype_aG_ci2: repeat type level ambivalence group ID (based on ci2 k-mers)
    - ci2kcount_inread: number of ci2 k-mers in read
    - reptype_aG_unq: repeat type level ambivalence group ID (based on unq k-mers)
    - unqkcount_inread: number of unq k-mers in read
    - ag_reptypes_ci2: repeat types corresponding to ci2 ag ID (col 3)
    - ag_reptypes_unq: repeat types corresponding to unq ag ID (col 5)
    - ag_reptypes_str: resolved repeat types resulting from comparison between unq and ci2 k-mers

2. Resolved ag stats - contains the number of ecah resolution that occured in resolve_reptypes function. See calc_enrichment.py for details.

    ```
    #Saved in: 
    /cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File Name:
    {ag_side}_resolved_ag_stats_{RNA_type}_{level}.txt
    
    #Example: 
    dna_resolved_ag_stats_exons_class.txt
    ```

    Columns:
    
    - res_type: resolving type (ex. subset; see resolve_reptypes in calc_enrich_doublesided.py for details)
    - number of reads resolved for each type

3. Combined annotations file - contains read level DNA and RNA annotation information for all reads including CASK annotated reads and aligned reads.

    ```
    #Saved in:
    /cask_parent_dir_base/{sample_name}/cask_subdir

    #File Name:
    cenDNAxcenRNA_withgenes_{g}gbins_combined_anno_{RNA_type}_{level}_df.tsv

    #Example: 
    cenDNAxcenRNA_withgenes_100kb1-22XMgbins_combined_anno_exons_class_df.tsv

    ```
4. Escore dataframe file contains counts, lengths, RNA total, E score, repeat types etc for each RNA-DNA contact. 

    ```
    #Saved in: 
    cask_parent_dir_base/{sample_name}/cask_subdir 

    #File Name:
    cenDNAxcenRNA_withalign_{g}gbin_ci2_and_unqv2.{level}.{RNA_type}.Escore

    #Exampe:
    cenDNAxcenRNA_withalign_100kb1-22XMgbin_ci2_and_unqv2.class.exons.Escore
    ```
    
    Columns:
    
    1. DNA_annot: string of repeat types/classes (sep by comma) or genomic bin if no repeat annotation
    2. DNA_chr: string of chromosomes occupied by DNA_annot (if more than one, chromosomes are seperated by comma)
    3. RNA_annot: string of repeat types/classes (sep by comma) or gene_name if no repeat annotation
    4. RNA_chr: string of chromosomes occupied by RNA_annot (if more than one, chromosomes are seperated by comma)
    5. RNA_type: class of RNA (ex. snoRNA, rep, cen)
    6. N_tot: total number of times RNA_annot was found in the entire dataset
    7. logNtot: np.log(N_tot)
    8. type: type of non-ag side annotation (ex. snoRNA)
    9. N_in: the number of reads corresponding to this contact (DNA_annot_y contacting RNA_annot_x)
    10. L_out: the length of regions outside DNA_annot (genome size - domain size)
    11. L_in: domain size
    12. raw_E: N_in/(N_tot-N_in)
    13. Escore: raw_E*(L_out/L_in)
    14. logEscore: np.log2(Escore)

4. Residual scores dataframe containing all the info from Escore df in output #3 as well as predicted Escores and R scores. 

    ```
    #Saved In:
    /cask_parent_dir_base/{sample_name}/cask_subdir 

    #File Name:
    cenDNAxcenRNA_withalign_{g}gbin_ci2_and_unqv2.{level}.{RNA_type}.E.GAM.resid

    #Example:
    cenDNAxcenRNA_withalign_100kb1-22XMgbin_ci2_and_unqv2.class.exons.E.GAM.resid
    ```

    Columns:

    - Columns 1-14 same as .Escore file above
    
    15. DNA_chrs: same as DNA_chr but as a list
    16. predicted_logEscore: predicted E score based on GAM of trans contacting protein coding RNAs
    17. R: Residual Score = logEscore - predicted_logEscore

5. Model fit stats file. R squared and number of contacts used in background model for each chromosome. 

     ```
    #Saved In:
    /cask_parent_dir_base/{sample_name}/cask_subdir 

    #File Name:
    cenDNAxcenRNA_withalign_{g}gbin_ci2_and_unqv2.{level}.{RNA_type}.E.GAM.resid.fitstats

    #Example:
    cenDNAxcenRNA_withalign_100kb1-22XMgbin_ci2_and_unqv2.class.exons.E.GAM.resid.fitstats
    ```

     Columns:
     
     1. index
     2. chromosome or group of chromosomes
     3. R squared
     4. Number of contacts used to generate model

Run as:
```
Edit run file (run_double_sided_calc_enrichmentv2_example.py) to update input parameters. Copy run file and calc_enrich_doublesided.py to same directory then run as:

python run_double_sided_calc_enrichmentv2_example.py {sample_name} 
```

Example:

```
python /home/groups/astraigh/kelsey/scripts/run_double_sided_calc_enrichmentv2_example.py ESrep1
```

## Next: Process and plot residual scores 
Follow centrochar_heatmap_example.ipynb

    












