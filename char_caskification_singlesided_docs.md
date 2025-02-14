# Single sided CASK ChAR-seq analysis steps
Follow these steps to count and calculate enrichment scores for non-repetitve (aligned) RNAs contacting repetitive DNA sides annotated by CASK (single sided analysis). 


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

## Step 3: Merge aligned RNA side info and CASK annotated DNA side info (unq and ci2)
Join ENSG info from pairs bedfile with CASK annotation based on read ID. 


Inputs:

1. Full kmatch_dir (same as input #1 for step 2 in this doc)
    
    Points to "/kmatch_dir/dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt" which is created by Step 5 in cask_docs.md
2. Pairs bedfile (same as input #2 for step 2 in this doc)

Steps:

1. Joins .ag.txt files with pairs bedfile to merge RNA alignment and DNA CASK classification

Outputs:

1. Full_kmatch_dir/dnaside.bed.exons.CASK_mergev2.txt
2. Full_kmatch_dir/dnaside.bed.introns.CASK_mergev2.txt

Run as:

```
cask/scripts/pair_dnabed_kfv2.sh {full_kmatch_dir} {pairs_dna_bedfile}
```

Example:

```
cask/scripts/pair_dnabed_kfv2.sh /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep1/cen_kmatch/k25.minusBKG /scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/ESrep1/paired.dna.bed.gz

```
## Step 4: Generate genome bin file

Generate genomic bins of various sizes based on chromosome coordinates.

Inputs:

1. Chromosome sizes file with columns [chromosome, start_coordinate, stop_coordinate]
Start coordinate for each chromosome should be 0
2. Output file path and prefix. Output file for 100kb will be saved as: 

    prefix + 100kb + _genome_bins.bed

Steps:

1. Creates a bed file of genomic coordinates of equal sized bins for each chromosome.


Outputs:

1. Genome bin bed file for each bin size (currently includes 1kb, 10kb, 100kb). Modify bin sizes by changing the bin_size_dict keys and values in generate_gbins.py

    Output file is saved as:

    "output_file_prefix" + binsize + "kb_genome_bins.bed"

Run as:
```
python generate_gbins.py my_chrom_sizes_file output_file_prefix
```

Example:

```
python generate_gbins.py /oak/stanford/groups/astraigh/T2T/CURRENT/T2T-CHM13v2.0-hs1_jan2022/hs1.chromosomes.bed /home/groups/astraigh/kelsey/chm13_genome/v2/hs1_

outputs:
/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_1kb_genome_bins.bed
/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_10kb_genome_bins.bed
/home/groups/astraigh/kelsey/chm13_genome/v2/hs1_100kb_genome_bins.bed

OR include only chromosomes 1-22,X,M

python generate_gbins.py /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/grch38_foralign/chrNameLength_1-22XM.sorted.txt /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_

outputs:
/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_1kb_genome_bins.bed
/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_10kb_genome_bins.bed
/home/groups/astraigh/kelsey/hg38/hg38_1-22XM_100kb_genome_bins.bed

```


## Step 5: Intersect pairs file with genome bins bed file 

Use bedtools to assign each read to a genomic bin based on DNA side alignment. This is to include reads that are not classified as repetitive so that they can be included in background contact rate calculations.

Inputs:

1. Pairs file output from charseq pipeline 

    /data/{sample_name}/pairs/gencondeV29_hg38/all/dna.bed.gz

2. Genome bins bedfile from Step 4

Steps:

Annotates each read with the genomic bin overlapping at least 50% of the read

Outputs:

1. Genome bin annotated pairs file

Run as:

```
bedtools intersect -a paired_dna.bed.gz -b genome_bin_file.bed -wo > output_file
```

Example:
```
bedtools intersect -a /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun/data/TFONM2_ES/pairs/gencondeV29_hg38/all/dna.bed.gz -b /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_10kb_genome_bins.bed -wo > /scratch/groups/astraigh/kelsey/charseq/diffchar/data/ESrep1/alignments/pairs/hg38/dna_10kbgbins_1-22XM.bed

OR parallelize 


parallel -j6 'bedtools intersect -a /oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun/data/TFONM2_ES/pairs/gencondeV29_hg38/all/dna.bed.gz -b /home/groups/astraigh/kelsey/hg38/hg38_1-22XM_{1}kb_genome_bins.bed -wo > /scratch/groups/astraigh/kelsey/charseq/diffchar/data/ESrep1/alignments/pairs/hg38/dna_{1}kbgbins_1-22XM.bed' ::: 10 100 1000
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

1. sample name

Run File Inputs:

1. cask_parent_dir_base: base directory of annotated reads

2. cask_subdir: subdirectory of annotated reads (this is also the output dir)

3. gbin_parent_dir_base: base directory of binned paired reads (output from Step 4). This is the file path prefix prior to sample name

4. gbin_suffix: subdirectory and file name of binned paired reads. This is the file path suffix after the sample name 

5. chrom_sizes_file: Chromosome sizes file used to generate genomic bins in step 4

6. mybedfile: Repeat annotation bed file used to classify reads using CASK

7. gbin_dpnII_file: file containing the number of dpnII sites found in each genomic bin (output from Step 6)

8. rep_dpnII_file: file containing the number of dpnII sites found in each repeat (output from Step 6)

9. RNA_types: a list of RNA types (ex. ```["exons", "introns"]```)

10. levels: a list of levels (ex. ```["class", "reptype"]```)

11. ag_side: which side of charseq read was annotated with CASK (ex. dna)

12. genome_bin_size: size of genomic bins in bp (ex. ```10000``` for 10kb)

13. g: Name of genomic bin size (this was included to differentiate between 10kb bins including all chromosomes and 10kb bins including just chromosomes 1-22,X, and M) (ex ```10kb``` OR ```10kb1-22XM```)

14. max_ag_chrs: The maximum number of chromosomes a signle amibivalence group can span. This is to exclude reads classified as ambivalence groups encompassing repeats on every chromosome which slows down enrichment calculations and complicates background modeling. 

Internal Inputs:

1. ci2_pickle: cask_parent_dir_base + i + cask_subdir + ag_side + ".fastq.gz.ANNOTATED.ci2.aglut.pickle"
   
   (Generated by CASK from Step 4: Process Reads)
2. unq_pickle: cask_parent_dir_base + i + cask_subdir + ag_side + ".fastq.gz.ANNOTATED.unq.aglut.pickle"

    (Generated by CASK from Step 4: Process Reads)
3. tot_RNA_counts_file: cask_parent_dir_base + i + cask_subdir + "dnaside.bed." + RNA_type +".countallRNAs.txt"

    (Generated by Step 2 in this doc: Count All RNAs)

Steps:

1. Loads ambivalence group assignments from ag_file where i is the sample name input with command and r is the RNA type from RNA_types list:

    ```
    #Saved in:
    cask_parent_dir_base/{sample_name}/cask_subdir
    
    #File name:
    {ag_side}side.bed.{RNA_type}.CASK_mergev2.txt

    #Example:
    dnaside.bed.exons.CASK.mergev2.txt
    
    ```

2. Maps ambivalence groups to repeat or class types based on pickle file output from CASK. CASK output file still has repeat types annotated as their unique IDs. This step looks up the ambiv group ID and fills in the repeat types that are possible for that ambiv group.

3. Resolves reptypes by comparing possible repeat types according to unq k-mer and ci2 k-mers. This step takes assignments from both ci2 and unq k-mer databased (see CASK docs for explanation of ci2 and unq k-mers) to more specifically assign reads.

4. Filters gbin file for CASK annotated reads. This step removes all reads that were assigned by CASK so that no reads are duplicated when the background genomic bin reads are combined with repeat annotated reads. 

5. Counts number of reads for each RNA-DNA contact (ex. DNA = [repeat_1,repeat_2] RNA = 'SNORD3A')

6. Calculates an enrichment score (E score) for each contact. See calc_enrichment_v2.py for details

7. Models background contacts for each chromosome or group of chromosomes using all protein coding RNAs derived from all other chromosomes

8. Uses model of background contacts to predict E score for each contact based on the N_tot of the RNA side and the domain length of the DNA side. Subtracts predicted_Escore from Escore to get residual score (R)

Outputs:

1. Counts of each RNA-DNA contact (after ambivalence group mapping, reptypes resolving, and combining with genomic_bin info).
    ```
    #Saved in:
    /cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File name:
    dnaside.bed.{RNA_type}_{level}.CASK_mergev2_{g}gbin_counts.txt

    #File name if combine_gbin=False:
    dnaside.bed.{RNA_type}_{level}.CASK_mergev2_counts.txt
    ```

2. Resolving stats file. This contains the number of ecah resolution that occured in resolve_reptypes function. See calc_enrichment.py for details. 

    ```
    #Saved in:
    /cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File name:
    resolved_ag_stats_{RNA_type}_{level}.txt
    ```


3. Escore dataframe as csv which contains counts, lengths, RNA total, E score,repeat types etc for each RNA-DNA contact. Saved as:

    ```
    #Saved in:
    /cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File name:
    dna.fastq.gz.ANNOTATED.ci2_and_unqv2.{level}.{RNA_type}.{g}.Escore

    #Example:
    dna.fastq.gz.ANNOTATED.ci2_and_unqv2.class.exons.100kb1-22XM.Escore

    ```
    
    Columns:
    
    1. ag_reptypes - list of possible repeat types (or repeat classes)
    2. ag_reptypes_str - possible repeat types/classes as a string (comma sep)
    3. ag_chrs - list of chromosomes the possible ag_reptypes span
    4. ag_chrs_str - list of chromosomes as a string (comma sep)
    5. gene_ID - ENSG ID for non-ag side (ex. for RNA side that is found contacting DNA ags in col 1)
    6. gene_name - gene name for non-ag side
    7. chrom - chromosome of non-ag side
    8. type - type of non-ag side annotation (ex. snoRNA)
    9. N_in - the number of reads corresponding to this contact (ex. DNA from ag_reptypes contacting RNA from gene_name)
    10. N_tot - total number of times gene_ID was found in the entire dataset
    11. dpnII - number of dpnII sites found in DNA domain
    12. logNtot - np.log(N_tot)
    13. L_out - Ltot (genome size) - L_in (domain size)
    14. L_in - domain size
    15. raw_E - N_in/(N_tot-N_in)
    16. Escore - raw_E*(L_out/L_in)
    17. logEscore - np.log2(Escore)

4. Residual scores dataframe containing all the info from Escore df in output #3 as well as predicted Escores and R scores. Saved as:

    ```
    #Saved in:
    /cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File name:
    {sample_name}_{RNA_type}_{level}_ci2_and_unqv2_{g}.E.GAM.resid

    #Example:
    ESrep1_exons_class_ci2_and_unq_100kb1-22XM.E.GAM.resid
    ```

    Columns:

    - Columns 1-17 same as .Escore file above
    
    18. predicted_logEscore = predicted E score based on GAM of trans contacting protein coding RNAs
    19. R score = logEscore - predicted_logEscore

5. Model fit stats file. R squared and number of contacts used in background model for each chromosome. Saved as:

    ```
    #Saved in:
    /cask_parent_dir_base/{sample_name}/cask_subdir 
    
    #File name:
    {sample_name}_{RNA_type}_{level}_ci2_and_unqv2_{g}.E.GAM_fit.stats

    #Example:
    ESrep1_exons_class_ci2_and_unq_100kb1-22XM.E.GAM_fit.stats
    ```

     Columns:
     
     1. index
     2. chromosome or group of chromosomes
     3. R squared
     4. Number of contacts used to generate model

Run as:
```
Edit run file (run_calc_enrichmentv2_example.py) to update input parameters. Copy run file and calc_enrichment_v2.py to same directory then run as:

python run_calc_enrichmentv2_example.py sample_or_rep_name 
```

Example:

```
python /home/groups/astraigh/kelsey/scripts/run_calc_enrichmentv2_example.py ESrep1
```

## Next: Process and plot residual scores 
Follow centrochar_heatmap_example.ipynb

    












