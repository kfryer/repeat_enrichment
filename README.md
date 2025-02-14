# Repeat enrichment analysis of ChAR-seq data

CASK (classification of ambivalent sequences using k-mers) to annotate RNA and DNA side ChAR-seq (chromatin associated RNA-sequencing) reads to identify RNAs enriched at repetitive regions of the genome

## Analysis Overview

1. Process raw ChAR-seq reads

2. Classify repetitive ChAR-seq reads using CASK 
    - Generate k-mer databases for repeat types
    - Identify repeat k-mers in sequencing reads 
    - Annotate sequencing reads with repeat annotation(s)
3. Count contact frequency for each RNA-DNA pair

4. Calculate enrichment values for each RNA-DNA contact
    
    Enrichment score normalizes read counts by RNA abundance, contact domain length, and DpnII site frequency to identify RNAs associated with repetitve regions of the genome above non-specific background interaction levels. 
   
    - Single sided analysis: enrichment of non-repetitive RNAs at repetitive genomic regions. This analysis uses CASK to annotate DNA side reads that are repetitive and transcriptome alignment of RNA side reads to annotate non-repetitive RNAs.
    - Double sided analysis: enrichment of repetitive RNAs at repetitive genomic regions. This analysis uses CASK to annotate both DNA and RNA side reads.

5. Data visualization

    - Read abundance comparison across repeat types
    - Summary of classes of RNAs enriched at repetitive regions
    - Heatmap of non-repetitive RNAs enriched at repetitive regions (class and repeat type level)
    - Heatmap of repetitive RNAs enriched at repetitive regions (class and repeat type level)




