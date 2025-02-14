from calc_enrich_doublesided import *
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
import ast
import seaborn as sns
from pygam import LinearGAM, s, f
import scipy.stats as stats
import gc
import sys


cask_parent_dir_base = "/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/"
cask_subdir = "/cen_kmatch/k25.minusBKG/"

chrom_sizes_file = "/oak/stanford/groups/astraigh/T2T/CURRENT/T2T-CHM13v2.0-hs1_jan2022/hs1.chromsizes.txt"
mybedfile = "/home/groups/astraigh/kelsey/chm13_genome/v2/censat.clean.chrclass.bed"


genome_bin_size = 100000
g = '100kb1-22XM'
gbin_parent_dir_base = "/scratch/groups/astraigh/kelsey/charseq/diffchar/data/"
gbinpath_suffix = "/alignments/dna/hg38/all_dna_" + g + "gbins.bed"


pairsfile_dict = {'ESrep1':'/oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun/data/TFONM2_ES/pairs/gencondeV29_hg38/all/dna.bed.gz','ESrep2':'/oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar2/CL_2020-05-06_rerun/data/TMF2_Rep2B_ES/pairs/gencondeV29_hg38/all/dna.bed.gz','DErep1':'/oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar1/NOVAseq_06-27-2019/CL_2020-05-06_rerun/data/TGNQ5_DE/pairs/gencondeV29_hg38/all/dna.bed.gz', 
'DErep2':'/oak/stanford/groups/astraigh/differentiation_paper_data/charles/libraries/novchar2/CL_2020-05-06_rerun/data/TMF3_Rep2A_DE/pairs/gencondeV29_hg38/all/dna.bed.gz'} 

gbin_dpnII_f = "/scratch/groups/astraigh/kelsey/centrochar/dpnii/hg38_1-22XM_100kb_genome_bins_dpnII_cut.txt"
rep_dpnII_f=  "/scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_withClass_dpnII_cut.txt"


RNA_types=["exons","introns"]

levels=["reptype","class"]


max_ag_chrs = 5

####################################################################################################


i = sys.argv[1]


print("Starting:", i)
i_cask_parent_dir = cask_parent_dir_base + i + cask_subdir
gbin_file = gbin_parent_dir_base + i + gbinpath_suffix
dna_cask_file = i_cask_parent_dir + "dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt"
rna_cask_file = i_cask_parent_dir + "rna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt"
pairs_bedfile = pairsfile_dict[i]
for l in levels:
	print("Handling ags for:", i, l)
	resolved_dna_df, dna_res_stats_df  = load_ag_dat(dna_cask_file, l, 'dna', save=True)
	resolved_rna_df, rna_res_stats_df  = load_ag_dat(rna_cask_file, l, 'rna', save=True)
	#calculate Ntot for rep RNAs:
	RNAag_totcounts_dict = count_rep_RNAs(resolved_rna_df)
	ccmerged_df = ccmerge(resolved_dna_df, resolved_rna_df)


	for r in RNA_types:
		print("Starting calculations for", i, l ,r)

		i_Edf_file = i_cask_parent_dir + "cenDNAxcenRNA_withalign_" + g + "gbin_ci2_and_unqv2." + l + "." + r + ".Escore"

		if os.path.exists(i_Edf_file):
			print(f"Loading exisiting Escore df from {i_Edf_file}")
			Edf = pd.read_csv(i_Edf_file, sep="\t", header=0, index_col=0) 
		else:
			combined_anno_df_file = i_cask_parent_dir + "cenDNAxcenRNA_withgenes_" + g + "gbins_combined_anno_" + r + "_" + l + "_df.tsv"
			if os.path.exists(combined_anno_df_file):
				print("Loading combined annotations file:", combined_anno_df_file)
				combined_annots = pd.read_csv(combined_anno_df_file, sep="\t", header=0, index_col=0)
			
			else:
				combined_annots = combine_annotations(ccmerged_df, gbin_file, pairs_bedfile, r)

				comb_out = combined_annots.copy()
				comb_out.to_csv(combined_anno_df_file, sep="\t")

			combined_ann_counts = count_contacts(combined_annots)
			Edf = calculate_Escore(r, i_cask_parent_dir, l, chrom_sizes_file, mybedfile, combined_ann_counts, RNAag_totcounts_dict, g, genome_bin_size, rep_dpnII_f, gbin_dpnII_f, save=True)


		if r == 'exons':
			Rdf, model_fit_df = compute_GAM_residual(Edf, Edf, max_ag_chrs)
		else:
			exon_Edf_file = i_cask_parent_dir + "cenDNAxcenRNA_withalign_" + g + "gbin_ci2_and_unqv2." + l + ".exons.Escore"
			
			if os.path.exists(exon_Edf_file):
				print("Loading existing Escores for exons from ", exon_Edf_file)
				ex_Edf=pd.read_csv(exon_Edf_file, sep="\t", header=0, index_col=0)
			else:
				ex_combined_anno_df_file = i_cask_parent_dir + "cenDNAxcenRNA_withgenes_" + g + "gbins_combined_anno_exons_" + l + "_df.tsv"
				if os.path.exists(ex_combined_anno_df_file):
					print("Loading combined annotations file:", ex_combined_anno_df_file)
					ex_combined_annots = pd.read_csv(ex_combined_anno_df_file, sep="\t", header=0, index_col=0)
					ex_combined_ann_counts = count_contacts(ex_combined_annots)
				else:
					ex_combined_annots = combine_annotations(ccmerged_df, gbin_file, pairs_bedfile, "exons")
					ex_combined_annots.to_csv(ex_combined_anno_df_file, sep="\t")
					ex_combined_ann_counts = count_contacts(ex_combined_annots)
				ex_Edf = calculate_Escore("exons", i_cask_parent_dir, l, chrom_sizes_file, mybedfile, ex_combined_ann_counts, RNAag_totcounts_dict, g, genome_bin_size, rep_dpnII_f, gbin_dpnII_f, save=True)
			
			Rdf, model_fit_df = compute_GAM_residual(Edf, ex_Edf, max_ag_chrs)


		R_df_file = i_cask_parent_dir + "cenDNAxcenRNA_withalign_" + g + "gbin_ci2_and_unqv2." + l + "." + r + ".E.GAM.resid"
		R_model_stats_file = i_cask_parent_dir + "cenDNAxcenRNA_withalign_" + g + "gbin_ci2_and_unqv2." + l + "." + r + ".E.GAM.resid.fitstats"
		Rdf.to_csv(R_df_file, sep="\t")
		model_fit_df.to_csv(R_model_stats_file, sep="\t")

		print("Saved residuals to:", R_df_file)
		print("Saved model fit stats to:", R_model_stats_file)




#RUN as:
#parallel -j2 'python run_double_sided_calc_enrichment.py {1}' ::: WT_rep1_1 WT_rep1_2 WT_rep1_3 WT_rep2_1 WT_rep2_2