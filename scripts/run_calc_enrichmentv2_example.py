from calc_enrichment_v2 import *
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
import sys


############### Edit input files and parameters below ######################

#base and subdirectories for CASK output
cask_parent_dir_base = "/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/"
cask_subdir = "/cen_kmatch/k25.minusBKG/"
#Base directory of genome bin size 
gbin_parent_dir_base = "/scratch/groups/astraigh/kelsey/centrochar/charseq_reads/decon/WT_K562/"
#genome bin file subdir and name, this is what comes after sample or rep name dir
gbin_suffix = "/paired_dna_100kbgbins_1-22XM.bed"
#chromosome sizes bed file
chrom_sizes_file = "/oak/stanford/groups/astraigh/T2T/CURRENT/T2T-CHM13v2.0-hs1_jan2022/hs1.chromsizes.txt"
#repeat annotation bedfile with length
mybedfile = "/home/groups/astraigh/kelsey/chm13_genome/v2/censat.clean.chrclass.bed"

#dpnII sites per genome bin
gbin_dpnII_file = "/scratch/groups/astraigh/kelsey/centrochar/dpnii/hg38_1-22XM_100kb_genome_bins_dpnII_cut.txt"
#dpnII sites per repeat
rep_dpnII_file =  "/scratch/groups/astraigh/kelsey/centrochar/dpnii/cen_only_repeats_withClass_dpnII_cut.txt"


RNA_types=["exons"]
#RNA_types=["introns"]
levels=["class"]
#levels=["class"]

counts_df_dict={}
res_stats_df_dict={}

ag_side = 'dna'

genome_bin_size = 100000
g = '100kb1-22XM'

max_ag_chrs = 5

################################################################################

i = sys.argv[1]


def generate_Edf(cask_parent_dir, gbin_file, rep_bedfile, chromsizes_file, gbin_dpnII_f, rep_dpnII_f, ag_side, gbin_size, rep, rna_type, lev):
	print("Starting calculations for", rep, rna_type, lev)
	

	i_ag_file = cask_parent_dir + ag_side + "side.bed." + r + ".CASK_mergev2.txt"
	print("genome_bin_size",genome_bin_size)

	resolved_RNAcounts_df, res_stats_df  = load_ag_dat(rna_type, i_ag_file, lev, ag_side, gbin_file, g, combine_gbin=True, save=True)

	df = calculate_Escore(rna_type, cask_parent_dir, lev, chromsizes_file, rep_bedfile, resolved_RNAcounts_df, g, gbin_size, rep_dpnII_f, gbin_dpnII_f, save=True)

	return(df)



i_cask_parent_dir = cask_parent_dir_base + i + cask_subdir
i_gbin_file = gbin_parent_dir_base + i + gbin_suffix
for r in RNA_types:
	
	for l in levels:
		i_Edf_file = i_cask_parent_dir + "dna.fastq.gz.ANNOTATED.ci2_and_unqv2." + l + "." + r + "." + g + ".Escore"
		if os.path.exists(i_Edf_file):
			print(f"Loading exisiting Escore df from {i_Edf_file}")
			Edf = pd.read_csv(i_Edf_file, sep="\t", header=0, index_col=0) 
			Edf['ag_reptypes']=Edf['ag_reptypes'].apply(ast.literal_eval)
			Edf['ag_chrs']=Edf['ag_chrs'].apply(ast.literal_eval)
		else:
			Edf = generate_Edf(i_cask_parent_dir, i_gbin_file, mybedfile, chrom_sizes_file, gbin_dpnII_file, rep_dpnII_file, ag_side, genome_bin_size, i, r, l)


		if r == "exons":
			res_df, model_fit_df = compute_GAM_residual(Edf, Edf, max_ag_chrs)
		else:
			exon_Edf_file = i_cask_parent_dir + "dna.fastq.gz.ANNOTATED.ci2_and_unqv2." + l + ".exons." + g + ".Escore"
			if os.path.exists(exon_Edf_file):
				print(f"Loading exon Escore df from {exon_Edf_file}")
				ex_Edf = pd.read_csv(exon_Edf_file, sep="\t", header=0, index_col=0) 
				ex_Edf['ag_reptypes']=ex_Edf['ag_reptypes'].apply(ast.literal_eval)
				ex_Edf['ag_chrs']=ex_Edf['ag_chrs'].apply(ast.literal_eval)
				res_df, model_fit_df = compute_GAM_residual(Edf, ex_Edf, max_ag_chrs)
			else:
				print(f"File {exon_Edf_file} does not exist. Generating now.")

				ex_Edf = generate_Edf(i_cask_parent_dir, i_gbin_file, mybedfile, chrom_sizes_file, gbin_dpnII_file, rep_dpnII_file, ag_side, genome_bin_size, i, "exons", l)
				res_df, model_fit_df = compute_GAM_residual(Edf, ex_Edf, max_ag_chrs)		

		resid_outfile = i_cask_parent_dir + i + "_" + r + "_" + l + "_ci2_and_unq_" + g + ".E.GAM.resid"
		res_df['ag_reptypes']=res_df['ag_reptypes'].apply(lambda x: str(x))
		res_df['ag_chrs']=res_df['ag_chrs'].apply(lambda x: str(x))
		res_df.to_csv(resid_outfile, sep="\t")

		model_fit_stats_file = i_cask_parent_dir + i + "_" + r + "_" + l + "_ci2_and_unq_" + g + ".E.GAM_fit.stats"
		model_fit_df.to_csv(model_fit_stats_file, sep="\t")

		print("Finished calculating residuals. Saved fit stats as:", resid_outfile)
		print("Saved Rdf as:", resid_outfile)

#run as 
#parallel -j6 'python run_calc_enrichmentv2_example.py {1}' ::: WT_rep1_1 WT_rep1_2 WT_rep1_3 WT_rep2_1 WT_rep2_2