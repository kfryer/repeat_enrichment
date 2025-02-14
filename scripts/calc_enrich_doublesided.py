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
import os
import gc


def map_ag(parent_dir, level, df, side, ci2_col, unq_col):
	""" 
	Maps ambiv group ID to repeat and class types. 
	Inputs: 
		ci2 and unq pickle file from CASK
		read level ambiv group assignments as a pandas dataframe (from loag_ag_dat function)
	Outputs:
		read level ambiv group assignments including repeat and class level information 
	"""
	df = df.copy()
	ci2_pickle = parent_dir + side + ".fastq.gz.ANNOTATED.ci2.aglut.pickle"
	unq_pickle = parent_dir + side + ".fastq.gz.ANNOTATED.unq.aglut.pickle"

	if level == "class":
		ag_dict_key = 'ag_to_chr_LUD_class'		
	else:
		ag_dict_key = 'ag_to_chr_LUD'
		
	with open(ci2_pickle, 'rb') as ci2p:
		ci2_ag_dict=pickle.load(ci2p)
		ci2agdict = ci2_ag_dict[ag_dict_key]
	with open(unq_pickle, 'rb') as unqp:
		unq_ag_dict=pickle.load(unqp)
		unqagdict = unq_ag_dict[ag_dict_key]

	#map the names of the repeat types to the ambiv groups
	df["ag_reptypes_ci2"]=df[ci2_col].map(ci2agdict)
	df["ag_reptypes_unq"]=df[unq_col].map(unqagdict)

	return df


def resolve_reptypes(row):
	""" 
	Compares repeat assignments according to ci2 k-mers (k-mers that were present at least 2 times in the full repeat database) and unq k-mers (only found once in the entire repeat database) to resolve ambivalence groups where possible.
	Example:
		A read that is assigned to an ambiv group encompassing {repeat_1, repeat_2, repeat_3} based on k-mers that were found more than once in the entire repeat database, and assigned {repeat_1} based on k-mers that were only found once (in this case they were found in repeat_1). As long as there are at least 2 unq k-mers found in the read, the read will be resolved to {repeat_1} assignment.
	Inputs: 
		Row from repeat/class level mapped ag data from load_ag_dat function (each row is one read)
	Outputs:
		1. Resolved repeat type list
		2. Type of resolution 
	Resolution Types:
		ci2z_unqres				No ci2 k-mers found (ag ID = 0) or repeat types assigned based on ci2 k-mers were conflicting (ag ID = -1). Resolved repeat types is based only on the unq k-mer assignment

		unqz_ci2res				No unq k-mers found (ag ID = 0) or repeat types assigned based on unq k-mers were conflicting (ag ID = -1). Resolved repeat types is based only on the ci2 k-mer assignment

		low_kcount_nores		No assignment (ag ID = 0) or conlficting assignment (ag ID = -1) from either unq or ci2 k-mers and too few k-mers from the other set (unq or ci2). No resolution could occur and returns NA

		same 					Assignment based on unq and ci2 k-mers is identical. Returns the assignment from unq k-mer assignment

		issubset				Types assigned based on unq k-mers is a subset of those assigned based on ci2 k-mers. Returns the assignment from unq k-mer assignment (as long as there are at least 2 unq k-mers present in the read). Example: ci2 = {repeat_1, repeat_2, repeat_3}, unq = {repeat_1}, resolved = {repeat_1}
		
		ci2_wins				Types assigned based on unq k-mers is not a subset of those assigned based on ci2 k-mers. There are more ci2 k-mers in the read than unq k-mers, so resolved repeat types = ci2 repeat types. Example: ci2 = {repeat_1, repeat_2, repeat_3},ci2_kcount = 8; unq = {repeat_6}, unq_kcount=2; resolved = {repeat_1, repeat_2, repeat_3}.

		unq_wins				Same as ci2_wins, but there were more unq k-mers than ci2 k-mers. Outputs repeat types based on unq k-mers

		conflicting				There are equal number of k-mers from unq and ci2 and the repeat types are not a subset of eachother. No resolution can occur. Outputs NA

	Note:
		Resolutions only occur if there are at least 2 k-mers on the side being used to resolve repeat types. This is to prevent the presence of a single k-mer (which could be the result of a sequencing error) from dominating the read assignment.
		CASK ambiv group ID = 0 --> No repeat k-mers were found in the read
		CASK ambiv group ID = -1 --> The possible repeat types for the read were conflicting. Example: 4 kmers found in the read are found in the repeat database (k1, k2, k3, k4). k1 and k2 are found in {repeat_1} k3 is found in {repeat_1, repeat_3, repeat_5}, and k4 is found in {repeat_5}. Since k1 and k2 can only be found in repeat_1 and k4 can only be found in repeat_5, CASK can not discern which repeat it comes from. This could be due to sequencing errors, reads overlapping multiple repeats, or genome differences (SNPs).

	"""

	ci2_reptypes_list, unq_reptypes_list = row['ag_reptypes_ci2'], row['ag_reptypes_unq']
	ci2kcount, unqkcount = row["ci2kcount_inread"], row['unqkcount_inread']

	if ci2_reptypes_list in [['0'], ['-1']] and unq_reptypes_list not in [['0'], ['-1']]:
		if unqkcount > 1:
			return unq_reptypes_list, "ci2z_unqres"
		else:
			return ["NA"], "low_kcount_nores"

	elif unq_reptypes_list in [['0'], ['-1']] and ci2_reptypes_list not in [['0'], ['-1']]:
		if ci2kcount > 1:
			return ci2_reptypes_list, "unqz_ci2res"
		else:
			return ["NA"], "low_kcount_nores"
	elif ci2_reptypes_list in [['0'], ['-1']] and unq_reptypes_list in [['0'], ['-1']]:
		return ["NA"], "both_z"
	else:
		ci2_set, unq_set = set(ci2_reptypes_list), set(unq_reptypes_list)
		if ci2_set == unq_set:
			return unq_reptypes_list, "same"
		elif unq_set.issubset(ci2_set):
			if unqkcount > 1:
				return unq_reptypes_list, "subset"
			elif ci2kcount > 1:
				return ci2_reptypes_list, "low_kcount_nores"
			else:
				return ["NA"], "low_kcount_nores"
		else:
			if ci2kcount > unqkcount and ci2kcount > 1:
				return ci2_reptypes_list, "ci2_wins"
			elif unqkcount > ci2kcount and unqkcount > 1:
				return unq_reptypes_list, "unq_wins"
			else:
				return ["NA"], "conflicting"




def load_ag_dat(ag_file, level, ag_side, save=False):
	"""
	Loads read level ambiv group assignment file from CASK, maps the corresponding repeat types based on their IDs using the function map_ag, and resolves possible repeat types using the function resolve_reptypes.

	Inputs:
		1. Read level ambic group assignment file from CASK (i.e. dna.fastq.gz.ANNOTATED.ci2_and_unq.ag.txt)
		2. class or repeat type (level)
		3. ag_side - rna or dna side of charseq reads being assigned based on CASK. This is to load the appropriate LUT .pickle file
	Outputs:
		Dataframe with each row as an ag mapped, resolved read. This is also saved as a parquet file in the same directory as the CASK files (Input 1)
		Resolution stats file (how many reads were resolved based on each resolution type). Saved as .txt file in the same directory as above. 
	"""

	cask_parent_dir = "/".join(ag_file.split("/")[0:-1]) + "/"
	out_file = cask_parent_dir + ag_side + "_" + level + ".resolved_agmapped_reads.parquet"
	out_stats_file = cask_parent_dir + ag_side + "_" + "resolved_ag_stats_" + level + ".txt"

	if os.path.exists(out_file):
		print("Loading existing ag_dat:", out_file)
		#out_df = pd.read_csv(out_file, sep="\t", header=0, index_col=0)
		#out_df["ag_reptypes"] = out_df["ag_reptypes"].apply(ast.literal_eval)
		out_df = pd.read_parquet(out_file)
		out_df["ag_reptypes"] = out_df["ag_reptypes_str"].str.split(',')
		res_stats_df = pd.read_csv(out_stats_file, sep="\t", header=0, index_col=0)
	else:

		if level == "class":
			ci2_col = "class_aG_ci2"
			unq_col = "class_aG_unq"
		else:
			ci2_col = "reptype_ag_ci2"
			unq_col = "reptype_ag_unq"

		ag_df = pd.read_csv(ag_file, sep="\t", usecols=[0,1,2,3,4,5,6], names=['readid','reptype_ag_ci2','class_aG_ci2','ci2kcount_inread','reptype_ag_unq','class_aG_unq','unqkcount_inread'])
		level_df = ag_df[['readid',ci2_col,'ci2kcount_inread', unq_col,'unqkcount_inread']]
		mapped_ag_df = map_ag(cask_parent_dir, level, level_df, ag_side, ci2_col, unq_col)
		print("Finished mapping ags")

		# Resolve reptypes based on both ci2 and unq
		mapped_ag_df[['ag_reptypes','res_type']] = mapped_ag_df.apply(resolve_reptypes, axis=1, result_type='expand')
		res_stats_df = mapped_ag_df.groupby(['res_type']).size()
		print("Finished resolving reptypes")
		
		
		mapped_ag_df['ag_reptypes_str'] = (mapped_ag_df['ag_reptypes'].transform(lambda x: ",".join(map(str,x))))

		#filter conflicts and zeros
		mapped_ag_df=mapped_ag_df[mapped_ag_df['ag_reptypes_str']!= "NA"]

		res_stats_df["Total_resolved"]=len(mapped_ag_df)

		out_df = mapped_ag_df.drop(columns=['res_type'])

		
		
		if save:
			out_df = out_df.drop(columns=['ag_reptypes'])
			out_df.to_parquet(out_file)
			#out_df['ag_reptypes']=out_df['ag_reptypes'].apply(lambda x: str(x))
			print("Saved counts as: ", out_file)

			
			res_stats_df.to_csv(out_stats_file, sep="\t", header=True)
			print("Saved stats as: ", out_stats_file)

	return out_df, res_stats_df


def count_rep_RNAs(res_RNA_agmapped_df):
	"""
	Counts the number of RNA sides assigned to each ambiv group. This serves as the N_tot for repeat annotated RNAs to be used to normalize enrichment scores later.
	Inputs:
		RNA side output from load_ag_dat
	Output:
		Dictionary with the ambiv group  (a comma sep string of the repeat types possible for that amibv group) as keys and the total number of RNA side reads assigned to that ambiv group as the value
		{repeat_1,repeat_2: 4, repeat_5: 2, ...}
	"""
	RNAag_counts = res_RNA_agmapped_df.groupby('ag_reptypes_str').size()
	RNAag_counts_dict = RNAag_counts.to_dict()
	return RNAag_counts_dict


def ccmerge(dnaside_ag_df, rnaside_ag_df):
	"""
	Merge rna and dna side ag mapped dfs to get repeat annotations for both sides of the read.
	Inputs:
		dna side ag mapped, resolved dataframe from load_ag_dat
		rna side ag mapped, resolved dataframe from load_ag_dat
	Output:
		Single dataframe where each read is a row. Contains dna, then rna side repeat annotation information. All reads that have a repeat annotation on either side are included. Reads with repeat annotation only a single side have NA on the other side.

	"""
	#combine caskrna and cask dna
	dnaside_ag_df.set_index('readid', inplace=True)
	rnaside_ag_df.set_index('readid', inplace=True)
	
	#combine r & d, keep all reads from both sides
	caskmerge_df = dnaside_ag_df.merge(rnaside_ag_df, left_index=True, right_index=True, suffixes=("_dna","_rna"), how='outer')
	caskmerge_df = caskmerge_df.fillna("NA")
	#print("caskmerge_df:",caskmerge_df)

	return caskmerge_df


def combine_annotations(caskmerge_df, gbin_file, pairs_bedfile, RNA_type):
	"""
	Combines annotations to get RNA and DNA side information for every read.
	Read types:
		Both repeat:			Both RNA and DNA sides have repeat annotations (that are not NA)

		RNA repeat only:		RNA side has repeat annotation and DNA side did not, DNA side assigned to genomic bin based on gbin file

		DNA repeat only:		DNA side has repeat annotation and RNA side did not, RNA side assigned gene name (based on charseq alignment and paired file output from charseq pipeline)

		Neither repeats			RNA assigned gene name & DNA side assigned genomic bin

		Multiple annotations 	RNA side has both repeat and gene (alignment) assignment. Multimapping reads (map qual < 255) are assigned the repeat annotation. Unique alignments (map qual >= 255) are assigned alignment annotation. DNA side reads with multiple annotations were just assigned the repeat annotation since the alignment annotation is just a genomic bin.
	"""

	gbin_df = pd.read_csv(gbin_file, sep="\t", usecols=[3,9], names=['readid','gbin'])

	#Load rna alignment info from pairs file
	#{charseq_data_dir}/data/{sample_name}/pairs/gencondeV29_hg38/all/dna.bed.gz
	align_df = pd.read_csv(pairs_bedfile, sep="\t", usecols=[3,6,9,11,13,14,16], names=['readid','chrom_rna','mapq','ENST','gene_name','type', 'gene_ID'])
	
	#filter for exons 
	if RNA_type == "exons":
		#filter for RNA sides that are from exons (have ENSTs)
		filt_align_df = align_df[align_df['ENST'].str.startswith('ENST')]
	else:
		#if not filtering for exons, use ENSG to map genes
		filt_align_df = align_df[align_df['ENST'].str.startswith('ENSG')]

	filt_align_df.set_index('readid', inplace=True)
	filt_align_df = filt_align_df.copy()
	filt_align_df.replace('*',"NA", inplace=True)
	
	#merge pairs file info with cask annotations df
	cask_pairs_merge_full = caskmerge_df.merge(filt_align_df, left_index=True, right_index=True, how='outer')
	#print("cask_pairs_merge_full:", cask_pairs_merge_full)
	#print("cask_pairs_merge_full cols:", cask_pairs_merge_full.columns.tolist())

	del align_df
	del filt_align_df
	del caskmerge_df
	gc.collect()

	#merge cask, pairs info with genomic bins to get all layers of annotation
	cask_pairs_merge_full = cask_pairs_merge_full.reset_index()
	cask_pairs_merge_full_g = cask_pairs_merge_full.merge(gbin_df, on='readid', how='left')
	#print("cask_pairs_merge_full_g cols:", cask_pairs_merge_full_g.columns.tolist())
	#filter for relevant columns
	cask_pairs_merge_cut = cask_pairs_merge_full_g[['readid','ag_reptypes_dna', 'ag_reptypes_str_dna', 'ag_reptypes_rna', 'ag_reptypes_str_rna', 'chrom_rna', 'mapq', 'ENST', 'gene_name', 'type', 'gene_ID', 'gbin']]
	
	
	cask_pairs_merge_df = cask_pairs_merge_cut.fillna("NA").copy()

	del cask_pairs_merge_full
	del cask_pairs_merge_cut
	del gbin_df
	gc.collect()

	#convert the genome bin annotaiton to a list (it's just a single item in the list) so it can be processed like the other annotation types
	cask_pairs_merge_df['gbin_l'] = cask_pairs_merge_df['gbin'].str.split(',')

	def extract_chrs(s):
		"""
		Extracts the chromosome or chromosomes from the genome bin name or the list of repeat type name(s) 
		"""
		if s != "NA":
			if not s[0].startswith('chr'):
				return list(set(['chr'+x.split('_')[1] for x in s]))
			else:
				return list(set([x.split('_')[0] for x in s]))
		else:
			return ["NA"]



	#filter out reads with no DNA side annotation. These are reads that were not given a repeat annotation and did not align to the genomic regions included in the gbin file (ie they aligned to unassembled scaffolds or chromosome Y). This is dependent on which regions were included when gbin file was made.
	DNAfilt_merge = cask_pairs_merge_df[~((cask_pairs_merge_df['ag_reptypes_str_dna'] == "NA") & (cask_pairs_merge_df['gbin'] == "NA"))].copy()

	#create final DNA annotation (DNA_annot and DNA_chrs) for reads that have repeat assignment on the DNA side
	repD = DNAfilt_merge[DNAfilt_merge['ag_reptypes_str_dna'] != "NA"]
	repDNA = repD.copy()
	repDNA['DNA_annot'] = repDNA['ag_reptypes_str_dna']
	repDNA['DNA_chrs'] = repDNA['ag_reptypes_dna'].apply(extract_chrs)

	#create final DNA annotation (DNA_annot and DNA_chrs) for reads that have no repeat assignment on the DNA side (they get genome bin assignments)
	norepD = DNAfilt_merge[DNAfilt_merge['ag_reptypes_str_dna'] == "NA"]
	norepDNA = norepD.copy()
	norepDNA['DNA_annot'] = norepDNA['gbin']
	norepDNA['DNA_chrs'] = norepDNA['gbin_l'].apply(extract_chrs) 

	#Concat repeat and non-repeat DNA side dataframes
	resDNA = pd.concat([repDNA,norepDNA], ignore_index=True)

	#clean up
	del DNAfilt_merge
	del cask_pairs_merge_df
	del repD
	del repDNA
	del norepD
	del norepDNA
	gc.collect()

	#resolve RNA side annotations

	#filter out reads with no repeat annotation and no gene_name 
	resDNA_RNAfilt = resDNA[(resDNA['ag_reptypes_str_rna'] != "NA") | (resDNA['gene_name'] != "NA")]

	#create final RNA side annotations for reads with only repeat annotation on RNA side (no gene_name assigned)
	reponly_RNA = resDNA_RNAfilt[(resDNA_RNAfilt['ag_reptypes_str_rna'] != "NA") & (resDNA_RNAfilt['gene_name'] == "NA")].copy()
	reponly_RNA['RNA_annot'] = reponly_RNA['ag_reptypes_str_rna']
	reponly_RNA['RNA_annot_ID'] = reponly_RNA['ag_reptypes_str_rna']
	reponly_RNA['RNA_chrs'] = reponly_RNA['ag_reptypes_rna'].apply(extract_chrs)
	reponly_RNA['RNA_type'] = 'rep'

	
	#create final RNA side annotations for reads with gene_name annotations that don't have repeat annotations
	alignonly_RNA = resDNA_RNAfilt[(resDNA_RNAfilt['ag_reptypes_str_rna'] == "NA") & (resDNA_RNAfilt['gene_name'] != "NA")].copy()
	alignonly_RNA['RNA_annot'] = alignonly_RNA['gene_name']
	alignonly_RNA['RNA_annot_ID'] = alignonly_RNA['gene_ID']
	alignonly_RNA['RNA_chrs'] = alignonly_RNA['chrom_rna'].str.split(',')
	alignonly_RNA['RNA_type'] = alignonly_RNA['type']
	
	
	
	#Filter for reads that have multiple annotations (repeat and alignment). 
	RNAags = resDNA_RNAfilt[resDNA_RNAfilt['ag_reptypes_str_rna'] != "NA"].copy()
	multi_annot_RNAs = RNAags[RNAags['gene_name'] != "NA"]

	#Set gene_name as RNA annotation for reads that have repeat and gene annotations, but are unique alignments (mapq = 255)
	q255_multiRNA = multi_annot_RNAs[multi_annot_RNAs['mapq'] >= 255].copy()
	q255_multiRNA['RNA_annot'] = q255_multiRNA['gene_name']
	q255_multiRNA['RNA_annot_ID'] = q255_multiRNA['gene_ID']
	q255_multiRNA['RNA_chrs'] = q255_multiRNA['chrom_rna'].str.split(',')
	q255_multiRNA['RNA_type'] = q255_multiRNA['type']
	#print("q255_multiRNA", q255_multiRNA.shape)
	
	#Set repeat annotation as RNA annotation for reads that have repeat and gene annotations, and are multimappers (mapq < 255)
	qlo_multiRNA = multi_annot_RNAs[multi_annot_RNAs['mapq'] < 255].copy()
	qlo_multiRNA = qlo_multiRNA[qlo_multiRNA['ag_reptypes_str_rna'] != "NA"]
	qlo_multiRNA['RNA_annot'] = qlo_multiRNA['ag_reptypes_str_rna']
	qlo_multiRNA['RNA_annot_ID'] = qlo_multiRNA['ag_reptypes_str_rna']
	qlo_multiRNA['RNA_chrs'] = qlo_multiRNA['ag_reptypes_rna'].apply(extract_chrs)
	qlo_multiRNA['RNA_type'] = 'rep'
	

	#combine resolved RNA annotations dataframes
	RNA_dfs = [reponly_RNA, alignonly_RNA, q255_multiRNA, qlo_multiRNA]
	DNAfilt_RNAres = pd.concat(RNA_dfs)	

	#clean up
	del reponly_RNA
	del alignonly_RNA
	del q255_multiRNA
	del qlo_multiRNA
	del multi_annot_RNAs
	del RNAags
	gc.collect()


	#convert list cols to strings
	DNAfilt_RNAres['DNA_chr']=(DNAfilt_RNAres['DNA_chrs'].transform(lambda x: ",".join(map(str,x))))
	DNAfilt_RNAres['RNA_chr'] = (DNAfilt_RNAres['RNA_chrs'].transform(lambda x: ",".join(map(str,x))))

	#filter for relevant columns. Excluding list columns because they mess with dataframe operations
	DNAfilt_RNAres = DNAfilt_RNAres[['readid','DNA_annot','DNA_chr','RNA_annot','RNA_annot_ID','RNA_chr','RNA_type']]

	print("Finished combining annotations")
	return DNAfilt_RNAres


def count_contacts(comb_annot_df):
	"""
	Count the number of reads with each combo of RNA and DNA annotations. This is to count the number of times a particular RNA-DNA contact was observed in the dataset. Set this count at N_in (for use in enrichment calculation)
	"""

	sum_cols = ['DNA_annot','DNA_chr','RNA_annot','RNA_annot_ID','RNA_chr','RNA_type']

	contact_counts_df = comb_annot_df.groupby(sum_cols, as_index=False).size()
	contact_counts_df.rename(columns={"size":"N_in"}, inplace=True)

	return contact_counts_df



def calculate_Escore(RNA_type, parent_dir, level, chrom_sizes_filepath, bedfile, resmapped_ag_df, RNAag_totcounts_dict, g, gbin_size, rep_dpnII_file, gbin_dpnII_file,save=False):
	"""
	Calculate an enrichment score for each RNA-DNA contact observed. The number of times a contact was captured (N_in) is normalized by the total number of times the RNA was found in the entire dataset (N_tot), the length of the DNA contact domain (L_in) which is either the length of all the repats in the ambiv group or the size of the genomic bin, and the number of DpnII sites found in the contact domain(s).
	
	N_out = N_tot - N_in
	L_out = L_tot - L_in
	((N_in/D_in)/N_out)((L_out/L_in))
	"""

	resolved_RNAcounts_df = resmapped_ag_df.copy()
	
	#import the sizes of each chromosome & calculate the total length of the genome in basepairs
	chrom_sizes_df = pd.read_csv(chrom_sizes_filepath, sep="\t", usecols=[0,1], names=["chrom","size"])
	Ltotal=chrom_sizes_df["size"].sum()

		
	#get repeat type and classinfo as well as length from bedfile
	bed_df = pd.read_csv(bedfile, sep="\t", usecols=[0,1,2,3,4,6], names=["chrom","start","stop","reptype","length","class"])
	class_length_s=bed_df.groupby(["class"])["length"].sum()
	class_length_df=class_length_s.to_frame().reset_index()
	
	#get total count for each non repeat derived RNA from CASK output .countallRNAs.txt (N_tot) & convert to dictionary
	tot_RNA_counts_file = parent_dir + "dnaside.bed." + RNA_type +".countallRNAs.txt"
	tot_RNA_counts_df = pd.read_csv(tot_RNA_counts_file, sep="\t", usecols=[0,1], names=["N_tot","gene_ID"])
	tot_RNA_counts_df=tot_RNA_counts_df.set_index("gene_ID")
	tot_RNAcounts_dict = tot_RNA_counts_df.to_dict(orient="dict")["N_tot"]

	#combine non-repetitive RNA counts dict with repeat RNA counts dict (output from count_rep_RNAs function)
	all_tot_RNAcounts_dict = tot_RNAcounts_dict | RNAag_totcounts_dict

	#get dpnII counts for each bin or reptype
	rep_dpnII_df = pd.read_csv(rep_dpnII_file, sep=" ", usecols = [0,1,2], names=['reptype','dpnII_count','class'])
	class_dpnII_s = rep_dpnII_df.groupby(["class"])["dpnII_count"].sum()
	class_dpnII_df = class_dpnII_s.to_frame().reset_index()
	gbin_dpnII_df = pd.read_csv(gbin_dpnII_file, sep="\t", usecols = [0,1], names=['gbin','dpnII_count'])

	
	#add Ntot from tot_RNAcounts_dict to resolved counts dataframe (from count_contacts function)
	resolved_RNAcounts_df["N_tot"]=resolved_RNAcounts_df["RNA_annot_ID"].map(all_tot_RNAcounts_dict)
	#print(resolved_RNAcounts_df[['DNA_annot','RNA_annot','RNA_annot_ID','N_tot','N_in']])

	#collapse duplicate gene_names
	all_columns = resolved_RNAcounts_df.columns
	groupby_columns = ['DNA_annot', 'RNA_annot']
	sum_columns = ['N_in','N_tot']

	# Create the aggregation dictionary
	agg_dict = {col: 'first' for col in all_columns if col not in groupby_columns + sum_columns}
	for col in sum_columns:
		agg_dict[col] = 'sum'

	# Aggregate duplicate gene_names
	resolved_RNAcounts_df = resolved_RNAcounts_df.groupby(groupby_columns).agg(agg_dict).reset_index()

	#remove rows where N_in = N_tot (these will result in invalid value for log)
	resolved_RNAcounts_df = resolved_RNAcounts_df[resolved_RNAcounts_df["N_tot"]!=resolved_RNAcounts_df["N_in"]]

	#split DNA annot into list to be used for domain length calculation
	resolved_RNAcounts_df['DNA_annot_l'] = resolved_RNAcounts_df['DNA_annot'].str.split(',')

	#count the number of dpnII sites and then total length for each. Example: Ambiv group that contains repeat_1, repeat_2, repeat_3 will have a dpnII count equal to the sum of the number of dpnII sites found in repeat_1, repeat_2, and repeat_3. Total length is calculated in the same manner. 
	#Note: 1 is added to the dpnII count for each domain so domains with no DpnII sites are not 0.
	def calc_dpnII_count_rep(row):
		if row['DNA_annot'].startswith('chr'):
			curr_gbin = row['DNA_annot']
			tot_dpnII_count = sum(gbin_dpnII_df[gbin_dpnII_df['gbin'] == curr_gbin]['dpnII_count']) + 1
			return tot_dpnII_count
		else:
			reptypes = row['DNA_annot_l']
			tot_dpnII_count = sum(rep_dpnII_df[rep_dpnII_df['reptype'].isin(reptypes)]['dpnII_count']) + 1
			return tot_dpnII_count

	def calc_dpnII_count_class(row):
		if row['DNA_annot'].startswith('chr'):
			curr_gbin = row['DNA_annot']
			tot_dpnII_count = sum(gbin_dpnII_df[gbin_dpnII_df['gbin'] == curr_gbin]['dpnII_count']) + 1
			return tot_dpnII_count
		else:
			reptypes = row['DNA_annot_l']
			tot_dpnII_count = sum(class_dpnII_df[class_dpnII_df['class'].isin(reptypes)]['dpnII_count']) + 1
			return tot_dpnII_count
	
	#calculate Lin
	def calculate_total_length(row):
		if row['DNA_annot'].startswith('chr'):
			return gbin_size
		else:
			reptypes = row['DNA_annot_l']
			total_length = sum(bed_df[bed_df['reptype'].isin(reptypes)]['length'])
			return total_length
	
	def calculate_total_length_class(row):
		if row['DNA_annot'].startswith('chr'):
			return gbin_size
		else:
			reptypes = row['DNA_annot_l']
			total_length = sum(class_length_df[class_length_df["class"].isin(reptypes)]["length"])
			return total_length
	
	if level == "reptype":
		resolved_RNAcounts_df['L_in'] = resolved_RNAcounts_df.apply(calculate_total_length, axis=1)
		resolved_RNAcounts_df['dpnII'] = resolved_RNAcounts_df.apply(calc_dpnII_count_rep, axis=1)
	if level == "class":
		resolved_RNAcounts_df['L_in'] = resolved_RNAcounts_df.apply(calculate_total_length_class, axis=1)
		resolved_RNAcounts_df['dpnII'] = resolved_RNAcounts_df.apply(calc_dpnII_count_class, axis=1)
	
	#calculate L_out by subtracting L_in from L_tot
	resolved_RNAcounts_df['L_out'] = Ltotal - resolved_RNAcounts_df['L_in']

	#dpnII norm
	resolved_RNAcounts_df["N_in_dnorm"] = resolved_RNAcounts_df["N_in"] / resolved_RNAcounts_df['dpnII']
	
	#calculate raw E_score
	resolved_RNAcounts_df["raw_E"] = resolved_RNAcounts_df["N_in_dnorm"] / (resolved_RNAcounts_df["N_tot"] - resolved_RNAcounts_df["N_in"])
	
	#calculate E score
	resolved_RNAcounts_df["Escore"] = resolved_RNAcounts_df["raw_E"] * (resolved_RNAcounts_df["L_out"]/resolved_RNAcounts_df["L_in"])

	#print negative and inf and Na Tscores
	neg_E = resolved_RNAcounts_df[resolved_RNAcounts_df["Escore"] < 0]
	if not neg_E.empty:
		print("Negative Escores:", neg_E[['DNA_annot','RNA_annot','L_in','L_out','N_in','N_tot','raw_E','Escore']])

	rows_with_inf = resolved_RNAcounts_df[np.isinf(resolved_RNAcounts_df["Escore"])]
	if not rows_with_inf.empty:
		print("Inf Escores:", rows_with_inf[['DNA_annot','RNA_annot','L_in','L_out','N_in','N_tot','raw_E','Escore']])

	rows_with_na = resolved_RNAcounts_df[resolved_RNAcounts_df["Escore"].isna()]
	if not rows_with_na.empty:
		print("NA Escores:", rows_with_na[['DNA_annot','RNA_annot','L_in','L_out','N_in','N_tot','raw_E','Escore']])


	#calculate log(Escore)
	resolved_RNAcounts_df["logEscore"] = np.log2(resolved_RNAcounts_df["Escore"])

	#add logNtot for plotting later on
	resolved_RNAcounts_df['logNtot'] = np.log(resolved_RNAcounts_df["N_tot"])

	#rearrange order of columns to make it pretty & output to a file
	resolved_RNAcounts_df = resolved_RNAcounts_df[['DNA_annot', 'DNA_chr', 'RNA_annot', 'RNA_chr','RNA_type','N_tot','logNtot','N_in','L_out','L_in','raw_E','Escore','logEscore']]
	
	if save:
		ofile_name = parent_dir + "cenDNAxcenRNA_withalign_" + g + "gbin_ci2_and_unqv2." + level + "." + RNA_type + ".Escore"
		df_out=resolved_RNAcounts_df.copy()
		df_out.to_csv(ofile_name, sep="\t")
		print("Saved Edf as",ofile_name)

	return resolved_RNAcounts_df



def model_and_predict(pc_df,all_df):
	"""
	Model background contact enrichment scores based on the abundance of the RNA (Ntot) and the length of the domain (Lin) for background contact data from pc_df. Use that model to predict enrichment scores for all_df based on their Ntot and Lin. 
	Output:
		all_df with predicted enrichment scores column
		r-squared value for all_df data fit to model of pc_df
	"""
	adf = all_df.copy()

	#X = pc_df['logNtot'].values
	X = pc_df[['logNtot', 'L_in']].values
	y = pc_df['logEscore'].values
	X_test = all_df['logNtot'].values
	

	gam = LinearGAM().fit(X,y)
	rsq = gam.statistics_['pseudo_r2']['explained_deviance']
	

	adf['predicted_logEscore'] = gam.predict(adf[['logNtot','L_in']])
	return adf, rsq



def gettransRNAs(chrom_list, exons_df):
	"""
	Filter for RNAs found contacting regions outside of their own chromosome - trans contacting RNAs. These will represent background contacts. 
	"""

	exons_dfc = exons_df.copy()
	exons_dfc['RNA_chrs'] = exons_dfc['RNA_chr'].str.split(',')

	notchrom_df = exons_dfc[exons_dfc['RNA_chrs'].apply(lambda x: not any(c in chrom_list for c in x))]

	return(notchrom_df)



def compute_GAM_residual(in_df, in_ex_df, max_ag_chr_val):
	"""
	Calculate a background model for each chromosome or group of chromosomes (for ambiv groups) and predict enrichment scores for all contacts on that chromosome. 
	Inputs:
		in_df 			df from calculate_Escore function that contains the enrichment scores for each RNA-DNA contact

		in_ex_df 		df from calculate_Escore function that contains the enrichment scores for each RNA-DNA contact where RNA alignment annotations are all exonic

		max_ag_chr_val	the maximum number of chromosomes a single ambivalence group can encompass to be considered in the model (this is to exclude ambiv groups that occupy almost all the chromosomes and thus the background model would be based on a tiny subset of the data)

	Outputs:
		1. Dataframe containing all the original info from in_df with predicted enrichment scores and residual scores
		2. Dataframe containing r-squared values and the number of data points included for each background model (for each chromosome)
	"""
	
	tmp_df=in_df.copy()
	#exclude contacts where N_in is equal to N_tot because these have inf enrichment scores
	df_out = tmp_df.drop(tmp_df[tmp_df['N_tot'] == tmp_df['N_in']].index)
	resid_dict={}

	
	df_out['DNA_chrs'] = df_out['DNA_chr'].str.split(',')

	#filter out contacts with ambiv groups that encompass more than max_ag_chr_val number of chromosomes
	df_out = df_out[df_out['DNA_chrs'].apply(lambda x: len(x) < max_ag_chr_val)]

	tmp_ex_df=in_ex_df.copy()
	ex_df=tmp_ex_df.drop(tmp_ex_df[tmp_ex_df['N_tot'] == tmp_ex_df['N_in']].index)
	ex_df['DNA_chrs'] = ex_df['DNA_chr'].str.split(',')
	ex_df = ex_df[ex_df['DNA_chrs'].apply(lambda x: len(x) < max_ag_chr_val)]

	model_stats = []
	#Calculate background model for each chromosome or group of ch romosomes that belong to an ambiv group
	#For each chromosome or group of chromosomes, i, select all protein coding exon derived RNAs contacting any chromosome except i. This represents the background data
	#use that background data to predict enrichment scores for all RNAs contacting i (including repeat derived RNAs)
	for i in df_out["DNA_chr"].unique():
		if i not in resid_dict.keys():
			#select all RNAS contacting the chromosomes in i & convert i to chrom_list
			i_df = df_out[df_out['DNA_chr']==i]
			chrom_list = i_df['DNA_chrs'].iloc[0]

			#select all exon derived RNAs contacting the chromosomes in i
			i_ex_df = ex_df[ex_df['DNA_chr']==i]
			
			#for ambiv groups that include multiple chromosomes, include only RNAs contacting chromosomes not in that ambiv group
			#filter for contacts including protein coding RNAs that have at least 2 reads 
			if len(chrom_list) > 1:
				coi = ex_df[ex_df['DNA_chrs'].apply(lambda x: any(c in chrom_list for c in x))]
				ex_pc_df = coi.query('(RNA_type == "protein_coding") & (N_in>1)')
			else:
				ex_pc_df = i_ex_df.query('(RNA_type == "protein_coding") & (N_in>1)')
			
			#return empty dataframe if the number of contacts meeting the above criteria is less than 10. This is to prevent trying to fit a model with a tiny subset of the data
			#filter for trans contacting RNAs using gettransRNAs function -- this is the background data used for model
			if len(ex_pc_df) < 10:
				print(i, "not enough data")
				resid_dict[i] = pd.DataFrame()
			else:
				ex_pc_offchr_df = gettransRNAs(chrom_list,ex_pc_df)

				if len(ex_pc_offchr_df) < 10:
					print(str(i) + " Does not contain enough samples. Using all protein coding RNAs as background")
					bkg_data = gettransRNAs([],ex_pc_df)
				else:
					bkg_data = ex_pc_offchr_df.copy()

				#Send background dataset and all contacts for given chromosome to model_and_predict function
				result, rsq = model_and_predict(bkg_data,i_df)
				result_df = result.copy()
				#calculate Residual (R) score by subtracting the observed Enrichment score from the predicted 
				result_df['R'] = result_df['logEscore'] - result_df['predicted_logEscore']
				#add all the data for the given chromosome to a dictionary
				resid_dict[i]=result_df
				#store fit stats and the number of samples used to generate the background model
				model_stats.append({'chr_group':i, 'rsquared':rsq, 'nSamples':ex_pc_offchr_df.size})
	
	#convert model stats list to dataframe including only chromosomes that had enough data to model
	model_stats_df = pd.DataFrame(model_stats)
	resid_dict_filt = {}
	for k in resid_dict.keys():
		df = resid_dict[k]
		if not df.empty:
			resid_dict_filt[k]=df
	

	#concatenate all contact data for all chromosomes 
	big_df = pd.concat(resid_dict_filt.values())



	return big_df, model_stats_df