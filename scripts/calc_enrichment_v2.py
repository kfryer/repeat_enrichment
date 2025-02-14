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
import sys


def map_ag(parent_dir, level, df, side, ci2_col, unq_col):
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



def filt_gbin(cask_df, paired_gbin_file, RNA_type):
	#cask_df = pd.read_csv(cask_merge_file, sep="\t", usecols=[0,1,2,3,4,5,6,7], names=['readid', 'ENSG', 'reptype_ag_ci2', 'class_aG_ci2', 'ci2kcount', 'reptype_ag_unq', 'class_aG_unq', 'unqkcount'])


	cask_readids = cask_df['readid'].tolist()

	#print("cask_reads:",len(cask_df))

	pairs_gbin_df = pd.read_csv(paired_gbin_file, sep="\t", usecols=[3,6,9,11,13,14,16,22], names=['readid','chrom','mapq','ENST','gene_name','type','gene_ID','gbin'])

	if RNA_type == "exons":
		#filter for RNA sides that are from exons (have ENSTs)
		pairs_gbin_df = pairs_gbin_df[pairs_gbin_df['ENST'].str.startswith('ENST')]
	else:
		pairs_gbin_df = pairs_gbin_df[pairs_gbin_df['ENST'].str.startswith('ENSG')]
	
	#print("post exon/intron filter:", pairs_gbin_df.shape)

	filt_pairs_gbin_df = pairs_gbin_df[~pairs_gbin_df['readid'].isin(cask_readids)]

	#print("post cask filter:",filt_pairs_gbin_df.shape)

	#count RNAs in filtered df
	sum_cols=['chrom','gene_name','type','gene_ID','gbin']
	filt_gbin_counts_df = filt_pairs_gbin_df.groupby(sum_cols, as_index=False).size()
	filt_gbin_counts_df.rename(columns={"size":"N_in","gbin":"ag_reptypes_str"}, inplace=True)
	#print("filt_gbin_counts_df", filt_gbin_counts_df)
	out_df = filt_gbin_counts_df.copy()
	#filt_gbin_df = filt_gbin_paired(mapped_ag_df, paired_gbin_file, RNA_type)
	#filt_gbin_df.rename(columns={"gbin":"ag_reptypes_str"}, inplace=True)
	#cask_df is resolved_RNAcounts_df from loag_ag_dat
	

	print("Finished combining gbin counts")
	
	return out_df
	

def load_ag_dat(RNA_type, ag_file, level, ag_side, paired_gbin_file, g, combine_gbin=True, save=False):
	#load counts for each RNA at each repeat type or class
	#RNAcounts_ag_file = parent_dir + "dnaside.bed." + RNA_type +".CASK_mergev2."+ level +"count_ci2_and_unqv2.txt"
	cask_parent_dir = "/".join(ag_file.split("/")[0:-1]) + "/"

	out_stats_file = cask_parent_dir + "resolved_ag_stats_" + RNA_type + "_" + level + ".txt"

	if combine_gbin:
		out_file = cask_parent_dir + "dnaside.bed." + RNA_type + "_" + level + ".CASK_mergev2_" + g + "gbin_counts.txt"
	else:
		out_file = cask_parent_dir + "dnaside.bed." + RNA_type + "_" + level + ".CASK_mergev2_counts.txt"

	if os.path.exists(out_file):
		print("Loading existing ag_dat:", out_file)
		#out_df = pd.read_csv(out_file, sep="\t", header=0, index_col=0)
		#out_df["ag_reptypes"] = out_df["ag_reptypes"].apply(ast.literal_eval)
		out_df = pd.read_csv(out_file, sep="\t", header=0, index_col=0)
		out_df["ag_reptypes"] = out_df["ag_reptypes_str"].str.split(',')
		res_stats_df = pd.read_csv(out_stats_file, sep="\t", header=0, index_col=0)
	else:
	
		if level == "class":
			ci2_col = "class_aG_ci2"
			unq_col = "class_aG_unq"
		else:
			ci2_col = "reptype_ag_ci2"
			unq_col = "reptype_ag_unq"
		
		#ag_file = parent_dir + "dnaside.bed." + RNA_type + ".CASK_mergev2.txt"
		ag_df = pd.read_csv(ag_file, sep="\t", usecols=[0,1,2,3,4,5,6,7,8,9,10], names=['readid','chrom','gene_name','type','gene_ID','reptype_ag_ci2','class_aG_ci2','ci2kcount_inread','reptype_ag_unq','class_aG_unq','unqkcount_inread'])
		ag_level_df = ag_df[['readid','chrom','gene_name','type','gene_ID', ci2_col,'ci2kcount_inread', unq_col,'unqkcount_inread']]

		mapped_ag_df = map_ag(cask_parent_dir, level, ag_level_df, ag_side, ci2_col, unq_col)
		

		print("Finished mapping ags")
		
		# Resolve reptypes based on both ci2 and unq
		mapped_ag_df[['ag_reptypes','res_type']] = mapped_ag_df.apply(resolve_reptypes, axis=1, result_type='expand')

		res_stats_df = mapped_ag_df.groupby(['res_type']).size()
		

		print("Finished resolving reptypes")
		
		
		mapped_ag_df['ag_reptypes_str'] = (mapped_ag_df['ag_reptypes'].transform(lambda x: ",".join(map(str,x))))

		#filter conflicts and zeros
		mapped_ag_df=mapped_ag_df[mapped_ag_df['ag_reptypes_str']!= "NA"]

		#combine counts for resolved ags
		sum_cols=['ag_reptypes_str','chrom','gene_name','type','gene_ID']
		resolved_RNAcounts_df = mapped_ag_df.groupby(sum_cols, as_index=False).size()
		
		resolved_RNAcounts_df.rename(columns={"size":"N_in"}, inplace=True)
		res_stats_df["Total_resolved"]=resolved_RNAcounts_df["N_in"].sum()
		#print(resolved_RNAcounts_df)

		if combine_gbin:
			#filter cask readids from gbin dat & count rnas
			filt_gbin_df = filt_gbin(mapped_ag_df, paired_gbin_file, RNA_type)
			#print(filt_gbin_df)
			
			#combine with ag counts
			out_df = pd.concat([resolved_RNAcounts_df,filt_gbin_df], ignore_index=True)
			#g = paired_gbin_file.split("/")[-1].split("_")[1].split("gbin")[0]
			#out_file = cask_parent_dir + "dnaside.bed." + RNA_type + "_" + level + ".CASK_mergev2_" + g + "gbin_counts.txt"
		else:
			out_df = resolved_RNAcounts_df.copy()
			#out_file = cask_parent_dir + "dnaside.bed." + RNA_type + "_" + level + ".CASK_mergev2_counts.txt"


		print("Finished counting RNAs for ags")
		
		if save:
			out_df.to_csv(out_file, sep="\t", header=True)
			print("Saved counts as: ",out_file)

			
			res_stats_df.to_csv(out_stats_file, sep="\t", header=True)
			print("Saved stats as: ",out_stats_file)

	return out_df, res_stats_df

def calculate_Escore(RNA_type, parent_dir, level, chrom_sizes_filepath, bedfile, resmapped_ag_df, g, gbin_size, rep_dpnII_file, gbin_dpnII_file, save=False):

	resolved_RNAcounts_df = resmapped_ag_df.copy()
	
	chrom_sizes_df = pd.read_csv(chrom_sizes_filepath, sep="\t", usecols=[0,1], names=["chrom","size"])
	Ltotal=chrom_sizes_df["size"].sum()

		
	#get repeat type and classinfo as well as length from bedfile
	bed_df = pd.read_csv(bedfile, sep="\t", usecols=[0,1,2,3,4,6], names=["chrom","start","stop","reptype","length","class"])
	class_length_s=bed_df.groupby(["class"])["length"].sum()
	class_length_df=class_length_s.to_frame().reset_index()
	
	#get total count for each RNA (Ntot)
	tot_RNA_counts_file = parent_dir + "dnaside.bed." + RNA_type +".countallRNAs.txt"
	tot_RNA_counts_df = pd.read_csv(tot_RNA_counts_file, sep="\t", usecols=[0,1], names=["N_tot","gene_ID"])
	tot_RNA_counts_df=tot_RNA_counts_df.set_index("gene_ID")
	tot_RNAcounts_dict = tot_RNA_counts_df.to_dict(orient="dict")["N_tot"]

	#get dpnII counts for each bin or reptype
	rep_dpnII_df = pd.read_csv(rep_dpnII_file, sep=" ", usecols = [0,1,2], names=['reptype','dpnII_count','class'])
	class_dpnII_s = rep_dpnII_df.groupby(["class"])["dpnII_count"].sum()
	class_dpnII_df = class_dpnII_s.to_frame().reset_index()

	gbin_dpnII_df = pd.read_csv(gbin_dpnII_file, sep="\t", usecols = [0,1], names=['gbin','dpnII_count'])



	#load counts for each RNA at each repeat type or class
	#RNAcounts_ag_file = parent_dir + "dnaside.bed." + RNA_type +".CASK_mergev2."+ level +"count_ci2_and_unqv2.txt"
	#RNAcounts_ag_df = pd.read_csv(RNAcounts_ag_file, sep="\t", usecols=[0,1,2,3,4,5,6], names=["N_in", "chrom", "gene_name", "type", "gene_ID", "ag_ci2", "ag_unq"])

	#calculate N_in
	#ag_df['ag_reptypes_str'] = (ag_df['ag_reptypes'].transform(lambda x: ",".join(map(str,x))))
	#sum_cols=['ag_reptypes_str','gene_ID', 'type', 'gene_name', 'chrom', 'N_tot']
	#RNAcounts_ag_df = ag_df.groupby(sum_cols, as_index=False)['N_in'].sum()
	#RNAcounts_ag_df['ag_reptypes'] = RNAcounts_ag_df['ag_reptypes_str'].str.split(',')

	#add Ntot from tot_RNAcounts_dict
	resolved_RNAcounts_df["N_tot"]=resolved_RNAcounts_df["gene_ID"].map(tot_RNAcounts_dict)

	#print("gene_name", resolved_RNAcounts_df[resolved_RNAcounts_df["gene_name"] == 'AC007495.1'])

	#print("after dict map", resolved_RNAcounts_df[resolved_RNAcounts_df["N_tot"] == 0])
	#collapse duplicate gene_names
	all_columns = resolved_RNAcounts_df.columns
	groupby_columns = ['ag_reptypes_str', 'gene_name']
	sum_columns = ['N_in','N_tot']

	# Create the aggregation dictionary
	agg_dict = {col: 'first' for col in all_columns if col not in groupby_columns + sum_columns}
	for col in sum_columns:
		agg_dict[col] = 'sum'

	# Aggregation
	resolved_RNAcounts_df = resolved_RNAcounts_df.groupby(groupby_columns).agg(agg_dict).reset_index()
	#print("after agg", resolved_RNAcounts_df[resolved_RNAcounts_df["N_tot"] == 0])
	#print("gene_name after agg", resolved_RNAcounts_df[resolved_RNAcounts_df["gene_name"] == 'AC007495.1'])

	#remove rows where N_in = N_tot (these will result in invalid value for log)
	resolved_RNAcounts_df = resolved_RNAcounts_df[resolved_RNAcounts_df["N_tot"]!=resolved_RNAcounts_df["N_in"]]

	resolved_RNAcounts_df['ag_reptypes'] = resolved_RNAcounts_df['ag_reptypes_str'].str.split(',')

	#get ag chroms
	# def extract_chrs(s):
	#	 return [x.split('_')[1] for x in s]

	def extract_chrs(s):
		if not s[0].startswith('chr'):
			return list(set(['chr'+x.split('_')[1] for x in s]))
		else:
			return list(set([x.split('_')[0] for x in s]))


	resolved_RNAcounts_df['ag_chrs'] = resolved_RNAcounts_df['ag_reptypes'].apply(extract_chrs)
	resolved_RNAcounts_df['ag_chrs_str'] = (resolved_RNAcounts_df['ag_chrs'].transform(lambda x: ",".join(map(str,x))))

	#print("after chr extraction", resolved_RNAcounts_df[resolved_RNAcounts_df["N_tot"] == 0])
	
	#calculate dpnII sites


	def calc_dpnII_count_rep(row):
		if row['ag_reptypes_str'].startswith('chr'):
			curr_gbin = row['ag_reptypes_str']
			tot_dpnII_count = sum(gbin_dpnII_df[gbin_dpnII_df['gbin'] == curr_gbin]['dpnII_count']) + 1
			return tot_dpnII_count
		else:
			reptypes = row['ag_reptypes']
			tot_dpnII_count = sum(rep_dpnII_df[rep_dpnII_df['reptype'].isin(reptypes)]['dpnII_count']) + 1
			return tot_dpnII_count

	def calc_dpnII_count_class(row):
		if row['ag_reptypes_str'].startswith('chr'):
			curr_gbin = row['ag_reptypes_str']
			tot_dpnII_count = sum(gbin_dpnII_df[gbin_dpnII_df['gbin'] == curr_gbin]['dpnII_count']) + 1
			return tot_dpnII_count
		else:
			reptypes = row['ag_reptypes']
			tot_dpnII_count = sum(class_dpnII_df[class_dpnII_df['class'].isin(reptypes)]['dpnII_count']) + 1
			return tot_dpnII_count


	#calculate Lin
	def calculate_total_length(row):
		if row['ag_reptypes_str'].startswith('chr'):
			return gbin_size
		else:
			reptypes = row['ag_reptypes']
			total_length = sum(bed_df[bed_df['reptype'].isin(reptypes)]['length'])
			return total_length
	
	def calculate_total_length_class(row):
		if row['ag_reptypes_str'].startswith('chr'):
			return gbin_size
		else:
			reptypes = row['ag_reptypes']
			total_length = sum(class_length_df[class_length_df["class"].isin(reptypes)]["length"])
			return total_length
	
	if level == "reptype":
		resolved_RNAcounts_df['L_in'] = resolved_RNAcounts_df.apply(calculate_total_length, axis=1)
		resolved_RNAcounts_df['dpnII'] = resolved_RNAcounts_df.apply(calc_dpnII_count_rep, axis=1)
	if level == "class":
		resolved_RNAcounts_df['L_in'] = resolved_RNAcounts_df.apply(calculate_total_length_class, axis=1)
		resolved_RNAcounts_df['dpnII'] = resolved_RNAcounts_df.apply(calc_dpnII_count_class, axis=1)
		
	#print("before L op", resolved_RNAcounts_df)
	#calculate L_out by subtracting L_in from L_tot
	resolved_RNAcounts_df['L_out'] = Ltotal - resolved_RNAcounts_df['L_in']

	#dpn norm
	resolved_RNAcounts_df["N_in_dnorm"] = resolved_RNAcounts_df["N_in"] / resolved_RNAcounts_df['dpnII']
	
	#calculate raw T_score
	resolved_RNAcounts_df["raw_E"] = resolved_RNAcounts_df["N_in_dnorm"] / (resolved_RNAcounts_df["N_tot"] - resolved_RNAcounts_df["N_in"])
	
	#calculate T score
	resolved_RNAcounts_df["Escore"] = resolved_RNAcounts_df["raw_E"] * (resolved_RNAcounts_df["L_out"]/resolved_RNAcounts_df["L_in"])


	#print negative and inf and Na Tscores
	neg_E = resolved_RNAcounts_df[resolved_RNAcounts_df["Escore"] < 0]
	if not neg_E.empty:
		print("Negative Escores:", neg_E[['ag_reptypes_str','gene_name','L_in','L_out','N_in','N_tot','raw_E','Escore']])

	rows_with_inf = resolved_RNAcounts_df[np.isinf(resolved_RNAcounts_df["Escore"])]
	if not rows_with_inf.empty:
		print("Inf Escores:", rows_with_inf[['ag_reptypes_str','gene_name','L_in','L_out','N_in','N_tot','raw_E','Escore']])

	rows_with_na = resolved_RNAcounts_df[resolved_RNAcounts_df["Escore"].isna()]
	if not rows_with_na.empty:
		print("NA Escores:", rows_with_na[['ag_reptypes_str','gene_name','L_in','L_out','N_in','N_tot','raw_E','Escore']])


	
	#calculate log(T score)
	resolved_RNAcounts_df["logEscore"] = np.log2(resolved_RNAcounts_df["Escore"])

	#add logNtot for resid & plot 
	resolved_RNAcounts_df['logNtot'] = np.log(resolved_RNAcounts_df["N_tot"])

	
	#rearrange order of columns to make it pretty & output to a file
	resolved_RNAcounts_df = resolved_RNAcounts_df[['ag_reptypes', 'ag_reptypes_str', 'ag_chrs', 'ag_chrs_str', 'gene_ID', 'gene_name', 'chrom','type','N_in','N_tot','dpnII','logNtot','L_out','L_in','raw_E','Escore','logEscore']]
	
	if save:
		#g = str(int(gbin_size / 1000)) + 'kb'
		ofile_name = parent_dir + "dna.fastq.gz.ANNOTATED.ci2_and_unqv2." + level + "." + RNA_type + "." + g + ".Escore"
		#stats_ofile_name= parent_dir + "dna.fastq.gz.ANNOTATED.ci2_and_unqv2" + "." + level + "." + RNA_type + "T.resolved_ag.stats"
		df_out=resolved_RNAcounts_df.copy()
		df_out['ag_reptypes']=df_out['ag_reptypes'].apply(lambda x: str(x))
		df_out['ag_chrs']=df_out['ag_chrs'].apply(lambda x: str(x))
		df_out.to_csv(ofile_name, sep="\t")
		#stats_df.to_csv(stats_ofile_name, sep="\t")
		print("Saved Edf as",ofile_name)
		#print("Saved resolved ag stats as ",stats_ofile_name)
	return resolved_RNAcounts_df



def model_and_predict(pc_df,all_df):
	adf = all_df.copy()

	X = pc_df[['logNtot', 'L_in', 'dpnII']].values
	y = pc_df['logEscore'].values
	#X_test = adf[['logNtot','L_in','dpnII']].values
	

	gam = LinearGAM().fit(X,y)
	#print(gam.summary())
	rsq = gam.statistics_['pseudo_r2']['explained_deviance']
	#print(f'R-squared value: {r2:.4f}')

	adf['predicted_logEscore'] = gam.predict(adf[['logNtot','L_in','dpnII']])
	return adf, rsq



#test model & predict with split dataframes
def gettransRNAs(chrom_list, exons_df):
	
	notchrom_df = exons_df[~exons_df["chrom"].isin(chrom_list)]
	
	
		#notchrom_df = full_df[full_df["chrom"] != chrom]
		#pd.concat([combined_df,notchrom_df]).drop_duplicates(keep='first')
	return(notchrom_df)


def compute_GAM_residual(in_df, in_ex_df, max_ag_chr_val):
	
	tmp_df=in_df.copy()
	df_out = tmp_df.drop(tmp_df[tmp_df['N_tot'] == tmp_df['N_in']].index)
	resid_dict={}

	df_out = df_out[df_out['ag_chrs'].apply(lambda x: len(x) < max_ag_chr_val)]

	tmp_ex_df=in_ex_df.copy()
	ex_df=tmp_ex_df.drop(tmp_ex_df[tmp_ex_df['N_tot'] == tmp_ex_df['N_in']].index)
	ex_df['ag_chrs'] = ex_df['ag_chrs_str'].str.split(',')
	ex_df = ex_df[ex_df['ag_chrs'].apply(lambda x: len(x) < max_ag_chr_val)]
	
	model_stats = []
	for i in df_out["ag_chrs_str"].unique():
		if i not in resid_dict.keys():
			i_df = df_out[df_out['ag_chrs_str']==i]
			chrom_list = i_df['ag_chrs'].iloc[0]

			i_ex_df = ex_df[ex_df['ag_chrs_str']==i]
			
			if len(chrom_list) > 1:
				coi = ex_df[ex_df['ag_chrs'].apply(lambda x: all(c in chrom_list for c in x))]
				ex_pc_df = coi.query('(type == "protein_coding") & (N_in>1)')
			else:
				ex_pc_df = i_ex_df.query('(type == "protein_coding") & (N_in>1)')
				
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


				result, rsq = model_and_predict(bkg_data,i_df)
				result_df = result.copy()
				result_df['R'] = result_df['logEscore'] - result_df['predicted_logEscore']
				resid_dict[i]=result_df
				model_stats.append({'chr_group':i, 'rsquared':rsq, 'nSamples':ex_pc_offchr_df.size})
	
	model_stats_df = pd.DataFrame(model_stats)
	resid_dict_filt = {}
	for k in resid_dict.keys():
		df = resid_dict[k]
		if not df.empty:
			resid_dict_filt[k]=df
	

	#resid_dict_filt = {k: v for k, v in resid_dict.items() if v is not None}
	big_df = pd.concat(resid_dict_filt.values())



	return big_df, model_stats_df


def compute_residuals_lr(in_df, in_ex_df, max_ag_chr_val):
	
	tmp_df=in_df.copy()
	df_out = tmp_df.drop(tmp_df[tmp_df['N_tot'] == tmp_df['N_in']].index)
	resid_dict={}

	df_out = df_out[df_out['ag_chrs'].apply(lambda x: len(x) < max_ag_chr_val)]

	tmp_ex_df=in_ex_df.copy()
	ex_df=tmp_ex_df.drop(tmp_ex_df[tmp_ex_df['N_tot'] == tmp_ex_df['N_in']].index)
	ex_df['ag_chrs'] = ex_df['ag_chrs_str'].str.split(',')
	ex_df = ex_df[ex_df['ag_chrs'].apply(lambda x: len(x) < max_ag_chr_val)]
	
	model_stats = []
	for i in df_out["ag_chrs_str"].unique():
		if i not in resid_dict.keys():
			i_df = df_out[df_out['ag_chrs_str']==i]
			chrom_list = i_df['ag_chrs'].iloc[0]

			i_ex_df = ex_df[ex_df['ag_chrs_str']==i]
			
			if len(chrom_list) > 1:
				coi = ex_df[ex_df['ag_chrs'].apply(lambda x: all(c in chrom_list for c in x))]
				ex_pc_df = coi.query('(type == "protein_coding") & (N_in>1)')
			else:
				ex_pc_df = i_ex_df.query('(type == "protein_coding") & (N_in>1)')
				
			if len(ex_pc_df) < 10:
				print(i, "not enough data")
				resid_dict[i] = "Drop"
			else:
				ex_pc_offchr_df = gettransRNAs(chrom_list,ex_pc_df)

				if len(ex_pc_offchr_df) < 10:
					print(str(i) + " Does not contain enough samples. Using all protein coding RNAs as background")
					bkg_data = gettransRNAs([],ex_pc_df)
				else:
					bkg_data = ex_pc_offchr_df.copy()


				#result, rsq = model_and_predict(bkg_data,i_df)
				result = linregress(np.log(bkg_data["N_tot"]), bkg_data["logEscore"])
				resid_dict[i] = [result.intercept, result.slope]
	drop_these = []
	for k in resid_dict.keys():
		if resid_dict[k] == "Drop":
			drop_these.append(k)
	for d in drop_these:
		resid_dict.pop(d, None)
	df_out = df_out[~df_out['ag_chrs_str'].isin(drop_these)]

	resid_df = pd.DataFrame.from_dict(resid_dict, orient = 'index', columns=["intercept","slope"])
	resid_df=resid_df.rename_axis('ag_chrs_str').reset_index()

	R_df=df_out.merge(resid_df, how="left", on="ag_chrs_str")
	R_df['predictedEscore'] = R_df['intercept'] + R_df['slope']  * np.log(R_df["N_tot"])
	R_df['R']= R_df["predictedEscore"] - R_df['logEscore'] 


	return R_df