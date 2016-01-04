import pandas as pd
import numpy as np
import warnings
from sklearn.ensemble import GradientBoostingRegressor
import math
import re
import matplotlib.pyplot as plt 
from scipy.stats import combine_pvalues
from collections import OrderedDict
import urllib,urllib2
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from matplotlib.pyplot import figure, show
import seaborn
warnings.filterwarnings('ignore')




class Report(object):
	def __init__(self,name):
		self.name = name
		self.histo_figure = plt.figure()
		self.scatter_figure = plt.figure()
		self.scatter_counter=1


	def plot_scatter(self,df,cols):
		fig = self.scatter_figure
		vals = df[cols].get_values()
		vals[~np.isnan(vals).any(axis=1)]
		ax = fig.add_subplot(111)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.axes.grid(False)
		ax = fig.add_subplot(int(math.sqrt(len(group_names)))+1,int(math.sqrt(len(group_names)))+1,self.scatter_counter);
		vals[:,0]=[~vals[:,0]>np.std(vals[:,0])]
		vals[:,1]=[~vals[:,1]>np.std(vals[:,1])]
		ax.scatter(vals[:,0],vals[:,1])
		self.scatter_counter+=1
		
		

	def plot_histo(self,df,cols):
		cntr=1;cntr2=1
		vals = df[cols].get_values()
		fig = self.histo_figure
		ax = fig.add_subplot(111)
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.axes.grid(False)
		for col in cols:
			ax = fig.add_subplot(int(math.sqrt(len(group_names)))+1,int(math.sqrt(len(group_names)))+1,cntr); cntr+=1
			ax.hist(df[col].values,bins=50)
			ax.set_title(col)
			ax.set_xlim([np.min(vals),np.max(vals)])
		fig.tight_layout()
		plt.show()










##########################BEGIN RUN
#list columns in a course to calculate the differential expression of
report_name = "CDK8_inhibited_run"
phospho_analysis = True
protein_analysis = False
silac = True
project_outpath = "/home/jjacob/python/zposs/Phospho_MQ_search/report/"
mqrun_path = '/home/jjacob/python/zposs/Phospho_MQ_search/'
ignore_missing_replicate_values = True
replicate_method = 'mean'
header_suffix = 'Ratio H/L normalized'
exp_design = [['1','2','3']]
inverse = [False,False,True]
group_names = ['treatment_1h']
experiments = {}





#list columns in a course to calculate the differential expression of
#report = Report("CDK8_inhibited_run")
#phospho_analysis = True
#protein_analysis = False
#silac = True
#project_outpath = "./report/"
#mqrun_path = '/home/jjacob/python/zposs/Proteome/txt/'
#ignore_missing_replicate_values = True
#replicate_method = 'mean'
#header_suffix = 'Ratio H/L normalized'
#exp_design = [['1','7'],['2','8'],['3','9'],['4','10'],['5','11'],['6','12']]
#inverse = [False,True]
#group_names = ['DMSO','treatment_1h','treatment_3h','treatment_6h','treatment_18h','treatment_24h']
#experiments = {}












class UniQuery(object):
	def __init__(self,name):
		self.name = name
	def queryUniprotGenes(self,df,accession_col):
		accessions = ""
		uniprot_accessions = df[accession_col].str.split("-")
		uniprot_accessions = [i[0] for i in uniprot_accessions]
		for accession in uniprot_accessions:
			accessions = accessions+"\t"+accession
		url = 'http://www.uniprot.org/mapping/'
		params = {}
		params['from'] = 'ACC'
		params['to'] = 'GENENAME'
		params['format'] = 'tab'
		params['query'] = accessions
	
		#params = {'from':'ACC','to':'P_REFSEQ_AC','format':'tab','query':'P13368 P20806 Q9UM73 P97793 Q17192'}
		data = urllib.urlencode(params)
		request = urllib2.Request(url, data)
		contact = "" # Please set your email address here to help us debug in case of problems.
		request.add_header('User-Agent', 'Python %s' % contact)
		response = urllib2.urlopen(request)
		page = response.read(200000)

		p=page.split("\n")
		map_out = {}
		for val in p:
			if len(val.split("\t"))>1:
				map_out[val.split("\t")[0]]=val.split("\t")[1]
		arr_out = []
		for acc in uniprot_accessions:
			if acc in map_out:
				arr_out.append(map_out[acc])	
			else:
				arr_out.append("no_match")
		df.insert(1,"gene_name",arr_out)
		return df












def impute(df,imp_val,headers):
	if np.isnan(imp_val):
		imp_val = -500	
	print "imputing..."

	model = GradientBoostingRegressor(loss='ls', learning_rate=0.1, n_estimators=100, subsample=1.0, min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0, max_depth=3, init=None, random_state=None, max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None, warm_start=False, presort='auto')
	data = np.array(df[headers].get_values())
	data[np.isnan(data)] = -500

	for col in range(0,len(headers)):
		#print "Working on column: "+str(col)
		##for the current column, remove rows where the current (row,column) value is not equal to zero
		##this way we are only training on data with non-zero target values
		reduced_data = data[np.logical_not(data[:,col] == imp_val)] #remove row if row,col_num value is zero
		target_set = reduced_data[:,col]
		training_set = np.delete(reduced_data,col,1)
		model.fit(training_set,target_set)
		row_num=0
		for row in data:
			remaining = np.delete(row,col,0)
			if data[row_num,col] == imp_val:
				data[row_num,col] = model.predict(remaining)
			row_num+=1
	cntr=0
	for h in headers:
		df[h] = data[:,cntr];cntr+=1
	return df










def isFloat(val):
	if float(val):
		return True
	else:
		return False










def normcdf(x, mu, sigma):
	cdf = 0.5 * (1 + math.erf((x - mu)/math.sqrt(2 * sigma**2)))
	return min(cdf,1-cdf)










def bh_correct(mapping): #apply benjamini-hochberg correction to values in a dict
	ordered_mapping = OrderedDict(sorted(mapping.items(),key=lambda t: t[1]))
	cntr=1
	output = []
	for key in ordered_mapping:
		mapping[key] = mapping[key]*len(mapping)/cntr
		cntr+=1
	return mapping








def l2fc_transform(df,columns):
	vals = df[columns].get_values()
	vals[~np.isnan(vals)] = np.log2(vals[~np.isnan(vals)])	
	for c in range(len(columns)):
		df[columns[c]] = vals[:,c]
	return df









def append_pval(df,columns,multitest):
	pvalcols = []
	for col in columns:
		pvals = [normcdf(val,np.mean(df[col].values),np.std(df[col].values)) for val in df[col].values]
		df.insert(len(df.columns),col+"_pvals",pvals)
		pvalcols.append(col+"_pvals")
	if multitest:
		combined=[]
		for row in df[pvalcols].get_values():
			combined.append(combine_pvalues(row,method='fisher', weights=None)[1])
		df.insert(len(df.columns),'fisher_combined_pval',combined)
		bh_corrected = bh_correct(dict(zip(df.index.values,df['fisher_combined_pval'].values)))
		corrected_vals = []
		for k in df.index.values:
			corrected_vals.append(bh_corrected[k])
		df.insert(len(df.columns),'benj_hoch_corrected_pval',corrected_vals)
	return df








#this applies a filter to the entire dataset based on the headers argument
def apply_filter(ftype,f_threshold,headers,df): ## all filters are inclusive
	if ftype == 'remove_nan':
		df = df[pd.notnull(df[headers[0]])]
	if ftype == 'impute_nan':
		all_vals = df[headers].get_values()
		all_vals_no_nans = all_vals[~np.isnan(all_vals)]
		for header in headers:
			vals = df[header].values
			vals_no_nans = vals[~np.isnan(vals)]
			samples = np.random.normal(np.mean(all_vals_no_nans),np.std(all_vals_no_nans)+.0001,vals[np.isnan(vals)].shape[0])
			neg_samples=True
			while neg_samples:
				if samples[samples<0].shape[0]==0:
					neg_samples=False
				else:
					samples[samples<0] = np.random.normal(np.mean(all_vals_no_nans),np.std(all_vals_no_nans)+.0001,samples[samples<0].shape[0])
	
			vals[np.isnan(vals)] = samples
			df[header] = vals
	if ftype == 'highpass':
		for header in headers:
			df = df[df[header]>=f_threshold]
	if ftype == 'lowpass':
		for header in headers:
			df = df[df[header]<=f_threshold]
	if ftype == "contaminants":
		df[headers[0]] = df[headers[0]].apply(str)
		df = df[~df[headers[0]].str.contains("CON__")]
	if ftype == "reverse":
		df[headers[0]] = df[headers[0]].apply(str)
		df = df[~df[headers[0]].str.contains("REV__")]
	if ftype == 'collapse':
		df = df[df[headers].isnull().sum(axis=1)<=f_threshold]

	return df




def generate_protein_series(headers,df,additional_headers):
	pids = ['Majority protein IDs']	
	df = df[pids+headers+additional_headers]
	return df




def expand_table(df,col,sep):
	exp_vals = df[col].values
	everything_else = df.drop(col,axis=1).get_values()
	out = []
	names = []
	h2idx = {}
	for c in range(len(df.columns)):
		h2idx[df.columns[c]] = c

	#make numpy matrix for value storage
	cntr=0
	for vals in exp_vals:
		for name in vals.split(";"):
			names.append(name)
			out.append(everything_else[cntr,:])
		cntr+=1
	dfout = pd.DataFrame(data=np.array(out),columns=df.drop(col,axis=1).columns)
	dfout.insert(h2idx[col],col,names)
	return dfout








			


def generate_phospho_series(headers,df): #phospho specific
	modseqs = df['Modified sequence'].values
	nsites = df['Number of Phospho (STY)'].values
	positions = df['Positions'].values
	AA = df['Amino acid'].values
	pgroups = df['Leading proteins'].values
	ids = df['id'].values
	uniprot_accessions = []
	phospho_sequence = []
	phospho_site = []
	multiplicity = []
	aa = []
	phospho_ids = []

	series_responses = [[] for g in group_names]
	additional_columns = [[] for c in range(4)]
	phosphosites=[]
	psite2vals = {}
	for i in range(0,len(modseqs)):
		all_sites = nsites[i].split(";")
		for site in all_sites:
			if int(site)>3:
				site=3
			time_set = []
			cntr=0
			for g in group_names:
				series_responses[cntr].append(df[g+"_site_"+str(site)].values[i])
				cntr+=1
			uniprot_accessions.append(pgroups[i])
			phospho_sequence.append(modseqs[i])
			phospho_site.append(str(positions[i]))
			aa.append(AA[i])
			multiplicity.append(str(site))
			phosphosites.append(pgroups[i]+modseqs[i]+AA[i]+str(positions[i])+"_"+str(site))
			phospho_ids.append(ids[i])
	
	series_df = pd.DataFrame()
	series_df.insert(0,"id",phospho_ids)
	series_df.insert(1,"Phosphosites",phosphosites)
	series_df.insert(2,"uniprot_accession",uniprot_accessions)
	series_df.insert(3,"phospho_sequence",phospho_sequence)
	series_df.insert(4,"amino_acid",aa)
	series_df.insert(5,"phospho_site",phospho_site)
	series_df.insert(6,"multiplicity",multiplicity)
	cntr=0
	for g in group_names:
		series_df.insert(series_df.shape[1],g,series_responses[cntr-1]); cntr+=1
	series_df = series_df.drop_duplicates(subset = ["Phosphosites"])
	series_df.drop('Phosphosites',axis=1,inplace=True)
	return series_df
		

















def generate_phospho_path_inp(df,phosphosite_col,series_cols):
	dic = {}
	nodes = set()
	for site in df[phosphosite_col].values:
		split = site.split("_")
		if "-" in split[0]:
			continue
		
		name = split[0]+"-"+split[2]+"-"+split[3]
		nodes.add(split[0]+"\t"+split[0]+"-"+split[2]+"\t"+split[2]+"\t"+split[3]+"\n")

		for values in df[df[phosphosite_col]==site][series_cols].values:
			cntr=1
			for val in values:
				dic[name+"-"+str(cntr)] = val;cntr+=1
	fout = open('./phospho_path_nodes.txt','w')
	for val in nodes:
		fout.write(val)
	fout.close()
	df_pp=pd.DataFrame()
	df_pp.insert(0,'ID',dic.keys())
	df_pp.insert(1,'Ratio',dic.values())
	df_pp.to_csv("./phospho_path_time_series.txt",sep="\t",index=False)








def columns_normalize(df,cols):
	for col in cols:
		coldata = np.array(df[col].get_values().astype(float))
		##..unity normalize
		if coldata.min()==coldata.max():
			continue
		normalized = np.divide(np.subtract(coldata,coldata.min()),(coldata.max()-coldata.min()))
		df[col] = normalized
	return df





def dfappend(df,series_df,cols):
	for col in cols:
		series_df.insert(series_df.shape[1],col,df[col][series_df['id'].values].values)
	return series_df









if phospho_analysis:
	report = Report(report_name)
	df = pd.read_csv(mqrun_path+"Phospho (STY)Sites.txt", sep = "\t")
	all_site_cols = []
	std_err_cols = []

	for exp_idx in range(0,len(exp_design)):
		all_sites = []
		for multiplicity in range(1,4):
			site = []
			for rep in exp_design[exp_idx]:
				site.append(header_suffix+' '+rep+'___'+str(multiplicity))
			all_sites.append(site)
		group_vals = []
		site_cols = []
		site_counter = 1
		for site in all_sites: ## invert sites where ratio is reversed
			cntr = 0
			for col in site:
				if inverse[cntr]:	
					df[col] = 1/df[col].get_values()
				cntr+=1
			if replicate_method=='mean':
				df.insert(df.shape[1],group_names[exp_idx]+"_site_"+str(site_counter),df[site].mean(axis=1).values)
				df.insert(df.shape[1],group_names[exp_idx]+"_site_cov_"+str(site_counter),np.divide(df[site].std(axis=1).values,df[site].mean(axis=1).values))
				std_err_cols.append(group_names[exp_idx]+"_site_cov_"+str(site_counter))
				site_cols.append(group_names[exp_idx]+"_site_"+str(site_counter))
				all_site_cols.append(group_names[exp_idx]+"_site_"+str(site_counter))
			site_counter+=1

	df = apply_filter('remove_nan',None,['Number of Phospho (STY)'],df) 
	df = apply_filter('contaminants',None,['Leading proteins'],df) 
	df = apply_filter('reverse',None,['Leading proteins'],df) 

	df = l2fc_transform(df,all_site_cols)
	series_df = generate_phospho_series(all_site_cols,df)	
	series_df = apply_filter('impute_nan',None,group_names,series_df)	
	series_df = append_pval(series_df,group_names,multitest=True)

	series_df = apply_filter('lowpass',0.05,['benj_hoch_corrected_pval'],series_df)
	series_df = dfappend(df,series_df,std_err_cols)
	series_df = dfappend(df,series_df,['Gene names'])

	report.plot_histo(series_df,group_names)

	uni_query = UniQuery("gene_query")
	series_df = uni_query.queryUniprotGenes(series_df,'uniprot_accession')
	series_df = apply_filter('highpass',1.0,group_names,series_df)
	#generate_phospho_path_inp(series_df,'Phosphosites',group_names[1:])
	series_df.to_csv(project_outpath+"phosphosite_diffexp.txt",sep="\t",index=False)













if protein_analysis:
	df = pd.read_csv(mqrun_path+"proteinGroups.txt", sep = "\t")
	all_protein_ratio_cols = []
	std_errs_cols = []
	
	for exp_idx in range(0,len(exp_design)):
		groups = []
		for h in range(0,len(exp_design[exp_idx])):
			if inverse[h]==True:
				df[header_suffix+' '+exp_design[exp_idx][h]] = 1/df[header_suffix+' '+exp_design[exp_idx][h]].values
			groups.append(header_suffix+' '+exp_design[exp_idx][h])
		if replicate_method=='mean':
			df.insert(df.shape[1],group_names[exp_idx]+"_Ratio_averaged",df[groups].mean(axis=1).values)
			df.insert(df.shape[1],group_names[exp_idx]+"_Ratio_stderr",df[groups].std(axis=1).values)
			all_protein_ratio_cols.append(group_names[exp_idx]+"_Ratio_averaged")
			std_errs_cols.append(group_names[exp_idx]+"_Ratio_stderr")	
			report.plot_scatter(df,groups)


	df.insert(df.shape[1],'err_sum',df[std_errs_cols].sum(axis=1).values)
	std_errs_cols.append('err_sum')
	df = apply_filter('collapse',1,all_protein_ratio_cols,df) #only allow 1 points to be missing (from both replicates)
	df = apply_filter('contaminants',None,['Majority protein IDs'],df) 
	df = apply_filter('reverse',None,['Majority protein IDs'],df) 
	series_df = generate_protein_series(all_protein_ratio_cols,df,std_errs_cols)
	series_df = impute(series_df,np.nan,all_protein_ratio_cols+std_errs_cols)
	report.plot_histo(series_df,all_protein_ratio_cols)
	series_df = append_pval(series_df,all_protein_ratio_cols,multitest=True)
	series_df = append_max_lf2c(series_df,all_protein_ratio_cols)
	series_df = apply_filter('lowpass',0.05,['benj_hoch_corrected_pval'],series_df)
	series_df = apply_filter('highpass',1.5,['max_log2_change'],series_df)
	#series_df = expand_table(series_df,'Majority protein IDs',sep=";")
	#uni_query = UniQuery("gene_query")
	#series_df = uni_query.queryUniprotGenes(series_df)
	#series_df = series_df.drop_duplicates(subset = ["Majority protein IDs"])
	#correl_mat(series_df)
	series_df.to_csv(project_outpath+"/protein_diffexp.txt",sep="\t",index=False)







			



