##This version is specific to the 6 month review for analysis of EGF 13 time point Silac data 


import warnings
warnings.filterwarnings('ignore')

import sys
from scipy.stats import ranksums
from scipy.stats import wilcoxon
from scipy.stats import sem
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.neighbors import NearestNeighbors
import ast
import math
import re
import matplotlib.pyplot as plt 
from scipy.stats import combine_pvalues
from collections import OrderedDict
import urllib,urllib2
from scipy.stats import pearsonr
from sklearn.cluster import KMeans
from matplotlib.pyplot import figure, show
import os
import seaborn
from collections import OrderedDict
from sklearn.preprocessing import Imputer
from scipy.stats import ttest_ind
warnings.filterwarnings("ignore", category=DeprecationWarning) 





class DesignParser:
	def __init__(self,path):
		self.design_path = path
		self.parse(path)
	def parse(self,path):
		tree = ET.parse(path)
		root = tree.getroot()
		runtype = root.iter('Run').next().get('type')
		runfile = root.iter('Run').next().get('file')
		suffix = root.iter('Run').next().get('suffix')
		spec = root.iter('Run').next().get('spec')

		self.run = Run(root.get('Description'),runtype,runfile,suffix,spec) #instantiate run object and append description
		for cond in root.iter('Condition'):  #initiate all condition objects
			cond_name = cond.attrib.get('name')
			cond_type = cond.attrib.get('type')
			cond_desc = cond.attrib.get('description')
			condition = self.run.condition_map[cond_name] = Condition(cond_name,cond_type,cond_desc)
			self.run.conditions.append(cond_name)
			self.run.condition_map[cond_name] = condition
			
			for rep in cond.iter('Replicate'):
				rep_name = rep.attrib.get('name')
				rep_type = rep.attrib.get('type')
				rep_desc = rep.attrib.get('description')
				inv_bool = rep.attrib.get('inverted')
				supergroup = rep.attrib.get('supergroup')
				new_replicate = Replicate(rep_name,rep_type,rep_desc,inv_bool,supergroup)


				if supergroup in self.run.supergroups:
					self.run.supergroups[supergroup].update_condmap(cond_name,new_replicate)
				else:
					self.run.supergroups[supergroup] = SuperGroup(supergroup,cond_name,new_replicate)


				condition.replicates.append(new_replicate)
				self.run.replicates[rep_name]=new_replicate
				for trep in rep.iter('TReplicate'):
					trep_name = trep.attrib.get('name')
					rep_type = rep.attrib.get('type')
					rep_desc = rep.attrib.get('description')
					new_trep = TReplicate(rep_name,rep_type,rep_desc)
					self.run.replicates[rep_name].treps[trep_name] = new_trep

		for ratios in root.iter("Ratios"):
			for ratio in ratios.iter("Ratio"):
				self.run.ratios.append(ratio.attrib.get("des"))

		for f in root.iter("filter"):
			subtype = f.attrib.get('subtype')
			ftype = f.attrib.get('ftype')
			threshold = f.attrib.get('threshold')
			headers = f.attrib.get('headers')
			function = f.attrib.get('function')
			self.run.filters.append(Filter(subtype,ftype,threshold,headers,function))
				
		for header in root.iter("header"):
			self.run.headers.append(header.attrib.get('name'))
		for t in root.iter("normalization"):
			self.run.norm = Norm(t.attrib.get('type'),t.attrib.get('scope'))
		for t in root.iter("appendGenes"):
			self.run.appendGenes = t.attrib.get('val')

	
		for T in root.iter("Tests"):
			for t in T.iter("test"):
				ttype=t.attrib.get("type")
				c1=t.attrib.get("Control")
				c2=t.attrib.get("Treatment")
				self.run.tests.append(Test(ttype,c1,c2))



class Run:
	def __init__(self,description,runtype,runfile,suffix,spec):
		self.description = description
		self.runtype = runtype
		self.runfile = runfile
		self.suffix = suffix
		self.spec = spec
		self.condition_map = OrderedDict() #head node for run traversal
		self.ratios = []
		self.conditions  = []
		self.condition_map = {}
		self.supergroups = {}
		self.replicates = {}
		self.filters = []
		self.headers = []
		self.c2RNames = {}
		self.norm = None
		self.appendGenes = "False"
		self.tests = []

	def c2repnames(self):
		for c in self.conditions:
			reps = self.condition_map[c].replicates
			self.c2RNames[c] = [(self.suffix+rep.name) for rep in reps]
			

class Test:
	def __init__(self,ttype,c1,c2):
		self.ttype=ttype
		self.c1=c1
		self.c2=c2

class Norm:
	def __init__(self,t,scope):
		self.type=t
		self.scope=scope

class Filter:
	def __init__(self,subtype,ftype,threshold,headers,function):
		self.subtype = subtype
		self.ftype = ftype
		self.threshold = threshold
		self.headers = headers
		self.function = function

class Xml_obj():
	def build(self,name,Ttype,description):
		self.name = name
		self.type = Ttype
		self.description = description
		self.attributes = {'name':self.name,'type':self.type,'desc':self.description}



class Condition(Xml_obj):
	def __init__(self,name,Ttype,description):
		self.build(name,Ttype,description)
		self.replicates = []
	

class Replicate(Xml_obj):
	def __init__(self,name,Ttype,description,inverted,sg):
		self.build(name,Ttype,description)
		self.name = name
		self.inverted = ast.literal_eval(inverted)
		self.supergroup = sg
		self.treps = {}
		

class TReplicate(Xml_obj):
	def __init__(self,name,Ttype,description):
		self.build(name,Ttype,description)

class SuperGroup():
	def __init__(self,name,condition_name,rep):
		self.name = name
		self.repmap = {}
		self.update_condmap(condition_name,rep)
	def update_condmap(self,condition_name,rep):
		self.repmap[condition_name] = rep







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





def log(string,mtype):
	if mtype==1:
		print "## "+string
	if mtype==2:
		print "!! "+string



def isFloat(val):
	if float(val):
		return True
	else:
		return False






#this applies a filter to the entire dataset based on the headers argument
def apply_filter(ftype,f_threshold,headers,df): ## all filters are inclusive
	#if f_threshold == "None":
	#	print "here"
	#	f_threshold = None
	if ftype == 'by_matching':
		for header in headers:
			df = df[df[header]!="By matching"]

	if ftype == 'remove_nan':
		for header in headers:
			df = df[pd.notnull(df[header])]
	if ftype == 'impute_nan':
		all_vals = df[headers].get_values()
		all_vals_no_nans = all_vals[~np.isnan(all_vals)]
		for header in headers:
			vals = df[header].values
			notnan = ~np.isnan(vals)
			vals_no_nans = vals[~np.isnan(vals)]
			vals_no_nans_sum = notnan.shape
			samples = np.random.normal(np.mean(all_vals_no_nans),np.std(all_vals_no_nans)+.0001,vals[np.isnan(vals)].shape[0])
			vals[np.isnan(vals)] = samples
			df[header] = vals
	if ftype == 'highpass':
		for header in headers:
			df = df[df[header]>=f_threshold]
	if ftype == 'lowpass':
		for header in headers:
			df = df[df[header]<=f_threshold]
	if ftype == 'abspass_low':
		for header in headers:
			df = df[abs(df[header])<=f_threshold]
	if ftype == 'abspass_high':
		for header in headers:
			df = df[abs(df[header])>=f_threshold]
	if ftype == "contaminants":
		log("Filtering Contaminants.",1)
		df[headers] = df[headers].apply(str)
		df = df[~df[headers].str.contains("CON__")]
	if ftype == "reverse":
		log("Filtering Reverse Sequences.",1)
		df[headers] = df[headers].apply(str)
		df = df[~df[headers].str.contains("REV__")]
	if ftype == 'collapse_nans':
		df = df[df[headers].isnull().sum(axis=1)<=f_threshold]
	if ftype == 'collapse_zeros':
		df = df[(df[headers]==0.0).sum(axis=1)<=f_threshold]
	if ftype == 'collapse_on_value':
		for header in headers:
			df = df[df[header]<=f_threshold]	
	if ftype == 'zeros_to_nans':
		for header in headers:
			df[header] = df[header].replace(['0'], np.nan)
	return df
	












def generate_series(analysis_type,df):
	log("Generating Expression Series",1)
	runtype = design.run.runtype
	group=[]	
	ids = ['id']
	leading_accessions = []	
	ratios = [ratcol+"^ratio" for ratcol in design.run.ratios]

	##attempt to parse out accessions
	if "__" not in df['Leading proteins'].values[0]:
		if analysis_type=='mod':
			try:
				leading_accessions = [val.split("|")[1].split("-")[0] for val in df['Leading proteins'].values]
			except:
				print "Could not parse leading accessions.  Either they are already parsed, or contaminants/reversed/nan accession(s) are present."
		else:
			try:
				leading_accessions.append(val.split("|")[1].split("-")[0])
			except:
				print "Could not parse leading accessions.  Either they are already parsed, or contaminants/reversed/nan accession(s) are present."

	

	series_df = pd.DataFrame()
	series_df.insert(0,"leading_accession",leading_accessions)
	[series_df.insert(series_df.shape[1],header,df[header].values) for header in design.run.headers]
	[series_df.insert(series_df.shape[1],col,df[col].values) for col in ids+group]
	[series_df.insert(series_df.shape[1],design.run.suffix+rep.name,df[design.run.suffix+rep.name].values) for c in design.run.conditions for rep in design.run.condition_map[c].replicates] # insert raw measurements
	raw_measurement_labels = [design.run.suffix+rep.name for c in design.run.conditions for rep in design.run.condition_map[c].replicates]

	'''
	Count up number of missing replicates in each condition
	'''
	for c in design.run.conditions: 
			series_df.insert(series_df.shape[1],c+"^rep_observation_count",(series_df[[design.run.suffix+rep.name for rep in design.run.condition_map[c].replicates]]).notnull().sum(axis=1))

	if runtype=="Silac":
		'''
		Append some columns that we will return from this method in order to perform later operations such as filtering and calculation of p-values
		'''
		n_missing = [c+"^rep_observation_count" for c in design.run.conditions]
		errcols = [c+"^stderr" for c in design.run.conditions]
		log2groups = {}
		for c in design.run.conditions:
			log2vals = []
			for rep in design.run.condition_map[c].replicates:
				log2vals.append(design.run.suffix+rep.name+"^log2vals")
			log2groups[c] = log2vals

		'''
		Invert replicate if appropriate
		'''
		for rep in design.run.replicates:
			if design.run.replicates[rep].inverted:
				series_df[design.run.suffix+rep] = 1/(df[design.run.suffix+rep].values)
		[series_df.insert(series_df.shape[1],design.run.suffix+rep.name+"^log2vals",series_df[design.run.suffix+rep.name].astype(float).apply(np.log2)) for c in design.run.conditions for rep in design.run.condition_map[c].replicates]
		log2_measurement_labels = [design.run.suffix+rep.name+"^log2vals" for c in design.run.conditions for rep in design.run.condition_map[c].replicates]

		
		#rep1headers = [("Ratio H/L normalized "+str(v)) for v in range(1,27)]
		#rep2headers = [("Ratio H/L normalized "+str(v)) for v in range(14,27)]
		#series_df=adv_impute(series_df,[(h+"^log2vals")for h in rep1headers])
		#series_df=adv_impute(series_df,[(h+"^log2vals")for h in rep2headers])
		


		
		##Calculate averages for each condition.  We want these to be mapped so the output is in block format
		##convert to log space first, otherwise the distribution is not normal and all other calculations will be unreliable
		for c in design.run.conditions:
			series_df.insert(series_df.shape[1],c+"^log2_average",(series_df[[design.run.suffix+rep.name+'^log2vals' for rep in design.run.condition_map[c].replicates]]).mean(axis=1).get_values())
		for c in design.run.conditions:
			series_df.insert(series_df.shape[1],c+"^log2_stdev",(series_df[[design.run.suffix+rep.name+'^log2vals' for rep in design.run.condition_map[c].replicates]]).std(axis=1).get_values())
		for c in design.run.conditions:
			series_df.insert(series_df.shape[1],c+"^log2_cov",(series_df[c+"^log2_average"].subtract(series_df[c+"^log2_stdev"]).multiply(series_df[c+"^log2_average"]).apply(np.abs).get_values()))
		for c in design.run.conditions:
			series_df.insert(series_df.shape[1],c+"^log2_stdev_lb",series_df[c+"^log2_average"].subtract(series_df[c+"^log2_stdev"].get_values()))
		for c in design.run.conditions:
			series_df.insert(series_df.shape[1],c+"^log2_stdev_ub",series_df[c+"^log2_average"].add(series_df[c+"^log2_stdev"].get_values()))





	
	else:
		n_missing = [c+"^rep_observation_count" for c in design.run.conditions]
		errcols = [c+"^stderr" for c in design.run.ratios]
		ratiogroups = {}
		log2groups = {}
		for ratio in design.run.ratios:
			divs = ratio.rstrip().split("|")
			numerator = divs[0]
			for denom in divs[1:]:
				replicates = []
				log2cols = []
				for sg in design.run.supergroups:
					divisors = [design.run.suffix+design.run.supergroups[sg].repmap[denom].name, design.run.suffix+design.run.supergroups[sg].repmap[numerator].name]
					series_df.insert(series_df.shape[1],ratio+"_"+design.run.supergroups[sg].name,series_df[divisors[1]].divide(series_df[divisors[0]]))
					replicates.append(ratio+"_"+design.run.supergroups[sg].name)
					series_df.insert(series_df.shape[1],ratio+"_"+design.run.supergroups[sg].name+"^log2vals",series_df[divisors[1]].divide(series_df[divisors[0]]).apply(np.log2).get_values())
					log2cols.append(ratio+"_"+design.run.supergroups[sg].name+"^log2vals")
			ratiogroups[ratio]=replicates
			log2groups[ratio]=log2cols

			for rg in ratiogroups:
				series_df.insert(series_df.shape[1],rg+"^repave",(series_df[[rep for rep in ratiogroups[rg]]]).mean(axis=1))
		
			series_df.insert(series_df.shape[1],rg+"^log2_average",(series_df[[c for c in log2cols]]).mean(axis=1))
			series_df.insert(series_df.shape[1],rg+"^log2_stdev",(series_df[[c for c in log2cols]]).std(axis=1))
			series_df.insert(series_df.shape[1],rg+"^log2_cov",series_df[rg+"^log2_average"].subtract(series_df[rg+"^log2_stdev"]).multiply(series_df[rg+"^log2_average"]).apply(np.abs))
			series_df.insert(series_df.shape[1],rg+"^log2_stdev_lb",(series_df[[c for c in log2cols]]).std(axis=1))
			series_df.insert(series_df.shape[1],rg+"^log2_lb",series_df[rg+"^log2_average"].subtract(series_df[rg+"^log2_stdev"]).get_values())
			series_df.insert(series_df.shape[1],rg+"^log2_ub",series_df[rg+"^log2_average"].add(series_df[rg+"^log2_stdev"]).get_values())

	return series_df,log2groups,n_missing,errcols











''' 20160627
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
		


def columns_normalize(df,cols):
	for col in cols:
		coldata = np.array(df[col].get_values().astype(float))
		##..unity normalize
		if coldata.min()==coldata.max():
			continue
		normalized = np.divide(np.subtract(coldata,coldata.min()),(coldata.max()-coldata.min()))
		df[col] = normalized
	return df

def dfappend_phos(df,series_df,cols):

	minidf = df[df['id'].isin(series_df['id'].values)]
	for col in cols:
		series_df.insert(series_df.shape[1],col,[v.split(";")[0] for v in minidf[col].astype(str).tolist()])
	return series_df
'''
		




def dfappend(df,series_df,cols):

	#append data from the untouched dataframe	

	df['id'] = df['id'].values.astype(str)
	series_df['id'] = series_df['id'].values.astype(str)

	errcols = [];
	errs_caught = False
	for col in cols:
		try:
   			vals = df[col]
		except KeyError:
			errs_caught=True
			errcols.append(col)
	[cols.remove(col) for col in errcols]

	ids = series_df['id']

	for col in cols:
		series_df.insert(series_df.shape[1],col,df[col][ids].values.astype(str))

	if errs_caught:	
		log("AppendWarning:: Could not append some column(s), as they are not present in original dataframe."+" :: "+str(errcols),2)
	else:
		log(str(cols) + " Appended successfully",1)


	return series_df















class Imputation:

	def medianFillImputation(df,cols):
		data = df[cols].get_values()
		imp = Imputer(missing_values='NaN', strategy='median', axis=1, verbose=0, copy=True)
		imputed = imp.fit_transform(data)
		df[cols] = imputed
		return df


	def return_non_missing_cols(row):
		indices = []
		for idx in range(row.shape[0]):
			if not np.isnan(row[idx]):
				indices.append(idx)
			return indices



	def adv_impute(df,headers):
		gbr_model = GradientBoostingRegressor(loss='ls', learning_rate=0.05, n_estimators=500, subsample=1.0, min_samples_split=2, min_samples_leaf=1, 		min_weight_fraction_leaf=0.0, max_depth=3, init=None, random_state=None, max_features=None, alpha=0.9, verbose=0, max_leaf_nodes=None, warm_start=False, presort='auto')  ##model for imputation

		data = df[headers].get_values()
		neigh = NearestNeighbors(100, 1.0, algorithm='auto',metric='euclidean') ##nn model
		imp_val = np.nan

	
		for col in range(data.shape[1]):
	
			row_num=0
			print "imputing column "+str(headers[col])
			for row in data:
				if np.isnan(row[col]):
					col_indices = return_non_missing_cols(row)	
					col_indices_plus = np.append(col,col_indices)
					#get all samples for col indices where there are no nans
					popblock = data[:,col_indices_plus][~np.isnan(data[:,col_indices_plus]).any(axis=1)]
					fullblock = data[:,:][~np.isnan(data[:,col_indices_plus]).any(axis=1)]				

					#get n nearest neighbors 
					model_fit = neigh.fit(popblock[:,1:popblock.shape[1]]) 
					knn_target = row[col_indices]
					nbrs_indices = model_fit.kneighbors(knn_target, return_distance=False)	
			

					tset = fullblock[nbrs_indices[0],:]				
					tset = tset[:,col_indices_plus]
				
					gbr_target_set = tset[:,0]
					gbr_training_set = tset[:,1:tset.shape[1]]

					normsum=0
					normfact = np.sum(np.arange(1,gbr_training_set.shape[0]+1,1))
					weights = [val/float(normfact) for val in reversed(np.arange(1,gbr_training_set.shape[0]+1,1))]
					gbr_model.fit(gbr_training_set,gbr_target_set,sample_weight=weights)
					row = row[col_indices]

					data[row_num,col] = gbr_model.predict(row)
					row_num+=1
				else:
					row_num+=1
			df[headers[col]] = data[:,col]
		return df












class StatisticalTest:
	@staticmethod
	def paired_ttest(self,df,C1,C2):
		G1vals = df[C1].get_values()
		G2vals = df[C2].get_values()
		pvals=[]
		stat=[]
		for ridx in range(G1vals.shape[0]):
			vals = ttest_ind(G1vals[ridx,:],G2vals[ridx,:])
			pvals.append(vals[1])
			stat.append(vals[0])
		df.insert(df.shape[1],'ttest_pvals',pvals)
		df.insert(df.shape[1],'ttest_stat',stat)
		return df


	def normcdf(x, mu, sigma):
		cdf = 0.5 * (1 + math.erf((x - mu)/math.sqrt(2 * sigma**2)))
		return min(cdf,1-cdf)



	def calculate_BH_correction(mapping): #apply benjamini-hochberg correction to values in a dict
		ordered_mapping = OrderedDict(sorted(mapping.items(),key=lambda t: t[1]))
		cntr=1
		output = []
		for key in ordered_mapping:
			mapping[key] = cntr*bh_FDR_cutoff/len(mapping)
			cntr+=1
		return mapping



	def BH_test(df,pvalcol,bhvalcol):
		passedlist = []
		df[pvalcol+'_passed_fdr'] = [False for x in range(df.shape[0])]
		for val in range(len(df[pvalcol].values)):
			if df[pvalcol].values[val]<df[bhvalcol].values[val]:
				passedlist.append(df[pvalcol].values[val])
	
		if len(sorted(passedlist))>0:
			maxval = sorted(passedlist)[-1]
		else:
			maxval = min(df[pvalcol].values)

	
		df[pvalcol+'_passed_fdr'][df[pvalcol]<maxval] = "True"
		return df



	@staticmethod
	def wmann_paired(self,df,c1,c2):
		log('Calculating Wilcox-Mann p-values via paired test',1)
		c1Vals = df[c1].get_values()
		c2Vals = df[c2].get_values()
		stats=[];pvals=[]
		for row in range(c2Vals.shape[0]):
			stat,pval = wilcoxon(c1Vals[row,:],c2Vals[row,:])
			stats.append(stat)
			pvals.append(pval)
		df.insert(df.shape[1],"^paired_U_test_stat",stats)
		df.insert(df.shape[1],"^paired_wmannpvals",pvals)
		return df

	




		'''
			bh_corrected = calculate_BH_correction(dict(zip(df.index.values,df[g+"^wmannpvals"].values)))
			corrected_vals = []
			for k in df.index.values:
				corrected_vals.append(bh_corrected[k])

			df.insert(df.shape[1],g+"^benj_hoch",corrected_vals)
			df = BH_test(df,g+"^wmannpvals",g+"^benj_hoch")
			pvalcols.append(g+"^wmannpvals")
			cntr+=1
		log("P-values calculated",1)



		if design.run.multitest:
			combined=[]
			for row in df[pvalcols].get_values():
				combined.append(combine_pvalues(row,method='fisher', weights=None)[1])
			log("Fisher combined p_value test completed",1)
			df.insert(len(df.columns),'fisher_combined',combined)
			bh_corrected = calculate_BH_correction(dict(zip(df.index.values,df['fisher_combined'].values)))
			corrected_vals = []
			for k in df.index.values:
				corrected_vals.append(bh_corrected[k])

			df.insert(df.shape[1],'benj_hoch',corrected_vals)
			df = BH_test(df,'fisher_combined','benj_hoch')
			log("Benjamini-hochberg FDR complete",1)
		'''
		



class Normalization:
	@staticmethod
	def quantileNormalize(df_input,columns):
	    df = df_input[columns].copy()
	    #compute rank
	    dic = {}
	    for col in df:
		dic.update({col : sorted(df[col])})
	    sorted_df = pd.DataFrame(dic)
	    rank = sorted_df.mean(axis = 1).tolist()
	    #sort
	    for col in df:
		t = np.searchsorted(np.sort(df[col]), df[col])
		df[col] = [rank[i] for i in t]
	    for col in df:
		df_input[col] = df[col]
	    return df_input














class DifferentialAnalysis:

	def __init__(self,mqrun_path,project_outpath):
		self.mqrun_path = mqrun_path
		self.project_outpath = project_outpath
		design.run.c2repnames()

	def get_runpath(self):
		return mqrun_path

	def get_outpath(self):
		return project_outpath

	def filter_replicates_by_group(self,df,n_reps,ctype):
		log("Filtering Missing Replicate Measurements by Condition",1)
		'''
		Allow for n missing replicates inside a condition
		'''
		if design.run.runtype=="Silac" and ctype =='collapse_zeros':
			ctype = 'collapse_nans'
			log("Converting to 'collapse nans' filter for SILAC run..",1)
		for c in design.run.conditions:
			cols = design.run.c2RNames[c]
			df = apply_filter(ctype,n_reps,cols,df)
		return df


	def filterZerosToNans(self,df,n_reps,ctype):
		log("Converting Zeros to Nans",1)
		for c in design.run.conditions:
			cols = design.run.c2RNames[c]
			df = apply_filter(ctype,None,cols,df)
		return df


	def filter_absent_conditions(self,df,n_conds,ctype):
		log("Filtering Absent Conditions",1)	
		zmap = np.zeros(df.shape[0],dtype=int)
		for c in design.run.conditions:
			groupheaders = design.run.c2RNames[c]
			missing_sum = df[groupheaders].isnull().sum(axis=1).values
			missing_sum[missing_sum < len(groupheaders)] = 0
			missing_sum[missing_sum != 0] = 1
			zmap = np.add(missing_sum,zmap)
		df.insert(df.shape[1],'sum_missing_conditions',zmap)
		df = apply_filter(ctype,n_conds,['sum_missing_conditions'],df)
		return df		



	def combine_treps(self,df):
		log("Combining Technical Replicates",1)
		for c in design.run.conditions:
			for rep in design.run.condition_map[c].replicates:
				if len(rep.treps)>0:
					trepheaders = [design.run.suffix+k for k in rep.treps.keys()]
					df[design.run.suffix+rep.name] = df[trepheaders].mean(axis=1)
		return df


	def apply_all_filters(self,df):
		for f in design.run.filters:
			thresh = f.threshold
			if f.threshold!='None':
				thresh=float(f.threshold)
			if f.subtype=='global':
				df = apply_filter(f.ftype,thresh,f.headers,df)
			if f.subtype=='conditional':
				m= getattr(self,f.ftype)
				df = m(df,thresh,f.function)	

		return df

	def normalize(self,df):
		N = Normalization
		if design.run.norm is None:
			return df
		scope = design.run.norm.scope
		t = design.run.norm.type
		if t=="quantile":
			log("Performing Quantile Normalization",1)
			for c in design.run.conditions:
				cols = design.run.c2RNames[c]
				df = N.quantileNormalize(df,cols)
		return df


	
	def runStatTests(self,df):
		Test = StatisticalTest
		for test in design.run.tests:
			func = getattr(Test,test.ttype)
			df = func(Test, df, design.run.c2RNames[test.c1],design.run.c2RNames[test.c2])
		return df



	def expand_site_table(self,df):
	    """	
		Method borrowed from PyMaxQuant 
	    """
	    log("Expanding Site Table.",1)
	    df = df.copy()

	    idx = df.index.names
	    df.reset_index(inplace=True)

	    def strip_multiplicity(df):
		df.columns = [c[:-4] for c in df.columns]
		return df
		
	    def strip_multiple(s):
		for sr in ['___1','___2','___3']:
		    if s.endswith(sr):
		        s = s[:-4]
		return s

	    base = df.filter(regex='.*(?<!___\d)$')
	    
	    # Remove columns that will match ripped multiplicity columns
	    for c in df.columns.values:
		if strip_multiple(c) != c and strip_multiple(c) in list(base.columns.values):
		    base.drop(strip_multiple(c), axis=1, inplace=True)

	    multi1 = df.filter(regex='^.*___1$')
	    multi1 = strip_multiplicity(multi1)
	    multi1['Multiplicity'] = '___1'
	    multi1 = pd.concat([multi1, base], axis=1)

	    multi2 = df.filter(regex='^.*___2$')
	    multi2 = strip_multiplicity(multi2)
	    multi2['Multiplicity'] = '___2'
	    multi2 = pd.concat([multi2, base], axis=1)

	    multi3 = df.filter(regex='^.*___3$')
	    multi3 = strip_multiplicity(multi3)
	    multi3['Multiplicity'] = '___3'
	    multi3 = pd.concat([multi3, base], axis=1)

	    df = pd.concat([multi1, multi2, multi3], axis=0)
	    df['id'] = ["%s%s" % (a, b) for a, b in zip(df['id'], df['Multiplicity'])]

	    if idx[0] is not None:
		df.set_index(idx, inplace=True)

	    df['multiplicity'] = df['id'].str.split("___").str.get(1)
	    df['id'] = df['id'].str.split("___").str.get(0)
	    return df

			

	
class ModificationAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis):
		log("Running ModSite differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(self.get_runpath()+design.run.runfile, sep = "\t")  
		df = self.expand_site_table(df)
		df = self.apply_all_filters(df)
		df = self.combine_treps(df)
		df = self.normalize(df)
		series_df,log2groups,missing_count,errcols = generate_series('mod',df)
		series_df = self.runStatTests(series_df)

		if design.run.appendGenes == "True":
			series_df = UniQuery("gene_query").queryUniprotGenes(series_df,'leading_accession')		

		series_df.to_csv(self.get_outpath()+"oxidation_diffexp_mf.txt",sep="\t",index=False)
		log("____ModSite differential analysis complete____\n",1)



class ProteinAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis):
		log("Running Protein differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(self.get_runpath()+design.run.runfile, sep = "\t") 
		df.rename(columns={'Majority protein IDs': 'Leading proteins'}, inplace=True) 
		df = self.apply_all_filters(df)
		df = self.combine_treps(df)
		df = self.normalize(df)
		series_df,log2groups,missing_count,errcols = generate_series('prot',df)
		series_df = self.runStatTests(series_df)

		if design.run.appendGenes == "True":
			series_df = UniQuery("gene_query").queryUniprotGenes(series_df,'leading_accession')	

		series_df.to_csv(self.get_outpath()+"protein_diffexp.txt",sep="\t",index=False)
		log("____Protein differential analysis complete____\n",1)


class PeptideAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis):
		log("Running Peptide differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(self.get_runpath()+design.run.runfile, sep = "\t") 
		df.rename(columns={'Proteins': 'Leading proteins'}, inplace=True) 
		df = self.apply_all_filters(df)
		df = self.combine_treps(df)
		df = self.normalize(df)
		series_df,log2groups,missing_count,errcols = generate_series('pep',df)
		series_df = self.runStatTests(series_df)
		
		#if design.run.appendGenes == "True":
		#	series_df = UniQuery("gene_query").queryUniprotGenes(series_df,'leading_accession')	

		series_df.to_csv(self.get_outpath()+"peptide_diffexp.txt",sep="\t",index=False)
		log("____Peptide differential analysis complete____\n",1)

	





if __name__ == "__main__":
	mqrun_path = './'
	project_outpath = mqrun_path+"report/"
	if not os.path.exists(project_outpath):
		os.makedirs(project_outpath)
	

	design = DesignParser(mqrun_path+"EGF_peptides.xml")
	analysis = DifferentialAnalysis(mqrun_path,project_outpath)
	#ModificationAnalysis(analysis)
	#ProteinAnalysis(analysis)
	PeptideAnalysis(analysis)

