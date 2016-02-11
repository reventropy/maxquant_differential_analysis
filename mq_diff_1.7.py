import warnings
warnings.filterwarnings('ignore')

from scipy.stats import ranksums
from scipy.stats import sem
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np
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
import os
import seaborn
from collections import OrderedDict


########################## setup experimental design
run_name = "Taxol Tymadine Run"
silac = False
mqrun_path = './dyrk1a/most_recent/'
project_outpath = mqrun_path+"report/"
exp_design_path = "./dyrk1a_expdesign.xml"
##########################








class DesignParser:
	def __init__(self,path):
		self.design_path = path
		self.parse(path)
	def parse(self,path):
		tree = ET.parse(path)
		root = tree.getroot()
		self.run = Run(root.get('Description')) #instantiate run object and append description
		


		for cond in root.iter('Condition'):  #initiate all condition objects
			cond_name = cond.attrib.get('name')
			cond_type = cond.attrib.get('type')
			cond_desc = cond.attrib.get('description')
			condition = self.run.condition_map[cond_name] = Condition(cond_name,cond_type,cond_desc)
			self.run.conditions.append(cond_name)
			for rg in cond.iter('ReplicateGroup'):  ## initiate replicate group objects
				rg_name = rg.attrib.get('name')
				rg_type = rg.attrib.get('type')
				rg_desc = rg.attrib.get('description')
				replicate_group = condition.replicate_groups[rg_name] = ReplicateGroup(rg_name,rg_type,rg_desc)
				self.run.groups.append(rg_name)
				for rep in rg.iter('Replicate'):
					rep_name = rep.attrib.get('name')
					rep_type = rep.attrib.get('type')
					rep_desc = rep.attrib.get('description')
					replicate_group.replicates[rep_name] = Replicate(rep_name,rep_type,rep_desc)
					self.run.replicates.append(rep_name)
					self.run.condition_map[cond_name].replicates.append(rep_name)

		for ratios in root.iter("Ratios"):
			for ratio in ratios.iter("Ratio"):
				self.run.ratios.append(ratio.attrib.get("des"))





class Run:
	def __init__(self,description):
		self.description = description
		self.condition_map = {} #head node for run traversal
		self.ratios = []
		self.conditions  = []
		self.groups = []
		self.replicates = []


	def print_design(self,):
		for condition in self.condition_map:
			cc = self.condition_map[condition]
			print cc.attributes

			for rg in cc.replicate_groups:
				group = cc.replicate_groups[rg]
				print '\t', group.attributes
				
				for rep in group.replicates:
					print '\t\t', group.replicates[rep].attributes




class Xml_obj():
	def build(self,name,Ttype,description):
		self.name = name
		self.type = Ttype
		self.description = description
		self.attributes = {'name':self.name,'type':self.type,'desc':self.description}


class Condition(Xml_obj):
	def __init__(self,name,Ttype,description):
		self.build(name,Ttype,description)
		self.replicate_groups = {}
		self.replicates = []
	
class ReplicateGroup(Xml_obj):
	def __init__(self,name,Ttype,description):
		self.build(name,Ttype,description)
		self.replicates = {}

class Replicate(Xml_obj):
	def __init__(self,name,Ttype,description):
		self.build(name,Ttype,description)

		
		


















class Report:
	def __init__(self,name):
		self.name = name
		self.subreports = {}

	def append_subreport(self,subreport,subreport_name):
		self.subreports[subreport_name] = subreport



class SubReport(Report):
	def __init__(self,repo_name):
		self.repo_name = repo_name
		self.figure = plt.figure()
		self.cnt=1
	def plot_histo(self,df,cols):

		cntr=1
		vals = df[cols].get_values()
		fig = self.figure
		ax = fig.add_subplot(1,2,self.cnt)
		self.cnt+=1
		ax.axes.get_xaxis().set_visible(False)
		ax.axes.get_yaxis().set_visible(False)
		ax.axes.grid(False)
		for col in cols:
			ax = fig.add_subplot(4,4,cntr); cntr+=1
			nbins = int(len(df[col].values)/10)
			ax.hist(df[col].values,bins=nbins)
			ax.set_title(col+" bins:["+str(nbins)+"]")
			ax.set_xlim([np.min(vals),np.max(vals)])
		fig.tight_layout()
		
	def plot_cor(self,df,cols):
		cntr=1
		vals = df[cols].get_values()
		cormat = np.corrcoef(vals.T)
		#for i in range(len(cols)):
		#	for j in range(len(cols)):
		#		cormat[i,j] = pearsonr(df[cols[i]].values,df[cols[j]].values)[1]
		#		print '%d , %d , %d' % (i,j,cormat[i,j])

		column_labels = cols
		row_labels = cols
		fig, ax = plt.subplots()
		heatmap = ax.pcolor(cormat, cmap=plt.cm.Blues)
		ax.set_xticks(np.arange(cormat.shape[0])+0.5, minor=False)
		ax.set_yticks(np.arange(cormat.shape[1])+0.5, minor=False)
		ax.invert_yaxis()
		ax.xaxis.tick_top()
		ax.set_xticklabels(row_labels, minor=False)
		ax.set_yticklabels(column_labels, minor=False)
		plt.show()
			

	def seal(self,):
		plt.show()
		#plt.savefig(project_outpath+'histo.pdf')









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


'''
def calculate_mean(df,array_of_headers,resulting_col_name):
	print str(array_of_headers)
	if header_suffix+" "+resulting_col_name+"_mean" not in df.columns:
		df.insert(df.shape[1],header_suffix+" "+resulting_col_name+"_mean",df[array_of_headers].mean(axis=1))
	return df


def calculate_std(df,array_of_headers,resulting_col_name):
	if header_suffix+" "+resulting_col_name+"_err" not in df.columns:
		df.insert(df.shape[1],header_suffix+" "+resulting_col_name+"_err",df[array_of_headers].std(axis=1))
	return df

'''







def impute(df,imp_val,headers):
	if np.isnan(imp_val):
		imp_val = -500	
	log("imputing...",1)

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
	log(str(columns)+" Converted to Log_2 space",1)
	return df









def append_pval(df,columns,multitest):
	pvalcols = []
	for col in columns:
		pvals = [normcdf(val,np.mean(df[col].values),np.std(df[col].values)) for val in df[col].values]
		df.insert(len(df.columns),col+"_pvals",pvals)
		pvalcols.append(col+"_pvals")
	log("P-values calculated for following groups: "+str(columns),1)
	if multitest:
		combined=[]
		for row in df[pvalcols].get_values():
			combined.append(combine_pvalues(row,method='fisher', weights=None)[1])
		log("Fisher combined p_value test completed",1)
		df.insert(len(df.columns),'fisher_combined_pval',combined)
		bh_corrected = bh_correct(dict(zip(df.index.values,df['fisher_combined_pval'].values)))
		corrected_vals = []
		for k in df.index.values:
			corrected_vals.append(bh_corrected[k])
		df.insert(len(df.columns),'benj_hoch_corrected_pval',corrected_vals)
		df.drop(pvalcols,axis=1,inplace=True)
		log("Benjamini-hochberg correction successfully applied to combined p-values",1)
	return df




def append_wilcoxmann(df,columns,multitest):
	pvalcols = []
	groups = [ratcol+"_ratio" for ratcol in design.run.ratios]
	cntr=0
	for col in columns:
		pvals = []
		data = df[col].values
		for vals in data:
			pvals.append(ranksums([v for v in vals],[item for sublist in data for item in sublist])[1])
		df.insert(len(df.columns),groups[cntr].replace("_ratio","")+"^wmannpvals",pvals)
		pvalcols.append(groups[cntr].replace("_ratio","")+"^wmannpvals")
		cntr+=1
	log("P-values calculated for following groups: "+str(columns),1)

	if multitest:
		combined=[]
		for row in df[pvalcols].get_values():
			combined.append(combine_pvalues(row,method='fisher', weights=None)[1])
		log("Fisher combined p_value test completed",1)
		df.insert(len(df.columns),'fisher_combined_wmannpval',combined)
		bh_corrected = bh_correct(dict(zip(df.index.values,df['fisher_combined_wmannpval'].values)))
		corrected_vals = []
		for k in df.index.values:
			corrected_vals.append(bh_corrected[k])
		df.insert(len(df.columns),'benj_hoch_corrected_wmannpval',corrected_vals)
		log("Benjamini-hochberg correction successfully applied to combined p-values",1)

	return df
	





''' deprecated
def float_convert(df,columns):
	print "converting"
	for col in columns:
		if df[col].dtypes == 'int64':
			df[col] = df[col].astype(float)
			print col,df[col].dtypes
	return df
'''




#this applies a filter to the entire dataset based on the headers argument
def apply_filter(ftype,f_threshold,headers,df): ## all filters are inclusive
	
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
		df[headers[0]] = df[headers[0]].apply(str)
		df = df[~df[headers[0]].str.contains("CON__")]
	if ftype == "reverse":
		df[headers[0]] = df[headers[0]].apply(str)
		df = df[~df[headers[0]].str.contains("REV__")]
	if ftype == 'collapse_nans':
		df = df[df[headers].isnull().sum(axis=1)<=f_threshold]
	if ftype == 'collapse_zeros':
		df = df[(df[headers]==0.0).sum(axis=1)<=f_threshold]
	return df
	
		









def generate_series(analysis_type,df,suffix):
	
	group=[]	
	ids = ['id']

	if analysis_type=='phospho':
		group = ['Leading proteins','Amino acid','Positions','Number of Phospho (STY)','Phospho (STY) Probabilities']
	else:
		group = ['Majority protein IDs']

		

	raw_ratios_headers = []
	raw_ratio_repgroups = []
	l2fc_ratios = []
	ratios = [ratcol+"^ratio" for ratcol in design.run.ratios]
	series_df = pd.DataFrame()


	[series_df.insert(series_df.shape[1],col,df[col].values) for col in ids+group]
	'''
	insert raw measurements
	'''
	[series_df.insert(series_df.shape[1],suffix+rep,df[suffix+rep].values) for c in design.run.conditions for rep in design.run.condition_map[c].replicates]
	'''
	count up number of missing values from each condition and append a column
	'''
	n_missing = [c+"^missing_count" for c in design.run.conditions]
	for c in design.run.conditions:
		series_df.insert(series_df.shape[1],c+"^missing_count",(df[[suffix+rep for rep in design.run.condition_map[c].replicates]].values==0).astype(int).sum(axis=1))
	'''
	Start by taking ratios for each replicate across all conditions.  if silac, then invert (1/x), values in relevant columns.
	'''
	if silac:
		for exp in design.run.conditions:
			for rep in range(len(exp)):
				if inverse[rep] == True:
					df[suffix+" "+exp[rep]] = 1/df[suffix+" "+exp[rep]].values		
	'''
	generate a ratio column for each replicate
	'''
	l2fcmap = OrderedDict({}); raw_ratio_map = OrderedDict({})

	for ratio in design.run.ratios: 
		numerator = ratio.split('|')[0];
		for denom in range(1,len(ratio.split("|"))):
			z = zip(design.run.condition_map[numerator].replicates, design.run.condition_map[ratio.split("|")[denom]].replicates)
			raw_ratios = [('|').join(x) for x in z]
			ratio_group = []
			l2fc_ratio_group = []
			for ratio in raw_ratios:
				ratio_arr = ratio.split('|')
				l2fcmap[ratio+'^raw ratio'] = np.divide(df[suffix+ratio_arr[0]].values.astype(float),df[suffix+ratio_arr[1]].values.astype(float))
				raw_ratio_map[ratio+'^l2fc ratio'] = np.log2(np.divide(df[suffix+ratio_arr[0]].values.astype(float),df[suffix+ratio_arr[1]].values.astype(float)))
				l2fc_ratio_group.append(ratio+'^l2fc ratio')
				raw_ratios_headers.append(ratio+'^raw ratio')
				ratio_group.append(ratio+'^raw ratio')
			l2fc_ratios.append(l2fc_ratio_group)
			raw_ratio_repgroups.append(ratio_group)



	'''
	Append columns just produced to dataframe
	'''
	[series_df.insert(series_df.shape[1],key,val) for key,val in l2fcmap.iteritems()]
	[series_df.insert(series_df.shape[1],key,val) for key,val in raw_ratio_map.iteritems()]


	series_df.replace([np.inf, -np.inf], np.nan,inplace=True) ##filter entire df
	series_df =  apply_filter('impute_nan',None,[c for outer in l2fc_ratios for c in outer],series_df)	

	'''
	output the mean, errors, l2fc for each replicate ratio
	'''
	[series_df.insert(series_df.shape[1],ratios[replicates].split('^')[0]+'^raw ratio mean',series_df[raw_ratio_repgroups[replicates]].mean(axis=1, skipna=True).values) for replicates in range(len(raw_ratio_repgroups))]
	[series_df.insert(series_df.shape[1],ratios[replicates].split('^')[0]+'^raw ratio stderr',series_df[raw_ratio_repgroups[replicates]].std(axis=1, skipna=True).values) for replicates in range(len(raw_ratio_repgroups))]
	[series_df.insert(series_df.shape[1],ratios[replicates].split('^')[0]+'^l2fc mean',series_df[l2fc_ratios[replicates]].mean(axis=1, skipna=True)) for replicates in range(len(l2fc_ratios))]
	[series_df.insert(series_df.shape[1],ratios[replicates].split('^')[0]+'^l2fc stderr',series_df[l2fc_ratios[replicates]].apply(np.nanstd,axis=1)) for replicates in range(len(l2fc_ratios))]
	[series_df.insert(series_df.shape[1],ratios[replicates].split('^')[0]+'^perc_err',series_df[l2fc_ratios[replicates]].apply(np.nanstd,axis=1).div(series_df[l2fc_ratios[replicates]].mean(axis=1, skipna=True), axis='index').apply(abs)) for replicates in range(len(l2fc_ratios))]
	[series_df.insert(series_df.shape[1],ratios[replicates].split('^')[0]+'^flagged',(series_df[ratios[replicates].split('^')[0]+'^perc_err']>0.5)) for replicates in range(len(l2fc_ratios))]
	errcols = [ratios[replicates].split('^')[0]+'^raw ratio stderr' for replicates in range(len(raw_ratio_repgroups))]

	

	return series_df,l2fc_ratios,n_missing,errcols












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
		



		


def generate_phospho_path_inp(df):
	nodes = set()
	time_map = {}
	accessions = df['uniprot_accession'].values
	amino_acids = df['amino_acid'].values
	sites = df['phospho_site'].values
	multiplicity = df['multiplicity'].values
	tcounter=1

	for g in group_names:
		expression = df[g].values
		for i in range(len(accessions)):
			acc=accessions[i].split(";")[0].split("-")[0]
			aa=amino_acids[i]
			s = sites[i].split(";")[0]
			line = acc+"\t"+acc+'-'+aa+s+"\t"+aa+s+"\t"+multiplicity[i]+"\n"
			time_map[acc+"-"+aa+s+"-"+multiplicity[i]+"-"+str(tcounter)] = expression[i] 
			nodes.add(line)
		tcounter+=1

	fout = open(project_outpath+'phospho_path_nodes.phos','w')
	for val in nodes:
		fout.write(val)
	fout.close()
	dfout = pd.DataFrame(columns = ['ID', 'Ratio'])
	dfout['ID'] = time_map.keys()
	dfout['Ratio'] = time_map.values()
	dfout.to_csv(project_outpath+'phospho_path_times.txt',sep="\t",index=False)







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
	if ids.str.split("__").shape[0]>0: ##then it's phospho data
		ids = [int(i.split("__")[0]) for i in series_df['id'].values]
		for col in cols:
			series_df.insert(series_df.shape[1],col,df[col][ids].values.astype(str))
	else:
		for col in cols:
			series_df.insert(series_df.shape[1],col,df[col][ids].values.astype(str))

	if errs_caught:	
		log("AppendWarning:: Could not append some column(s), as they are not present in original dataframe."+" :: "+str(errcols),2)
	else:
		log(str(cols) + " Appended successfully",1)


	return series_df




def output_for_chiayu_phospho(df):
	dfout = pd.DataFrame(columns=['gene_name','site','course'])
	df['Gene names'] = df['Gene names'].astype(str)
	genes = df['Gene names'].values
	aa = df['amino_acid'].values
	position = df['phospho_site'].values
	tcourse = df[group_names].values
	arrs_out = [[] for i in range(3)]	
		

	for idx in range(len(genes)):

		arrs_out[0].append(genes[idx].split(';')[0])
		arrs_out[1].append(aa[idx]+str(position[idx].split(';')[0]))
		tc_str = ""
		for timeval in tcourse[idx]:
			tc_str += str(timeval)+","
		arrs_out[2].append(tc_str[0:len(tc_str)-1])

	dfout['gene_name'] = arrs_out[0]
	dfout['site'] = arrs_out[1]
	dfout['course'] = arrs_out[2]
	fout = open(project_outpath+"/for_cy_header.txt",'w')
	fout.write(str(group_names))
	fout.close()
	dfout.to_csv(project_outpath+"/for_cy_exp_phospho.txt",sep="\t",index=False)
 


def output_for_chiayu_protein(df):
	dfout = pd.DataFrame(columns=['gene_name','site','course'])
	df['Gene names'] = df['Gene names'].astype(str)
	genes = df['Gene names'].values
	tcourse = df[group_names].values
	arrs_out = [[] for i in range(3)]	
		

	for idx in range(len(genes)):

		arrs_out[0].append(genes[idx].split(';')[0])
		arrs_out[1].append("")
		tc_str = ""
		for timeval in tcourse[idx]:
			tc_str += str(timeval)+","
		arrs_out[2].append(tc_str[0:len(tc_str)-1])

	dfout['gene_name'] = arrs_out[0]
	dfout['site'] = arrs_out[1]
	dfout['course'] = arrs_out[2]
	fout = open(project_outpath+"/for_cy_header_prot.txt",'w')
	fout.write(str(group_names))
	fout.close()
	dfout.to_csv(project_outpath+"/for_cy_exp_prot.txt",sep="\t",index=False)






def expand_site_table(df):
    """	
	Method borrowed from PyMaxQuant 
    """

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

    return df




		








class DifferentialAnalysis:

	def __init__(self,run_name):
		self.run_name = run_name
		self.report = Report(run_name)






	def filter_replicates_by_group(self,df,n_reps,ctype,suffix):
		'''
		Allow for n replicates inside a replicate group to be missing.  This is a more stringent, and specific filter
		'''
		conditions = design.run.condition_map
		for c in conditions:
			repgroups = conditions[c].replicate_groups
			for rg in repgroups:
				reps = repgroups[rg].replicates
				passed = True
				cols = [suffix+c for c in reps.keys()]
				df = apply_filter(ctype,n_reps,cols,df)
		return df





	def filter_replicates_by_condition(self,df,n_reps,target,ctype,suffix):
		'''
		Require at least n replicates inside a condition to be real numbers (not NANs)
		'''
		conditions = design.run.condition_map
		for c in conditions:
			repgroups = conditions[c].replicate_groups
			repheaders = []
			for rg in repgroups:
				reps = repgroups[rg].replicates
				passed = True
				cols = [suffix+c for c in reps.keys()]
				repheaders+=cols
			df = apply_filter(ctype,n_reps,repheaders,df)
		return df	
					
					
		
		
	def apply_global_filters(self,series_df,avcols,errcols):
		log('Applying global filters',1)
		series_df.replace([np.inf, -np.inf], np.nan,inplace=True) ##filter entire df
		series_df =  apply_filter('impute_nan',None,[c for outer in avcols for c in outer],series_df)  ##impute missing values for fold change
		series_df = apply_filter('impute_nan',None,errcols,series_df) ##impute missing values for standard error columns
		series_df = append_wilcoxmann(series_df,avcols,multitest=True)
		#series_df = apply_filter('lowpass',0.05,['benj_hoch_corrected_pval'],series_df)
		return series_df




	
class PhosphoAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis,suffix):
		log("Running phospho differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(mqrun_path+"Phospho (STY)Sites.txt", sep = "\t") 
		df = apply_filter('contaminants',None,['Leading proteins'],df) 
		df = apply_filter('reverse',None,['Leading proteins'],df) 
		df = expand_site_table(df)
		df = self.filter_replicates_by_group(df,1,'collapse_zeros',suffix)
		
		series_df,avcols,missing_count,errcols = generate_series('phospho',df,suffix)
		series_df = dfappend(df,series_df,['Gene names'])
		series_df = self.apply_global_filters(series_df,avcols,errcols)
		series_df.to_csv(project_outpath+"phosphosite_diffexp.txt",sep="\t",index=False)
		log("____Phosphopeptide differential analysis complete____\n",1)



class ProteinAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis,suffix):
		log("Running protein differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(mqrun_path+"proteinGroups.txt", sep = "\t")
		df = apply_filter('contaminants',None,['Majority protein IDs'],df) 
		df = apply_filter('reverse',None,['Majority protein IDs'],df) 
		df = self.filter_replicates_by_group(df,2,'collapse_zeros',suffix)
		series_df,avcols,missing_count,errcols = generate_series('protein',df,suffix)
		series_df = self.apply_global_filters(series_df,avcols,errcols)
		series_df = dfappend(df,series_df,['Gene names'])
		series_df.to_csv(project_outpath+"protein_diffexp.txt",sep="\t",index=False)
		log("____Protein differential analysis complete____\n",1)



	





if __name__ == "__main__":
	if not os.path.exists(project_outpath):
		os.makedirs(project_outpath)
	design = DesignParser(exp_design_path)
	analysis = DifferentialAnalysis(run_name)
	PhosphoAnalysis(analysis, 'Intensity ')
	ProteinAnalysis(analysis, 'iBAQ ')
