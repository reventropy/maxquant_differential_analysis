import warnings
warnings.filterwarnings('ignore')

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
import seaborn




########################## setup experimental design
#report_name = "CDK8_inhibited_run_phospho_enriched"
#silac = True
#mqrun_path = '/home/jjacob/python/zposs/Phospho_MQ_search/'
#project_outpath = mqrun_path+"report/"
#header_suffix = 'Ratio H/L normalized'
#exp_design = [['1','2','3']]
#inverse = [False,False,True]
#group_names = ['treatment_1h']
##########################begin differential analysis


########################## setup experimental design
#run_name = "CDK8_inhibited_run_proteome"
#silac = True
#mqrun_path = '/home/jjacob/python/zposs/Proteome/txt/'
#project_outpath = mqrun_path+"report/"
#header_suffix = 'Ratio H/L normalized'
#exp_design = [['1','7'],['2','8'],['3','9'],['4','10'],['5','11'],['6','12']]
#inverse = [False,True]
#group_names = ['DMSO','treatment_1h','treatment_3h','treatment_6h','treatment_18h','treatment_24h']
##########################

########################## setup experimental design
run_name = "Taxol Tymadine Run"
silac = False
mqrun_path = '/home/jjacob/python/tx293_synch/'
project_outpath = mqrun_path+"report/"
header_suffix = 'iBAQ'
biological_replicates = {}
conditions = {}

biological_replicates['DMSOM_brep1'] = ['DMSOM_brep1_trep1','DMSOM_brep1_trep2'] #DMSOM
biological_replicates['DMSOM_brep2'] = ['DMSOM_brep2_trep1','DMSOM_brep2_trep2'] 
biological_replicates['TaxolM_brep1'] = ['TaxolM_brep1_trep1','TaxolM_brep1_trep2'] #TaxolM
biological_replicates['TaxolM_brep2'] = ['TaxolM_brep2_trep1','TaxolM_brep2_trep2'] 
biological_replicates['TaxolG2_brep1'] = ['TaxolG2_brep1_trep1','TaxolG2_brep1_trep2'] #TaxolG2
biological_replicates['TaxolG2_brep2'] = ['TaxolG2_brep2_trep1','TaxolG2_brep2_trep2'] 
biological_replicates['DMSOG2_brep1'] = ['DMSOG2_brep1_trep1','DMSOG2_brep1_trep2'] #DMSOG2
biological_replicates['DMSOG2_brep2'] = ['DMSOG2_brep2_trep1','DMSOG2_brep2_trep2'] 
biological_replicates['Thymadine'] = ['Thymadine_Release_trep1','Thymadine_Release_trep2'] #Thymadine

conditions['DMSOM'] = ['DMSOM_brep1','DMSOM_brep2']
conditions['TaxolM'] = ['TaxolM_brep1','TaxolM_brep2']
conditions['TaxolG2'] = ['TaxolG2_brep1','TaxolG2_brep2']
conditions['DMSOG2'] = ['DMSOG2_brep1','DMSOG2_brep2']
conditions['Thymadine'] = ['Thymadine']


group_names = ['TaxolM|DMSOM|Thymadine','TaxolG2|DMSOG2|Thymadine']
##########################













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



def calculate_mean(df,array_of_headers,resulting_col_name):
	print str(array_of_headers)
	if header_suffix+" "+resulting_col_name+"_mean" not in df.columns:
		df.insert(df.shape[1],header_suffix+" "+resulting_col_name+"_mean",df[array_of_headers].mean(axis=1))
	return df


def calculate_std(df,array_of_headers,resulting_col_name):
	if header_suffix+" "+resulting_col_name+"_err" not in df.columns:
		df.insert(df.shape[1],header_suffix+" "+resulting_col_name+"_err",df[array_of_headers].std(axis=1))
	return df









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
		log("Benjamini-hochberg correction successfully applied to combined p-values",1)
	return df








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
			vals_no_nans = vals[~np.isnan(vals)]
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
	if ftype == 'collapse':
		df = df[df[headers].isnull().sum(axis=1)<=f_threshold]

	return df




def generate_protein_series(df):
	ids = ['id']
	pids = ['Majority protein IDs']	
	avcols = [gn+"_repave" for gn in group_names]
	errcols = [gn+"_err" for gn in group_names]
	cols = ids+pids+avcols+errcols
	series_df = pd.DataFrame(columns = cols)
	

	for exp in exp_design:
		for rep in range(len(exp)):
			if inverse[rep] == True:
				df[header_suffix+" "+exp[rep]] = 1/df[header_suffix+" "+exp[rep]].values

	for gn in range(len(group_names)):
		df[group_names[gn]+"_repave"] = df[[header_suffix+" "+rep for rep in exp_design[gn]]].mean(axis=1)
		df[group_names[gn]+"_err"] = df[[header_suffix+" "+rep for rep in exp_design[gn]]].std(axis=1)

	for col in cols:
		series_df[col] = df[col]
	
	return series_df,avcols










def generate_phospho_series(df): #phospho specific
	ids = ['id']
	group = ['Leading proteins']
	AAs = ['Amino acid']
	positions = ['Positions']
	multiplicity = ['Number of Phospho (STY)']
	avcols = [gn+"_repave" for gn in group_names]
	errcols = [gn+"_err" for gn in group_names]
	
	additional_cols = ids+group+AAs+positions+multiplicity
	sdf_cols = additional_cols+avcols+errcols
	series_df = pd.DataFrame(columns = sdf_cols)	

	#invert ratios where appropriate, and overwrite in original dataframe
	for exp in exp_design:
		for rep in range(len(exp)):
			if inverse[rep] == True:
				df[header_suffix+" "+exp[rep]] = 1/df[header_suffix+" "+exp[rep]].values
	#average across replicates and make appropriate calculations
	#for error.
	for gn in range(len(group_names)):
		series_df[group_names[gn]+"_repave"] = df[[header_suffix+" "+rep for rep in exp_design[gn]]].mean(axis=1)
		series_df[group_names[gn]+"_err"] = df[[header_suffix+" "+rep for rep in exp_design[gn]]].std(axis=1)
		
	for col in additional_cols:
		series_df[col] = df[col]
	return series_df,avcols












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
			

	def apply_global_filters(self,series_df,avcols):
		log('Applying global filters',1)
		series_df.replace([np.inf, -np.inf], np.nan,inplace=True)
		series_df = apply_filter('collapse',3,avcols,series_df) #only allow 1 points to be missing (from both replicates)
		series_df = apply_filter('impute_nan',None,avcols,series_df)	
		series_df = append_pval(series_df,avcols,multitest=True)
		series_df = apply_filter('lowpass',0.05,['benj_hoch_corrected_pval'],series_df)
		return series_df

	@staticmethod
	def generate_report(current_analysis,name,full_df,series_df,avcols):
		subreport = SubReport(name)
		#subreport.plot_histo(series_df,avcols)
		raw_columns = raw_columns = [header_suffix+" "+x for x in [e for exp in exp_design for e in exp]]
		subreport.plot_cor(apply_filter('remove_nan',None,raw_columns,full_df) ,raw_columns)
		current_analysis.report.append_subreport(name,subreport)
		subreport.seal()
	header_suffix+" "


	def combine_replicates(self,df,arr,appended_str):
		print 'here'
		for rep in arr:	
			df= calculate_mean(df,[header_suffix+" "+r+appended_str for r in arr[rep]],rep)
			df= calculate_std(df,[header_suffix+" "+r+appended_str for r in arr[rep]],rep)
		return df


class PhosphoAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis):
		log("Running phospho differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(mqrun_path+"Phospho (STY)Sites.txt", sep = "\t")
		df = apply_filter('remove_nan',None,['Number of Phospho (STY)'],df) 
		df = apply_filter('contaminants',None,['Leading proteins'],df) 
		df = apply_filter('reverse',None,['Leading proteins'],df) 
		df_expanded = expand_site_table(df)
		df_expanded.to_csv(project_outpath+"phosphosite_sty_appended.txt",sep="\t",index=False)
		series_df,avcols = generate_phospho_series(df_expanded)
		series_df = dfappend(df,series_df,['Gene names'])
		series_df = dfappend(df,series_df,['Phospho (STY) Probabilities'])
		series_df = self.apply_global_filters(series_df,avcols)
		series_df.to_csv(project_outpath+"phosphosite_diffexp.txt",sep="\t",index=False)
		subreport = self.generate_report(current_analysis,"Phospho Level Analysis Report",df,series_df,avcols)
		log("____Phosphopeptide differential analysis complete____\n",1)



class ProteinAnalysis(DifferentialAnalysis):
	def __init__(self,current_analysis):
		log("Running protein differential analysis",1)
		self.series_df = pd.DataFrame()
		series_df = self.series_df
		df = pd.read_csv(mqrun_path+"proteinGroups.txt", sep = "\t")
		df = apply_filter('contaminants',None,['Majority protein IDs'],df) 
		df = apply_filter('reverse',None,['Majority protein IDs'],df) 
		if not silac:
			df = self.combine_replicates(df,biological_replicates,"")
			df = self.combine_replicates(df,conditions,'_mean')
			df = self.combine_replicates(df,conditions,'_err')

		df.to_csv('out_test.csv',sep="\t",index=False)
		series_df,avcols = generate_protein_series(df)
		series_df = dfappend(df,series_df,['Gene names'])
		series_df = self.apply_global_filters(series_df,avcols)
		series_df.to_csv(project_outpath+"/protein_diffexp.txt",sep="\t",index=False)
		subreport = self.generate_report(current_analysis,"Protein Level Analysis Report",df,series_df,avcols)
		log("____Protein differential analysis complete____\n",1)



	





if __name__ == "__main__":
	analysis = DifferentialAnalysis(run_name)
	#PhosphoAnalysis(analysis)
	ProteinAnalysis(analysis)





			



