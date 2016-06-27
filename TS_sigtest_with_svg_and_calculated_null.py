import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy.stats import rankdata
import copy
from compiler.ast import flatten
from scipy.stats import chisquare
from scipy.stats import mode
import random
from scipy.signal import savgol_filter

'''
The purpose of this script is to use the method of calculating significance for time series analysis
By W. Allen Wallis and Geoffry H. Moore., with two modifications.

It ignores rank, and attempts to converge on a local minimum p-value within 1 std err.

'''

def rank_expressions(exp_vals):
	ranks = rankdata(exp_vals)
	return ranks
	





def rank_sign_diff(arr):
	reparr = []
	for idx in range(1,len(arr)):
		current=arr[idx-1];next=arr[idx]
		if current>next:
			reparr.append(0)
			continue
		if current<next:
			reparr.append(1)
			continue
		if current==next:
			reparr.append(2)
	return reparr



def calc_phases(ranksign):
	arrs = [ranksign]
	EqualConsecRanks = True
	while EqualConsecRanks:
		EqualConsecRanks = False

		for arr_idx in range(len(arrs)):
			for val_idx in range(len(arrs[arr_idx])):
				if arrs[arr_idx][val_idx]==2:
					EqualConsecRanks=True
			
					arrs[arr_idx][val_idx] = 1; arrs.append(copy.deepcopy(arrs[arr_idx]))
					arrs[arr_idx][val_idx] = 0; arrs.append(copy.deepcopy(arrs[arr_idx]))
					arrs.remove(arrs[0])
					break
	return arrs
				



def calc_phase_lengths(arrs):
	phase_arrs = []
	for subarr in arrs:
		phase_lengths = []
		startval = subarr[0]
		fcount=1
		for idx in range(1,len(subarr)):
			if subarr[idx]==startval:
				fcount+=1
			else:
				phase_lengths.append(fcount)
				startval = subarr[idx]
				fcount=1
		phase_lengths.append(fcount)
		if np.array(arrs).shape[0]==1:
			return phase_lengths		
		phase_arrs.append(phase_lengths)
	return phase_arrs



def normalize(v):
	nfact = 1/np.sum(v)
	return np.dot(v,nfact)




def calculate_frequencies(arrs,bins):
	bins2cnt = dict(zip(bins,np.zeros(len(bins))))
	if type(arrs[0]) is int:
		for val in arrs:
			bins2cnt[val] = bins2cnt[val]+1
	else:
		for arr in arrs:
			for val in arr:
				bins2cnt[val] = bins2cnt[val]+1
		for key in bins2cnt:
			bins2cnt[key] = bins2cnt[key]/len(arrs)
	return normalize(bins2cnt.values())







def returnNull(time_series):
	rank_arr = np.array(np.zeros((time_series.shape[0],time_series.shape[1]))) #copy for array dimension instantiation
	cntr=0
	for col in time_series.T:
		rank_arr[:,cntr]=rank_expressions(col); cntr+=1
	signarr = np.array(np.zeros((rank_arr.shape[0],rank_arr.shape[1]-1)))
	cntr=0
	for row in time_series:
		signarr[cntr,:] = rank_sign_diff(row); cntr+=1

	phase_arrs=[]
	for row in signarr:
		arr = calc_phases(row)
		phase_arrs.append(arr)

	phase_counts = []
	for arrs in phase_arrs:
		phase_lengths = calc_phase_lengths(arrs)
		phase_counts.append(phase_lengths)
	bins =  set(flatten(np.arange(1,time_series.shape[1],1)))
	frequencies = [];
	for arrs in phase_counts:
		freqs = calculate_frequencies(arrs,bins)
		frequencies.append(freqs)
	popdist = calculate_frequencies(flatten(phase_counts),bins)
	popdist[popdist==0] = 1e-100
	return popdist,bins



def calc_freq(arr,bins):
	signarr = rank_sign_diff(arr);
	phase_lengths = calc_phase_lengths([signarr])
	freqs = calculate_frequencies(phase_lengths,bins)
	freqs[freqs==0] = 1e-100
	return freqs,phase_lengths,signarr





def return_null(n_bins):
	arrs = np.zeros((20000,n_bins))
	for row_ix in range(arrs.shape[0]):
		arrs[row_ix] = np.random.random_integers(0,1,n_bins)
	return returnNull(arrs)









def RUN(series_df,timecols,err_ub,err_lb,EF,n_iter,gene_names):


	times = [-1,0,0.5,1,1.5,2,2.5,5,10,15,20,30,60]

	maxCol=lambda x: max(x.min(), x.max(), key=abs)  #function to return row maxval
	maxvals = series_df[timecols].apply(maxCol,axis=1).apply(np.abs).values
	nulldist,bins = return_null(len(timecols))

	t_series = series_df[timecols].get_values()
	err_lb_series = series_df[err_lb].get_values()
	err_ub_series = series_df[err_ub].get_values()


	current_lb = err_lb_series[0,:]
	current_ub = err_ub_series[0,:]
	current_t_series = t_series[0,:]

	fits = []
	pvals = []
	stats = []
	most_sig_fits = np.array(np.zeros(series_df[timecols].shape),dtype=float)
	print most_sig_fits.shape
	for e in range(t_series.shape[0]):
		all_curves = []	

		all_pvals = []
		freq,phase_lengths,signarr = calc_freq(current_t_series,bins)
		init = chisquare(freq, nulldist)[1]
		chisq_stat=0
		mincnt=0
		current_lb = err_lb_series[e,:]
		current_ub = err_ub_series[e,:]
		current_t_series = t_series[e,:]
		svg_filtered = savgol_filter(current_t_series,13,4)

		most_sig_freq = freq
		most_sig_phase = phase_lengths
		most_sig_signs = signarr

		freq,phase_lengths,signarr = calc_freq(svg_filtered,bins)	
		stat, pval = chisquare(freq, nulldist)
		

		pvals.append(pval)
		stats.append(chisq_stat)
		most_sig_fits[e,:] = np.array(svg_filtered,dtype=float)
		'''
		uppers = np.abs(np.subtract(current_ub,current_t_series))
		lowers = np.abs(np.subtract(current_t_series,current_lb))
		plt.tick_params(axis='both', which='major', labelsize=20)
		plt.tick_params(axis='both', which='minor', labelsize=20)
		plt.xlabel("Time (Minutes)",fontsize=25)
		plt.ylabel("Expression log2(ratio)",fontsize=25)
		plt.plot(times,svg_filtered,color='blue',linewidth=4,alpha=0.5, label="SVG filtered")
		plt.errorbar(times,current_t_series,yerr=[lowers,uppers],color='red',linewidth=3, alpha=0.5,label="Observed Data")
		plt.legend()
		plt.show()
		'''

	fitcols=[]
	for s in range(most_sig_fits.shape[1]):
		series_df.insert(series_df.shape[1],"fit_"+str(s),most_sig_fits[:,s])
		fitcols.append("fit_"+str(s))


	series_df.insert(series_df.shape[1],"maxvals",maxvals)
	series_df.insert(series_df.shape[1],"pvals",pvals)
	series_df.insert(series_df.shape[1],"stats",stats)
	return series_df,fitcols
	





standard_cols=['leading_accession','gene_name','id']

timecols=['neg1_min^log2_average','0_min^log2_average','0.5_min^log2_average','1.0_min^log2_average','1.5_min^log2_average','2.0_min^log2_average','2.5_min^log2_average','5_min^log2_average','10_min^log2_average','15_min^log2_average','20_min^log2_average','30_min^log2_average','60_min^log2_average']


df = pd.read_csv('/home/jjacob/python/mq_analysis/spring_sixmoth_2016/analysis_20160424/report/phosphosite_diffexp.txt',sep="\t",low_memory=False)
errcols_lb = [(t.split("^")[0]+'^log2_stdev_lb') for t in timecols]
errcols_ub = [(t.split("^")[0]+'^log2_stdev_ub') for t in timecols]
outpath = "sampling_out.txt"
err_fraction=0.5 ## err_fraction>0
n_iter=50


series_df,fitcols = RUN(df,timecols,errcols_ub,errcols_lb,err_fraction,n_iter,df['gene_name'].values)  ##########start the run
series_df.to_csv(outpath,index=False,sep="\t")




