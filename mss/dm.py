import pandas as pd
import numpy as np
from itertools import groupby
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
from sklearn import cluster,mixture
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
import scipy

def data_prep(d_input, blank_keyword, simp_summary = False,svb_thres=10, empty_thres=0,rt_range=[0, 30], mz_range=[0, 1200], sn_thres=3, score_thres=0, area_thres=5000):
    '''
    The function is used to clean the dataframe according to user setting
    blank_keyword: part of string from column that indicates the column is a blank sample
    svb_thres: sample vs blank thres
    empty_thres: empty cell thres in a row
    cv_thres: as all sample is in triplicate, calculate the CV for every triplicate sample set #Needs to be updated in case there is no triplicate samples
    rt_range: rt filter
    mz_range: mz filter
    sn_thres: signal/noise column thres
    score_thres: score column thres
    area_thres: count for max peak area from each row
    '''
    d_thres = d_input[d_input[d_input.columns[4:]].max(1) >= area_thres]
    
    d_thres = d_thres[(d_thres['Average RT (min)'] > rt_range[0]) & (d_thres['Average RT (min)'] < rt_range[1])]
    d_thres = d_thres[(d_thres['Average m/z'] > mz_range[0]) & (d_thres['Average m/z'] < mz_range[1])]
    d_thres = d_thres[d_thres['Average sn'] >= sn_thres]
    d_thres = d_thres[d_thres['Average score'] >= score_thres]
    d_thres.reset_index(inplace=True, drop=True)
    
    col_blank = []
    for key in blank_keyword:
        # Get column name if it contains blank indicating strings
        col_blank.extend([col for col in d_thres.columns if key in col])
        
    col_sample = [col for col in d_thres.columns if col not in col_blank]
    # Sample maximum area vs Blank average area to count for svb
    d_sample = d_thres[d_thres[col_sample[4:]].max(axis=1) / d_thres[col_blank].mean(axis=1) > svb_thres][col_sample] 
    d_sample.reset_index(inplace=True)
    d_sample.drop(columns=['index'],inplace=True)
    
    # Get a list of triplicate, every triplicate is in a sublist
    #Sample: [[a1,a2,a3],[b1,b2,b3]]
    #Note: the triplicate parsing is now only used '_' which needs update in the future
    #d_transpose['dilu_vol'] = d_transpose['dilu_vol'].apply(lambda x : x.replace('-','_')) in case people use '-' for parsing
    trip_list = [list(i) for j, i in groupby(d_sample.columns[4:], lambda a: a.split('_')[:-1])] 
    trip_list = [i for i in trip_list if len(i)>=2] #filter out columns that is not in triplicate -- sample naming issue

    for triplicate in tqdm(trip_list):
        for row in d_sample[triplicate].itertuples(): # Loop for every sets of triplicates
            if row[1:].count(0) > empty_thres:
                d_sample.loc[row.Index, triplicate] = 0 # if more than thres, then set all three values to 0
#             elif np.mean(row[1:]) != 0:
#                 if np.std(row[1:]) / np.mean(row[1:]) > cv_thres:
#                     d_sample.loc[row.Index, triplicate] = 0 #need verify, not work for now
            else:
                pass
            
    d_sample = d_sample[~(d_sample[d_sample.columns[4:]]==0).all(1)] #clean rows with all 0
    if simp_summary == True:
        simp_dict={}
        for i, column in enumerate(trip_list):
            avg = d_sample[column].mean(1)
            cv = d_sample[column].std(1) / d_sample[column].mean(1) #optional display CV
            simp_dict.update({column[0][:-2]:avg, ' CV #' + str(i):cv})
        d_result = pd.DataFrame(simp_dict)
        d_result = pd.concat([d_sample[d_sample.columns[:4]], d_result], axis=1)
    elif simp_summary == False:
        d_result = d_sample.copy()
    
    return d_result



def ms_cluster(d_input, select_keyword, normalization='linear', visual=False, d_reduce=True, d_reduce_method='tsne', perplexity=20, cluster_method='dbscan',eps=0.8,min_samples=10):
    '''
    Function for direct clustering:
    normalization method: linear, zscore, log
    d_reduce: if perform the dimension reduction algorithm, method: only tsne is avilable now
    perplexity: parameter for tsne
    cluster_method: dbscan, later will update optic and spectrum
    eps: parameter for dbscan, threshold of radius that used to count neighbours
    min_samples: general parameter for clustering, min neighbourhoods to be counted as a cluster
    '''
    col_select = []
    for key in select_keyword:
        col_select.extend([col for col in d_input.columns if key in col])
    d_clu = d_input[col_select]
    
    c_data = d_clu.values
    c_norm = []
    #Performs normalization
    np.seterr(divide='ignore', invalid='ignore') #silent the warning -- but divide by 0 still exist
    for row in c_data:
        if normalization == 'linear':
            c_norm.append(row/max(row))
        elif normalization == 'zscore':
            c_norm.append((row-np.mean(row))/np.std(row))
        elif normalization == 'log':
            row[row==0]=1
            c_norm.append(np.log10(row)/np.log10(max(row)))
        else:
            pass
    #Clean up dataframe
    c_norm = np.asarray(c_norm)
    d_norm = pd.DataFrame(c_norm)
    d_norm['index']=d_input.index
    d_norm.set_index('index',inplace=True)
    d_norm.dropna(how='all',inplace=True)
    
    if d_reduce == True:
        if d_reduce_method == 'tsne':
            model = TSNE(learning_rate=100,perplexity=50,n_iter=1000) #Tune perplexity and n_iter
            transformed = model.fit_transform(d_norm)
            d_feature = transformed.copy()
        else:
            pass
    elif d_reduce == False:
        d_feature = d_norm.copy()
    else:
        pass
    
    if cluster_method == 'dbscan':
        dbscan = cluster.DBSCAN(eps=eps, min_samples=min_samples).fit(d_feature)
        labels = dbscan.labels_
        unique_labels = set(dbscan.labels_)
        
        if visual == True:
            for i,k in enumerate(unique_labels):
                indexlist = list(np.argwhere(labels==k).reshape(1,-1)[0])
                sns.clustermap(d_norm.iloc[indexlist].values,cmap='Reds',col_cluster=True,yticklabels=False,xticklabels=False,figsize=(5,5))
                plt.title(str(dbscan)+'label='+ str(k))
                plt.show()
        else:
            pass
        d_init = d_input.copy()
        d_label = d_init.loc[d_norm.index] #Use the index to match back to the original datasheet
        d_label.insert(4,"label", dbscan.labels_.tolist())
    elif cluster_method == 'optics':
        optics = cluster.OPTICS(min_samples=min_samples).fit(d_feature)
        labels = optics.labels_
        unique_labels = set(optics.labels_)
        if visual == True:
            for i,k in enumerate(unique_labels):
                indexlist = list(np.argwhere(labels==k).reshape(1,-1)[0])
                sns.clustermap(d_norm.iloc[indexlist].values,cmap='Reds',col_cluster=True,yticklabels=False,xticklabels=False,figsize=(5,5))
                plt.title(str(optics)+'label='+ str(k))
                plt.show()
        else:
            pass
        d_init = d_input.copy()
        d_label = d_init.loc[d_norm.index] #Use the index to match back to the original datasheet
        d_label.insert(4,"label", optics.labels_.tolist())
    else:
        pass
    
    #Post filter -- filter out features that present in other sources but not SR520 -- keep it open for now
    #If activate add one more variable:source_keyword
#     col_source = []
#     for key in source_keyword:
#         col_app = [col for col in d_thres.columns if key in col]
#         col_source += col_app
#     col_rest = [col for col in d_label.columns if col not in source][5:]
#     d_label[col_app].max(1) / d_label[col_rest].max(1)
    
    return d_label


def trend_calc(d_input, select_keyword, min_size=5, normalization='linear', method='pearsonr',visual=True):

    """This function calculates clustering based on the pearson correlation.
    It takes in a dataframe and a user defined value for what qualifies as a cluster.
    User can choose whether or not to have a visual plot of the scatter with True/False."""
    col_select = []
    for key in select_keyword:
        col_app = [col for col in d_input.columns if key in col]
        col_select += col_app
    d_clu = d_input[col_select]
    
    c_data = d_clu.values
    c_norm = []
    for row in c_data:
        if normalization == 'linear':
            c_norm.append(row/max(row))
        elif normalization == 'zscore':
            c_norm.append((row-np.mean(row))/np.std(row))
        elif normalization == 'log':
            row[row==0]=1
            c_norm.append(np.log10(row)/np.log10(max(row)))
    c_norm = np.asarray(c_norm)
    d_norm = pd.DataFrame(c_norm)
    d_norm['index']=d_input.index
    d_norm.set_index('index',inplace=True)
    d_norm.dropna(how='all',inplace=True)
    
    #Post treatment to fit the d_norm into original codes
    d_norm.insert(0,"RT", d_input['Average RT (min)'].tolist())
    d_norm.insert(1,"MZ", d_input['Average m/z'].tolist())
    d_norm = d_norm.reset_index(drop=True)
    
    
    #Original codes
    cluster = [] # individual cluster holder
    cluster_sum = [] # total clusters
    drop_list = [] # rows that are dropped from the df
    noise = [] # list for containing noise features
    while len(d_norm) > 0:
        for row in range(len(d_norm)):
            feature_1 = d_norm.iloc[0]
            feature_2 = d_norm.iloc[row]
            if method == 'pearsonr':
                corr, p_val = scipy.stats.pearsonr(d_norm.iloc[0, 2:], d_norm.iloc[row, 2:]) 
            elif method == 'mannwhitneyu':
                corr, p_val = scipy.stats.mannwhitneyu(d_norm.iloc[0, 2:], d_norm.iloc[row, 2:]) 
            elif method == 'kruskal':
                corr, p_val = scipy.stats.kruskal(d_norm.iloc[0, 2:], d_norm.iloc[row, 2:]) 
            if p_val < 0.05:
                drop_list.append(row)
                cluster += [feature_2]
            else:
                pass
        if len(cluster) <= min_size:
            noise += [cluster]
            cluster = []
        else:
            cluster_sum += [cluster]
            cluster = []
        d_norm = d_norm.drop(drop_list)
        d_norm = d_norm.reset_index(drop=True)
        drop_list = []
    append_list = []
    for i in range(len(cluster_sum)):
        for j in range(len(cluster_sum[i])):
            cluster_sum[i][j].loc['Score']= i
            listing = np.array(cluster_sum[i][j])
            append_list.append(listing)
    cluster_df = pd.DataFrame(append_list) #Add columns use d_clu
    append_list2 = []
    for k in range(len(noise)):
        for l in range(len(noise[k])):
            noise[k][l].loc['Score']= -1
            listing2 = np.array(noise[k][l])
            append_list2.append(listing2)
    noise_df = pd.DataFrame(append_list2)
    final_df = pd.concat([cluster_df, noise_df])
    final_df = final_df.reset_index(drop=True)
    if visual == True:
        labels = final_df.iloc[:,-1:].values.reshape(1,-1)[0]
        unique_labels = set(labels)
        for i,k in enumerate(unique_labels):
            indexlist = list(np.argwhere(labels==k).reshape(1,-1)[0])
            sns.clustermap(final_df.iloc[indexlist,2:-1].values,cmap='Reds',col_cluster=True,yticklabels=False,xticklabels=False,figsize=(5,5))
            plt.title('trend'+'label='+ str(k))
            plt.show()
    else:
        pass
    return final_df

def source_label(d_input, sourcelist,area_thres=5000, concat = True): #noise removal only based on sourcelist cols
    np.seterr(divide='ignore', invalid='ignore')
    #source labeling
    d_result = d_input.copy()
    source_col=[]
    for s in sourcelist:
        source = [col for col in d_input.columns if s in col]
        source_col.append(source)
    simp_dict={}
    for i, column in enumerate(source_col):
        avg = d_result[column].mean(1)
        cv = d_result[column].std(1) / d_result[column].mean(1) #optional display CV
        cv_nan=np.isnan(cv)
        cv[cv_nan]=0.0 #replace nan with 0
        simp_dict.update({sourcelist[i]:avg, str(sourcelist[i])+' Cv':cv})
    d_summary = pd.DataFrame(simp_dict)
    d_summary['source']="NA"
    for row in d_summary.itertuples():
        sourcelabel = list(d_summary.columns[[col_index for col_index, peak_avg in enumerate(row[1:-1]) if peak_avg >= area_thres]])
        if len(sourcelabel) != 0:
            labelstr = ','.join(sourcelabel)
            d_summary.at[row.Index,'source'] = labelstr
    if concat == True:
        d_concat = pd.concat([d_result, d_summary], axis=1)
    elif concat == False:
        d_concat=d_result.copy()
        d_concat['source'] = d_summary['source']
    
    return d_concat

def source_report(d_input, source_key, mix_key, method='multiple', pa_thres=10000, CV_thres=2): #source key needs to be the same as source_label above
    #**only take the concat dataframe from labeling function
    #prefilter & dataframe arrangement
    d_mix = d_input[(d_input[[col for col in d_input.columns if 'Cv' in col]] <= CV_thres).all(1)]# all cv should below thres in order to be checked
    d_simp = d_mix[source_key]
    print('Threshold set to', pa_thres)
    mix_col = []
    for key in mix_key:
        mix_col.extend([col for col in d_mix.columns if 'Mix' in col])
    if len(mix_col) == 0:
        print("didn't find mixture by keyword!")

    d_st = pd.DataFrame(mix_col)

    c_name = ['sample']
    for source in source_key:
        result = []
        for col in mix_col:
            n_feature = sum(d_mix[d_mix[col] >= pa_thres]['source'].str.contains(source))
            cov_score = n_feature / sum(d_mix['source'].str.contains(source))
            if method == 'single':
                mix = d_mix[d_mix['source'] == source][col]
                s_simp = d_simp[d_simp['source'] == source][source]
            elif method == 'multiple':
                mix = d_mix[d_mix['source'].str.contains(source)][col]
                s_simp = d_simp.loc[d_mix[d_mix['source'].str.contains(source)].index][source]
            match_index = [i for i, j in enumerate(mix) if j >= pa_thres]
            dilu = mix.iloc[match_index] / s_simp.iloc[match_index]
            ratio_score = np.average(dilu[dilu<1])
            result.append([n_feature, cov_score, ratio_score])
        d_st = pd.concat([d_st, pd.DataFrame(result)], axis = 1)
        c_name.extend(['n_'+str(source), 'cover_s', 'ratio_s'])
    d_st.columns = c_name
    
    return d_st