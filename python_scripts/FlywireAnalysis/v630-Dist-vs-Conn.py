import pandas as pd
import numpy as np
import gc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix

sns.set_context("talk")

def gen_proximity_model(dev_lim):
    
    print('loading raw data..')
    
    fdir = "data/synapse_coordinates.csv"
    # load id as integers and allowing NaN
    syn_xyz = pd.read_csv(fdir, dtype={'pre_root_id': pd.Int64Dtype(), 'post_root_id': pd.Int64Dtype(), 
                                    'x': np.int32, 'y': np.int32, 'z': np.int32})

    fdir = "data/connections.csv"
    syn_table = pd.read_csv(fdir)

    syn_xyz.fillna(method='ffill', inplace=True)

    syn_xyz.pre_root_id = syn_xyz.pre_root_id.astype(np.int64)
    syn_xyz.post_root_id = syn_xyz.post_root_id.astype(np.int64)


    print('calculating centers...')
    pre_centers = syn_xyz.groupby('pre_root_id')[['x', 'y', 'z']].mean()
    post_centers = syn_xyz.groupby('post_root_id')[['x', 'y', 'z']].mean()

    del syn_xyz

    print('calculating pairwise distance...')
    pdist = pd.DataFrame(
        distance_matrix(pre_centers.loc[pre_centers.index[:dev_lim]], 
                                        post_centers.loc[post_centers.index[:dev_lim]],
                        threshold= 10000 * 10000 * 3
                        ),
        columns=post_centers.index[:dev_lim],
        index=pre_centers.index[:dev_lim]
        )

    print('aggragating connections...')
    dev_mask = syn_table.pre_root_id.isin(pre_centers.index[:dev_lim]) & syn_table.post_root_id.isin(post_centers.index[:dev_lim])
    pconn = syn_table[dev_mask].groupby(['pre_root_id', 'post_root_id']).sum()[['syn_count']]
    pconn['connected'] = 1
    del syn_table
    
    print('reformating distance matrix...')
    pdist.reset_index(inplace=True)
    pdist.reset_index(drop = True, inplace=True)

    gc.collect()
    pdist = pdist.melt(id_vars=['pre_root_id'], var_name='post_root_id', value_name='dist')

    # remove autapses since they are not in our dataset
    print('removing autapses...')
    mask_autapses = (pdist.pre_root_id == pdist.post_root_id)
    pdist = pdist[~mask_autapses]

    print('merging...')
    dist_vs_syn = pdist.merge(pconn, how='left', on=['pre_root_id', 'post_root_id']).fillna(0)

    n_bins = 50
    print(f'qcut data in {n_bins} bins...')
    dist_vs_syn['dist_bins'] = pd.qcut(dist_vs_syn.dist, q=np.linspace(0,1,n_bins))

    print(f'calculating for mean and error for each bin...')
    error = lambda x: np.std(x) / np.sqrt(len(dist_vs_syn) // n_bins)
    dist_vs_syn_fn =  dist_vs_syn.groupby('dist_bins').agg({'syn_count': [('mean', np.mean), ('error', error)], 
                                            'connected': [('mean', np.mean), ('error', error)]})
    # dist_vs_syn_fn                                        
    if dev_lim is None:
        print(f'saving binned stats...')
        dist_vs_syn_fn.to_pickle(f'saved/dist_vs_syn_fn.pkl')
        print(f'saving raw stats...')
        dist_vs_syn.to_pickle(f'saved/dist_vs_syn_raw.pkl')
    else:
        print(f'saving binned stats...')
        dist_vs_syn_fn.to_pickle(f'saved/dist_vs_syn_fn_n={dev_lim}.pkl')
        print(f'saving raw stats...')
        dist_vs_syn.to_pickle(f'saved/dist_vs_syn_raw_n={dev_lim}.pkl')






if __name__ == '__main__':
    dev_lim = 10000
    gen_proximity_model(dev_lim)    