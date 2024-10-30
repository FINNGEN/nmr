import os,logging,sys,errno,gzip,mmap,math
import numpy as np
from itertools import combinations
from color_dict import color_dict
from collections import defaultdict as dd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
sns.set(palette='Set2')
import matplotlib.ticker as ticker
import scipy.cluster.hierarchy as sch

def file_exists(fname):
    '''
    Function to pass to type in argparse
    '''
    if os.path.isfile(fname):
        return str(fname)
    else:
        print(fname + ' does not exist')
        sys.exit(1)
 
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def plot_dist(df,fig_dir,prefix):
    """
    Basic optional plotting
    """
    for column in df.columns:
         print(column)
         fig, axs = plt.subplots(3, 1)
         sns.histplot(df[column],kde=True,ax=axs[0])
         renorm = np.log1p(df[column])
         sns.histplot(renorm,kde=True,ax=axs[1])
         normalized = (renorm - np.mean(renorm)) / np.std(renorm)
         sns.histplot(normalized,kde=True,ax=axs[2])
         fig.savefig(os.path.join(fig_dir,f'{prefix}_{column}.png'))
         plt.close()


def get_corr_cluster(df):
    print('correlation')
    c = df.corr()
    c[c<.9]=0
    print('done')
    
    
    print('clustering')
    X = df.corr().values
    d = sch.distance.pdist(X)   # vector of ('55' choose 2) pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    print('done.')
    columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
    d = {elem:int(ind[i]) for i,elem in enumerate(df.columns.tolist())}
    col_cluster = sorted(d.items(), key=lambda kv: kv[1])
    df_clust = df.reindex(columns, axis=1)

    return c,df_clust,col_cluster
    
         
def plot_corr(corr,df_cluster,fig_dir,prefix):
    plt.figure()
    sns.heatmap(corr, annot = True,vmin=0,vmax=1)
    fig_path = os.path.join(fig_dir,f'{prefix}_correlation.png')
    plt.savefig(fig_path)
    print('done.')
    plt.close()
    
    #CLUSTERED
    plt.figure()
    print('plotting')
    sns.heatmap(df_cluster.corr(), annot = True,vmin=0,vmax=1)
    plt.savefig(os.path.join(fig_dir,f'{prefix}_correlation_cluster.png'))
    print('done.')
    plt.close()

    
def plot_2d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 6,label_fontsize = 5,tag_column = "TAG",max_size = 3000,max_map = None,axis_legend =2,legend_location = "lower left",rescale = 4):
    '''
     Inputs:
    -- pc_file : name of file where to fetch data
    -- out_file : name of figure for saving
    -- tag_column : columns of the dataframe to plot
    -- tags: list of tags to plot in tag_column
    -- pc_columns : name of the columns in the file
    -- pc_tags : name of the pcs to plt
    -- colors_map : annotation to color_dict
    '''

    if not pc_tags: pc_tags = pc_columns

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(3,1)

    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0],sharex = ax1)
    ax3 = fig.add_subplot(gs[2,0],sharex = ax1)
    axes=[ax1,ax2,ax3]

    
    # init tag sizes
    sizes = dd(lambda :.5)
    if size_map:
        for key in size_map:sizes[key] = size_map[key]
    alphas = dd(lambda:0.2)
    if alpha_map:
        for key in alpha_map:alphas[key] = alpha_map[key]
    max_sizes = dd(lambda:max_size)
    if max_map:
        for key in max_map: max_sizes[key] = max_map[key]

    if not color_map:
        color_maps = list(color_dict[len(tags)]['qualitative'].keys())
        cm = 'Set1' if 'Set1' in color_maps else color_maps[0]
        color_map = {tag:color_dict[len(tags)]['qualitative'][cm][i] for i,tag in enumerate(tags)}

    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        max_data = max_sizes[tag]
        n = min(max_data,len(tag_data))
        tag_data = tag_data.sample(n=n)
        color = color_map[tag]
        size = sizes[tag]
        alpha = alphas[tag]
        print(tag,len(tag_data),alpha,size,n)
        for i,pcs in enumerate(list(combinations(pc_columns,2)) ):
            ax = axes[i]
            plot_x,plot_y = tag_data[pcs[0]],tag_data[pcs[1]]
            ax.scatter(plot_x,plot_y, s= size,alpha = alpha,color = color,label = tag)        


    for i,tags in enumerate(list(combinations(pc_tags,2))):
        tag1,tag2 = tags
        ax = axes[i]
        ax.set_xlabel(tag1)
        ax.set_ylabel(tag2)
        trim_axis(ax)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
        ax.tick_params(axis='both', which='major', labelsize=3)
        ax.tick_params(axis='both', which='minor', labelsize=3)

    leg_ax = axes[axis_legend]
    leg = leg_ax.legend(loc=legend_location, numpoints=1, fancybox = True,prop={'size': legend_fontsize})
    for lh in leg.legend_handles:
        lh.set_alpha(1)
        if rescale:
            lh._sizes = [lh._sizes[0]*7]
        else:
            lh._sizes = [50]

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.savefig(out_file)
    fig.savefig(out_file.replace('.pdf','.png'),dpi=300)
    plt.close()


    

def trim_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    try:
        ax.get_yaxis().tick_left()
    except:
        pass
