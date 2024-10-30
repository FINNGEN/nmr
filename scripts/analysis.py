import numpy as np
import pandas as pd
import os,argparse,random
from utils import file_exists,make_sure_path_exists,plot_dist,plot_corr,plot_2d,get_corr_cluster
from sklearn.decomposition import PCA

TEST_COLS=['FINNGENID','SAMPLE_COLLECTION','ACETATE','ALA','ALBUMIN','SPECTROMETER']

def merge_data(file_list,pheno_file,test=True):
    print("MERGE DATA")
    df = pd.DataFrame()
    for f in file_list:
        nrows = 100 if args.test else None
        tmp= pd.read_csv(f,sep='\t',nrows = nrows)
        df = pd.concat((df,tmp),axis=0)

    print(df)
    df = df.copy()
    # extract TAG and YEAR from sample collection field
    df[['TAG','YEAR']] = df['SAMPLE_COLLECTION'].str.split(" ",n=1,expand=True)
    print(df)
    print('ADD SEX COLUMN')
    min_pheno=pd.read_csv(pheno_file,sep='\t',usecols=['FINNGENID','SEX','APPROX_BIRTH_DATE'],index_col='FINNGENID')
    min_pheno['APPROX_BIRTH_DATE'] = pd.to_datetime(min_pheno['APPROX_BIRTH_DATE'])
    df = df.merge(min_pheno,how='left',on='FINNGENID')
    #CALCULATED APPROX AGE
    df['AGE'] = round((pd.to_datetime(df.YEAR,format='%Y') - pd.to_datetime(df.APPROX_BIRTH_DATE)).dt.days/365,2)
    df =df.drop(columns=['YEAR','APPROX_BIRTH_DATE'])
    df = df.set_index(['FINNGENID','SAMPLE_COLLECTION'])
    print(df)
    return df


def pca(data_df,out_dir,n_components = 10,renorm = True):
    """
    Performs PCA and adds PC columns to data
    """
    
    print(data_df)
    pca = PCA(n_components=n_components)
    pcs =  pca.fit_transform(data_df)
    # saves weights
    pc_cols = [f"PC{i+1}" for i in range(n_components)]
    weight_df = pd.DataFrame()
    weight_df[pc_cols] = pca.components_.T
    weight_df.index = data_df.columns.tolist()
    weight_df.to_csv(os.path.join(out_dir,'weights.tsv'),sep='\t')
    # ADD PCs to df
    pca_df = pd.DataFrame(pcs, columns=pc_cols,index=data_df.index)
    data_df = pd.concat((data_df,pca_df),axis=1)
    print(data_df)
    return data_df



def main(args):

    # MERGE DATA
    with open(args.file_list) as i:file_list=[elem.strip() for elem in i.readlines()]
    df = merge_data(file_list,args.pheno_file,args.test)
    print(f"data dumped to {args.out_file}")
    df.to_csv(args.out_file,sep='\t')

    # read in analysis cols (i.e. float cols)
    with open(args.analysis_cols) as i:analysis_cols=  [elem.strip() for elem in i.readlines()]
    final_cols = sorted(list(set(analysis_cols) & set(df.columns)),key=str.lower)
    data_df = df[final_cols]
    meta_df = df[[col for col in df.columns if col not in final_cols]]

    # CORRELATION BETWEEN MAKERS
    corr,df_cluster,col_cluster = get_corr_cluster(data_df)
    with open(os.path.join(args.out,'nmr_correlation_clusters.txt'),'wt') as o:
        for elem in col_cluster:
            o.write('\t'.join(map(str,elem)) + '\n')

    # PLOT
    if args.plot:
        print("PLOT")
        plot_dist(data_df,args.fig_path,'nmr')
        plot_corr(corr,df_cluster,args.fig_path,'nmr')

    renorm = np.log1p(data_df)
    data_df = (renorm - renorm.mean()) / renorm.std()
    data_df = data_df.fillna(data_df.mean())
        
    print("PCA")
    data_df = pca(data_df,args.out,args.n_pca)
    plot_df = data_df.copy().assign(TAG=meta_df.TAG)
    plot_2d(plot_df,os.path.join(args.fig_path,'PCA_plot.pdf'),set(plot_df.TAG),label_fontsize=3)


    age_sex_df = df[['AGE','SEX']]
    final_df = data_df.join(age_sex_df)
    print(final_df)
    final_df.to_csv(args.out_file.replace('.tsv','_renormed.tsv'),sep='\t')

    
    
if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="Basic analysis of nmr data.")
    parser.add_argument("--file-list", type=file_exists, help="List of files to merge.", default='/mnt/disks/data/nmr/file_list.txt')
    parser.add_argument("--pheno-file", type=file_exists, help="Path to input raw file. File should be tsv.", default='/home/pete/r12/phenotype_internal_1.0/data/finngen_R12_minimum_internal_1.0.txt.gz')
    parser.add_argument("--test", action='store_true', help="Uses limited files and rows")
    parser.add_argument("--n_pca", type= int,default = 10, help='number of PCs')
    parser.add_argument("--plot", action='store_true', help="Plot or not distributions")
    parser.add_argument("--analysis-cols", type=file_exists, help="List of cols to analyze.", default='/mnt/disks/data/nmr/analysis_cols.txt')
    parser.add_argument('-o', "--out", type=str, help="Folder in which to save the results (default = current working directory)", default="/mnt/disks/data/nmr/data")
    args = parser.parse_args()

    #CREATE NECESSARY PATHS AND DEFINE OUT FILES
    print(args.out)
    args.fig_path = os.path.join(args.out,'figs/')
    make_sure_path_exists(args.fig_path)
    args.out_file =os.path.join(args.out,'nmr_data.tsv')
    main(args)
