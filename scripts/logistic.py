import numpy as np
import pandas as pd
import argparse,os
from utils import make_sure_path_exists,file_exists
from sklearn.linear_model import LogisticRegression


def read_df(args):
    " DO ONE BIG ROUND AND THEN FORGET ABOUT IT,IDEALLY"

    out_file = os.path.join(args.out,'pheno_cov_data.txt')
    if  os.path.isfile(out_file):
        merged_df = pd.read_csv(out_file,index_col='FINNGENID',sep='\t')

    else:
        # READ IN NMR DATA AND DROP DUPLICATES
        print("READING NMR")
        data_df = pd.read_csv(args.data_file,sep='\t',index_col=['FINNGENID','SAMPLE_COLLECTION'])
        var_cols = data_df.columns.tolist()
        print(data_df)
        print("done")
        
        print("READING PHENO")
        # EXCLUDE OVERLAPPING COLS
        cols = [col for col in list(pd.read_csv(args.pheno_file, nrows=1,sep='\t')) if col not in var_cols]
        # READ IN ONLY FINNENGID
        tmp_df = pd.read_csv(args.pheno_file,sep='\t',usecols =['FINNGENID'])
        # GET ROWS WITH SHARED INDEX
        shared_df = tmp_df.reset_index().merge(data_df,on='FINNGENID').set_index('index')
        skiplist = [elem +1  for elem in set(range(0, len(tmp_df))) - set(shared_df.index)]
        pheno_df = pd.read_csv(args.pheno_file,sep='\t',index_col='FINNGENID',usecols =['FINNGENID']+ cols,skiprows=skiplist)
        print(pheno_df)
        
        print("MERGE")
        merged_df = data_df.merge(pheno_df,how='left',on='FINNGENID')
        merged_df = merged_df[~merged_df.index.duplicated(keep='first')]
        print(merged_df)
        merged_df.to_csv(out_file,sep='\t')
        
    return merged_df




def main(args):

    df = read_df(args)
                 



if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="Basic analysis of nmr data.")
    parser.add_argument("--data_file", type=file_exists, help="List of files to merge.", default='/mnt/disks/data/nmr/data/data.tsv')
    parser.add_argument("--pheno-file", type=file_exists, help="Path to input raw file. File should be tsv.", default='/home/pete/r12/pheno/R12_COV_PHENO_V2.txt.gz')
    parser.add_argument("--test", action='store_true', help="Uses limited files and rows")
    parser.add_argument("--phenos", type=file_exists, help="List of cols to analyze.",default= '/home/pete/r12/pheno/R12_analysis_endpoints.txt') 
    parser.add_argument('-o', "--out", type=str, help="Folder in which to save the results (default = current working directory)", default='/mnt/disks/data/nmr/regression/')
    args = parser.parse_args()

    #CREATE NECESSARY PATHS AND DEFINE OUT FILES
    make_sure_path_exists(args.out)

    # READ IN 
    with open(args.phenos) as i:args.pheno_list=  [elem.strip() for elem in i.readlines()]
    main(args)
