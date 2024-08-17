import argparse
import pandas as pd
import os
from tqdm import tqdm
from netZooPy import condor
import numpy as np

def main(args):

    
    print(f"start:, {args.start}")
    print(f"end:, {args.end}")

    AD_filtered_data = pd.read_csv('AD_filtered_data.csv', index_col=0)
    genes= []
    TF = []
    for i in AD_filtered_data.index:
        TF.append(i.split('_')[0])
        genes.append(i.split('_')[1])
    AD_patients_list = AD_filtered_data.columns.tolist()
    AD_sample_1 = AD_filtered_data[AD_patients_list[0]].reset_index(drop= True)
    AD_sample_1 = pd.DataFrame(AD_sample_1)
    AD_sample_1[0] = TF
    AD_sample_1[1] = genes
    AD_sample_1 = AD_sample_1[[0,1,AD_patients_list[0]]]

    
    condor_object = condor.condor_object(dataframe = AD_sample_1)
    condor_object.initial_community()
    condor_object.brim()

    tf_index = condor_object.reg_memb['reg'].tolist()
    gene_index = condor_object.tar_memb['tar'].tolist() 

    rows = len(tf_index) # number of rowsprint
    cols = len(tf_index) # number of columns

    # Create a DataFrame with all zero values
    

    #gene_matrix = pd.DataFrame(np.zeros((rows, cols)), columns=[f'{i}' for i in gene_index])
    #gene_matrix.index = gene_index
    



    subset_AD_patients= AD_patients_list[args.start:args.end]


    for i in tqdm(subset_AD_patients):
        tf_matrix = pd.DataFrame(np.zeros((rows, cols)), columns=[f'{i}' for i in tf_index])
        tf_matrix.index = tf_index
        AD_sample_1 = AD_filtered_data[[i]].reset_index(drop=True)
        AD_sample_1[0] = TF 
        AD_sample_1[1] = genes
        AD_sample_1 = AD_sample_1[[0,1,i]]
        remove_genes = list(set(TF).intersection(set(genes)))
        AD_sample_1 = AD_sample_1[~AD_sample_1[1].isin(remove_genes)]
        for k in tqdm(range(10)):
            condor_object.initial_community()
            condor_object.brim()
            
        
            for index, row in tf_matrix.iterrows():
                    for col in tf_matrix.columns:
                        #print(f'Row: {index}, Column: {col}, Value: {row[col]}')
                        #print(index)
                        #print(col)
                        if condor_object.reg_memb[condor_object.reg_memb.reg == index].community.values == condor_object.reg_memb[condor_object.reg_memb.reg == col].community.values:
                            tf_matrix.loc[index,col] = tf_matrix.loc[index,col] + 1
        
        tf_matrix.to_csv('./tf_matrix_similarity/'+str(i)+'.csv', index= True)





    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-s', '--start', type=int, default=0)
    parser.add_argument('-e', '--end',type=int, default=91)
    
    args = parser.parse_args()
    main(args)
