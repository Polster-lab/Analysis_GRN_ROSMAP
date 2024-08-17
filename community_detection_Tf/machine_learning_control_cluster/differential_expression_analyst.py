import os
import pickle as pkl

import pandas as pd

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats



def differential_expression_analysis(output_dir, count_matrix, meta_data):

    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=count_matrix,
        metadata=meta_data,
        design_factors="condition",
        refit_cooks=True,
        inference=inference,
    )

    dds.deseq2()

    print(dds)
    print(dds.varm["dispersions"])
    print(dds.varm["LFC"])

    stat_res = DeseqStats(dds, inference=inference)
    print(stat_res.summary())
    stat_res.results_df.to_csv(os.path.join(output_dir, "results.csv"))



count_matrix = pd.read_csv('count_matrix_diff_analysis.csv',index_col=0)
meta_data = pd.read_csv('meta_data_diff_analysis.csv',index_col=0)
meta_data = meta_data[meta_data.condition.isin(['C','NCI'])]
#count_matrix = count_matrix.T
count_matrix = count_matrix[count_matrix.index.isin(meta_data.index)]
print(count_matrix)
print(meta_data)


differential_expression_analysis(output_dir='./differential_expression_output/C_NCI/', count_matrix=count_matrix, meta_data=meta_data)