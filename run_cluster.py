import initialization_paths
from functions import cluster_funcs

print('--- --- --- send on the cluster --- --- --- ')

#cluster_funcs.create_qsub('convert_data_to_BIDS', 'BIDS', 'BIDS', queue='NSpin_long')
cluster_funcs.create_qsub('test_SVM_decoder_on_each_sequence', 'SVM_proj_norm', 'SVM_proj_norm', queue='Nspin_bigM')

print('--- --- --- finished sending on the cluster --- --- --- ')

