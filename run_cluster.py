import initialization_paths
from functions import cluster_funcs

print('--- --- --- send on the cluster --- --- --- ')

cluster_funcs.create_qsub('convert_data_to_BIDS', 'BIDS', 'BIDS', queue='Nspin_run32')

print('--- --- --- finished sending on the cluster --- --- --- ')
