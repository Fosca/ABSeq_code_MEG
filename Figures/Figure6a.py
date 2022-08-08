"""
==========================================================
  PLOT GFP FOR DIFFERENT CATEGORIES OF TRIALS
  + CORRELATION WITH COMPLEXITY
==========================================================

1 - Compute (or load) GFP for each subject and 3 categories of trials (items):
habituation, standard, violation. Epochs are mapped to magnetometers data

2 - Plot data:
    - line plot 7 sequences
    - correlation with complexity: heatmap
"""

# ---- import the packages -------
import sys
sys.path.append("/neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/ABSeq_scripts/")
import numpy as np
import pickle
from functions import article_plotting_funcs, utils, epoching_funcs
import config
import os.path as op

#_______________________________________________________________________________________________________________________
def compute_GFP_asmag_allparticipants(results_path = op.join(config.result_path, 'Corr_GFPxComplexity', 'items')):
    # Empty dictionaries to fill
    gfp_data = {}
    for ttype in ['habituation', 'standard', 'violation']:
        gfp_data[ttype] = {}
        gfp_data[ttype] = {}
        for seqID in range(1, 8):
            gfp_data[ttype][seqID] = []

    # Extract the data: subjects loop
    for subject in config.subjects_list:
        print('-- Subject -- \n' + subject)
        print(' -- LOAD THE EPOCHS and remap to only magnetometers -- #')
        epochs = epoching_funcs.load_epochs_items(subject, cleaned=True, AR_type='global')
        epochs = epochs.as_type('mag', mode="accurate")
        for ttype in ['habituation', 'standard', 'violation']:
            print('  ---- Trial type ' + str(ttype))
            for seqID in range(1, 8):
                print(' -- Computing GFP for sequence ' + str(seqID)+' --\n')
                if ttype == 'habituation':
                    epochs_subset = epochs['TrialNumber <= 10 and SequenceID == ' + str(seqID)].copy()
                elif ttype == 'standard':  # For standards, taking only items from trials with no violation
                    epochs_subset = epochs['TrialNumber > 10 and SequenceID == ' + str(seqID) + ' and ViolationInSequence == 0'].copy()
                elif ttype == 'violation':
                    epochs_subset = epochs['TrialNumber > 10 and SequenceID == ' + str(seqID) + ' and ViolationOrNot == 1'].copy()
                ev = epochs_subset.average()
                # Compute GFP
                gfp = np.sqrt(np.sum(ev.copy().pick_types(eeg=False, meg='mag').data ** 2, axis=0))
                # Store gfp each seq
                gfp_data[ttype][seqID].append(gfp)
    # Keep "times" in the dict
    gfp_data['times'] = ev.times
    utils.create_folder(results_path)
    # Save the data
    with open(op.join(results_path, 'gfp_each_seq_data.pickle'), 'wb') as f:
        pickle.dump(gfp_data, f, pickle.HIGHEST_PROTOCOL)

#_______________________________________________________________________________________________________________________
def plot_and_save_GFP_stats(results_path = op.join(config.result_path, 'Corr_GFPxComplexity', 'items')):

    with open(op.join(results_path, 'gfp_each_seq_data.pickle'), 'rb') as f:
        gfp_data = pickle.load(f, pickle.HIGHEST_PROTOCOL)

    f = open(op.join(config.fig_path, 'GFP','statistics.txt'), 'w')
    #  ============== HABITUATION PLOTS ============== #
    data_7seq = np.dstack(gfp_data['habituation']['mag'].values())
    data_7seq = np.transpose(data_7seq, (2, 0, 1))
    # Data line plot 7seq
    article_plotting_funcs.plot_7seq_timecourses(data_7seq, gfp_data['times'] * 1000, save_fig_path='GFP/', fig_name='GFPxComplexity_Habituation', suffix='',
                                                 pos_horizontal_bar=0.47, plot_pearson_corrComplexity=True, chance=None, xlims=[-50, 350], ymin=0, ylabel='GFP', logger=f)
    # Correlation with complexity heatmap
    pearsonr = article_plotting_funcs.compute_corr_comp(data_7seq)
    article_plotting_funcs.heatmap_avg_subj(pearsonr, gfp_data['times'] * 1000, xlims=[-50, 350], ylims=[-0.5, 0.5], fig_name=op.join(config.fig_path, 'GFP', 'GFPxComplexity_Habituation_heatmap_complexity_pearsonr.png'), figsize=(10, 0.5))
    article_plotting_funcs.heatmap_avg_subj(pearsonr, gfp_data['times'] * 1000, xlims=[-50, 350], ylims=[-0.5, 0.5], fig_name=op.join(config.fig_path, 'GFP', 'GFPxComplexity_Habituation_heatmap_complexity_pearsonr.svg'), figsize=(10, 0.5))

    #  ============== STAND PLOTS ============== #
    data_7seq = np.dstack(gfp_data['standard']['mag'].values())
    data_7seq = np.transpose(data_7seq, (2, 0, 1))
    # Data line plot 7seq
    article_plotting_funcs.plot_7seq_timecourses(data_7seq, gfp_data['times'] * 1000, save_fig_path='GFP/', fig_name='GFPxComplexity_Standard', suffix='',
                                                 pos_horizontal_bar=0.47, plot_pearson_corrComplexity=True, chance=None, xlims=[-50, 350], ymin=0, ylabel='GFP', logger=f)
    # Correlation with complexity heatmap
    pearsonr = article_plotting_funcs.compute_corr_comp(data_7seq)
    article_plotting_funcs.heatmap_avg_subj(pearsonr, gfp_data['times'] * 1000, xlims=[-50, 350], ylims=[-0.5, 0.5], fig_name=op.join(config.fig_path, 'GFP', 'GFPxComplexity_Standard_heatmap_complexity_pearsonr.png'), figsize=(10, 0.5))
    article_plotting_funcs.heatmap_avg_subj(pearsonr, gfp_data['times'] * 1000, xlims=[-50, 350], ylims=[-0.5, 0.5], fig_name=op.join(config.fig_path, 'GFP', 'GFPxComplexity_Standard_heatmap_complexity_pearsonr.svg'), figsize=(10, 0.5))

    #  ============== DEV PLOTS ============== #
    data_7seq = np.dstack(gfp_data['violation']['mag'].values())
    data_7seq = np.transpose(data_7seq, (2, 0, 1))
    # Data line plot 7seq
    article_plotting_funcs.plot_7seq_timecourses(data_7seq, gfp_data['times'] * 1000, save_fig_path='GFP/', fig_name='GFPxComplexity_Deviant', suffix='',
                                                 pos_horizontal_bar=0.47, plot_pearson_corrComplexity=True, chance=None, xlims=[-50, 600], ymin=0, ylabel='GFP', logger=f)
    # Correlation with complexity heatmap
    pearsonr = article_plotting_funcs.compute_corr_comp(data_7seq)
    article_plotting_funcs.heatmap_avg_subj(pearsonr, gfp_data['times'] * 1000, xlims=[-50, 600], ylims=[-0.5, 0.5], fig_name=op.join(config.fig_path, 'GFP', 'GFPxComplexity_viol_heatmap_complexity_pearsonr.png'), figsize=(10, 1))
    article_plotting_funcs.heatmap_avg_subj(pearsonr, gfp_data['times'] * 1000, xlims=[-50, 600], ylims=[-0.5, 0.5], fig_name=op.join(config.fig_path, 'GFP', 'GFPxComplexity_viol_heatmap_complexity_pearsonr.svg'), figsize=(10, 0.5))

    f.close()

# run !
compute_GFP_asmag_allparticipants()
plot_and_save_GFP_stats()