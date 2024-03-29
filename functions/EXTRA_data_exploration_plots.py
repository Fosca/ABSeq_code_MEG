from ABseq_func import *


def make_figures(subject):
    # ----------------------------------------------------------------------------------------------------------- #
    # PLOTS - Items epochs
    # ----------------------------------------------------------------------------------------------------------- #

    # # -- LOAD THE EPOCHS -- #
    epochs = epoching_funcs.load_epochs_items(subject, cleaned=True, AR_type='local')
    #
    # # -- PLOTS -- #
    evoked_funcs.plot_butterfly_items(epochs, subject, violation_or_not=1, apply_baseline=True)  # items epoch are not already baseline corrected
    evoked_funcs.plot_butterfly_items(epochs, subject, violation_or_not=0, apply_baseline=True)
    GFP_funcs.plot_gfp_items_standard_or_deviants(epochs, subject, h_freq=None, standard_or_deviant='standard')
    GFP_funcs.plot_gfp_items_standard_or_deviants(epochs, subject, h_freq=None, standard_or_deviant='deviant')

    # ----------------------------------------------------------------------------------------------------------- #
    # PLOTS - Full-sequence epochs
    # ----------------------------------------------------------------------------------------------------------- #
    # -- LOAD THE EPOCHS -- #
    epochs = epoching_funcs.load_epochs_full_sequence(subject, cleaned=True, AR_type='global')

    # -- PLOTS -- #
    evoked_funcs.plot_butterfly_first_item(epochs, subject, apply_baseline=False, ch_types = ['grad', 'mag'])  # fullseq epoch are already baseline corrected
    GFP_funcs.plot_gfp_full_sequence_standard(epochs, subject, h_freq=None, ch_types = ['grad', 'mag'])
    GFP_funcs.plot_gfp_full_sequence_deviants_4pos(epochs, subject, h_freq=None, ch_types = ['grad', 'mag'])

