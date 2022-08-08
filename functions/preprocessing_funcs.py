import os.path as op

import mne
import numpy as np
from warnings import warn
import os.path as op
import mne
from mne.report import Report
from mne.preprocessing import ICA
import config
import os.path as op

import mne
from mne.parallel import parallel_func
from mne.preprocessing import read_ica
from mne.preprocessing import create_eog_epochs, create_ecg_epochs
# from mne.report import Report
from mne import open_report
from mne.preprocessing import read_ica

import numpy as np
import config


#_______________________________________________________________________________________________________________________
def run_filter(subject):
    """
    ===========================
    01. Filter using MNE-python
    ===========================

    The data are bandpass filtered to the frequencies defined in config.py
    (config.h_freq - config.l_freq Hz) using linear-phase fir filter with
    delay compensation.
    The transition bandwidth is automatically defined. See
    `Background information on filtering
    <http://mne-tools.github.io/dev/auto_tutorials/plot_background_filtering.html>`_
    for more. The filtered data are saved to separate files to the subject's'MEG'
    directory.
    If config.plot = True plots raw data and power spectral density.
    """

    print("Processing subject: %s" % subject)

    meg_subject_dir = op.join(config.meg_dir, subject)

    n_raws = 0
    for run in config.runs_dict[subject]:

        # read bad channels for run from config
        if run:
            bads = config.bads[subject][run]
        else:
            bads = config.bads[subject]

        extension = run + '_raw'
        raw_fname_in = op.join(meg_subject_dir,
                               config.base_fname.format(**locals()))

        extension = run + '_filt_raw'
        raw_fname_out = op.join(meg_subject_dir,
                                config.base_fname.format(**locals()))

        print("Input: ", raw_fname_in)
        print("Output: ", raw_fname_out)

        if not op.exists(raw_fname_in):
            warn('Run %s not found for subject %s ' %
                 (raw_fname_in, subject))
            # continue  ## SyntaxError: 'continue' not properly in loop

        raw = mne.io.read_raw_fif(raw_fname_in,
                                  allow_maxshield=config.allow_maxshield,
                                  preload=True, verbose='error')

        # ============ NECESSARY? ===========#
        raw.set_eeg_reference('average', projection=True)  # set EEG average reference
        # ===================================#

        # ============ New ===========#
        print("Warning: notch filter 50Hz")
        raw.notch_filter(np.arange(50, 250, 50), fir_design='firwin')
        # ===================================#

        # add bad channels
        raw.info['bads'] = bads
        print("added bads: ", raw.info['bads'])

        if config.set_channel_types is not None:
            raw.set_channel_types(config.set_channel_types)
        if config.rename_channels is not None:
            raw.rename_channels(config.rename_channels)

        # # Interpolating bad channels // DONE BY MAXFILTER FOR MEG, AFTER MAGFILTER FOR EEG
        # print("WARNING - interpolating bad channels: ")
        # print(*bads, sep=", ")
        # raw.interpolate_bads()

        # Band-pass the data channels (MEG and EEG)
        print("Filtering data between %s and %s (Hz)" %
              (config.l_freq, config.h_freq))
        raw.filter(
            config.l_freq, config.h_freq,
            l_trans_bandwidth=config.l_trans_bandwidth,
            h_trans_bandwidth=config.h_trans_bandwidth,
            filter_length='auto', phase='zero', fir_window='hamming',
            fir_design='firwin')

        # ============ AVOID RESAMPLING BEFORE EPOCHING? ===========#
        # if config.resample_sfreq:
        #     print("Resampling data to %.1f Hz" % config.resample_sfreq)
        #     raw.resample(config.resample_sfreq, npad='auto')
        # ==========================================================#

        raw.save(raw_fname_out, overwrite=True)
        n_raws += 1

        config.plot = False
        if config.plot:
            # plot raw data
            raw.plot(n_channels=50, butterfly=True, group_by='position')

            # plot power spectral density
            raw.plot_psd(area_mode='range', tmin=10.0, tmax=100.0,
                         fmin=0., fmax=50., average=True)

    if n_raws == 0:
        raise ValueError('No input raw data found.')

#_______________________________________________________________________________________________________________________
def run_maxwell_filter(subject):
    """
    ===================================
    03. Maxwell filter using MNE-Python
    ===================================

    The data are Maxwell filtered using SSS or tSSS (if config.mf_st_duration
    is not None) and movement compensation.

    Using tSSS with a short duration can be used as an alternative to highpass
    filtering. For instance, a duration of 10 s acts like a 0.1 Hz highpass.

    The head position of all runs is corrected to the run specified in
    config.mf_reference_run.
    It is critical to mark bad channels before Maxwell filtering.

    The function loads machine-specific calibration files from the paths set for
    config.mf_ctc_fname  and config.mf_cal_fname.
    """
    print("Processing subject: %s" % subject)

    meg_subject_dir = op.join(config.meg_dir, subject)

    # To match their processing, transform to the head position of the
    # defined run
    extension = config.runs[config.mf_reference_run] + '_filt_raw'
    raw_fname_in = op.join(meg_subject_dir,
                           config.base_fname.format(**locals()))
    info = mne.io.read_info(raw_fname_in)
    destination = info['dev_head_t']

    for run in config.runs_dict[subject]:

        extension = run + '_filt_raw'
        raw_fname_in = op.join(meg_subject_dir,
                               config.base_fname.format(**locals()))

        extension = run + '_sss_raw'
        raw_fname_out = op.join(meg_subject_dir,
                                config.base_fname.format(**locals()))

        print("Input: ", raw_fname_in)
        print("Output: ", raw_fname_out)

        raw = mne.io.read_raw_fif(raw_fname_in, allow_maxshield=True)
        # Fix coil types (does something only if needed). See:
        # https://martinos.org/mne/stable/generated/mne.channels.fix_mag_coil_types.html  # noqa
        raw.fix_mag_coil_types()

        if config.mf_st_duration:
            print('    st_duration=%d' % (config.mf_st_duration,))

        raw_sss = mne.preprocessing.maxwell_filter(
            raw,
            calibration=config.mf_cal_fname,
            cross_talk=config.mf_ctc_fname,
            st_duration=config.mf_st_duration,
            origin=config.mf_head_origin,
            destination=destination)

        # Interpolating bad EEG channels (bad MEG channels were already recontructed by Maxfilter)
        print("Interpolating bad EEG channels: ")
        print(raw_sss.info['bads'])
        raw_sss.interpolate_bads(reset_bads=True, method=dict(eeg="spline"))

        # Save
        raw_sss.save(raw_fname_out, overwrite=True)

        config.plot = False
        if config.plot:
            # plot maxfiltered data
            raw_sss.plot(n_channels=50, butterfly=True, group_by='position')

#_______________________________________________________________________________________________________________________
def run_ica(subject, tsss=config.mf_st_duration):
    """
    ===========
    03. Run ICA
    ===========
    This fits ICA on the last 4 blocks of the experiment data filtered with 1 Hz highpass,
    for this purpose only using fastICA. Separate ICAs are fitted and stored for
    MEG and EEG data.
    To actually remove designated ICA components from your data, you will have to
    1- first run 04-identify_EOG_ECG_components_ica.py to automatically identify the components related to the EOG and ECG artefacts.
    2- Inspecting the report, confirm or correct the proposed components and mark them in config.rejcomps_man
    3 - Only once you did so, run 05-apply_ica.py
    """


    print("Processing subject: %s" % subject)

    meg_subject_dir = op.join(config.meg_dir, subject)

    raw_list = list()
    print("  Loading raw data")
    runs = config.runs_dict[subject]

    for run in runs[-4:-1]: # load four last runs
        if config.use_maxwell_filter:
            extension = run + '_sss_raw'
        else:
            extension = run + '_filt_raw'

        raw_fname_in = op.join(meg_subject_dir,
                               config.base_fname.format(**locals()))

        raw = mne.io.read_raw_fif(raw_fname_in, preload=True)
        raw_list.append(raw)

    print('  Concatenating runs')
    raw = mne.concatenate_raws(raw_list)
    if "eeg" in config.ch_types:
        raw.set_eeg_reference(projection=True)
    del raw_list


    # don't reject based on EOG to keep blink artifacts
    # in the ICA computation.
    reject_ica = config.reject
    if reject_ica and 'eog' in reject_ica:
        reject_ica = dict(reject_ica)
        del reject_ica['eog']

    # produce high-pass filtered version of the data for ICA
    raw_ica = raw.copy().filter(l_freq=1., h_freq=None)

    print("  Running ICA...")

    picks_meg = mne.pick_types(raw.info, meg=True, eeg=False,
                               eog=False, stim=False, exclude='bads')
    picks_eeg = mne.pick_types(raw.info, meg=False, eeg=True,
                               eog=False, stim=False, exclude='bads')
    all_picks = {'meg': picks_meg, 'eeg': picks_eeg}

    n_components = {'meg': 0.999, 'eeg': 0.999}

    ch_types = []
    if 'eeg' in config.ch_types:
        ch_types.append('eeg')
    if set(config.ch_types).intersection(('meg', 'grad', 'mag')):
        ch_types.append('meg')

    for ch_type in ch_types:
        print('Running ICA for ' + ch_type)

        ica = ICA(method='fastica', random_state=config.random_state,
                  n_components=n_components[ch_type])

        picks = all_picks[ch_type]

        ica.fit(raw, picks=picks, decim=config.ica_decim)

        print('  Fit %d components (explaining at least %0.1f%% of the'
              ' variance)' % (ica.n_components_, 100 * n_components[ch_type]))

        ica_fname = \
            '{0}_{1}_{2}-ica.fif'.format(subject, config.study_name, ch_type)
        ica_fname = op.join(meg_subject_dir, ica_fname)
        ica.save(ica_fname)

        # if config.plot:

        # plot ICA components to html report
        report_fname = \
            '{0}_{1}_{2}-ica.h5'.format(subject, config.study_name,
                                          ch_type)
        report_fname = op.join(meg_subject_dir, report_fname)
        report_fname_html = \
            '{0}_{1}_{2}-ica.html'.format(subject, config.study_name,
                                          ch_type)
        report_fname = op.join(meg_subject_dir, report_fname)
        report = Report(report_fname, verbose=False)

        for idx in range(0, ica.n_components_):
            figure = ica.plot_properties(raw,
                                         picks=idx,
                                         psd_args={'fmax': 60},
                                         show=False)

            report.add_figs_to_section(figure, section=subject,
                                       captions=(ch_type.upper() +
                                                 ' - ICA Components'))

        report.save(report_fname, overwrite=True)
        report.save(report_fname_html, overwrite=True)


#_______________________________________________________________________________________________________________________
def automatic_identification_of_components(subject):
    """
    ===============
    04- identify_EOG_ECG_components_ica.py
    ===============

    Blinks and ECG artifacts are automatically detected and the corresponding ICA
    components are removed from the data.
    This relies on the ICAs computed in 03-run_ica.py
    !! If you manually add components to remove (config.rejcomps_man),
    make sure you did not re-run the ICA in the meantime. Otherwise (especially if
    the random state was not set, or you used a different machine, the component
    order might differ).

    !! Inspect the .html report, confirm or correct the proposed components and mark them in config.rejcomps_man
    """
    if subject in config.rejcomps_man:
        print(subject)
        raise Exception(
            'The EOG and ECG components were already identified and hand written in the config.rejcomps_man\n Delete it if you want to rerun this part of the script.')
    else:
        ica = dict(meg=[], eeg=[])
        ica_reject = dict(meg=[], eeg=[])

        print("Identifying the components for subject: %s" % subject)
        meg_subject_dir = op.join(config.meg_dir, subject)

        # ==================================================
        # concatenate the four last runs to compute the correlation
        # between the ICA components determined in 03-run_ica with the epochs
        # on ECG and EOG events
        # ==================================================
        runs = config.runs_dict[subject]


        raw_list = list()
        for run in runs[-4:-1]: # load four last runs
            if config.use_maxwell_filter:
                extension = run + '_sss_raw'
            else:
                extension = run + '_filt_raw'

            raw_fname_in = op.join(meg_subject_dir,
                                   config.base_fname.format(**locals()))
            raw = mne.io.read_raw_fif(raw_fname_in, preload=True)
            raw_list.append(raw)

        print('  Concatenating runs')
        raw_ref_runs = mne.concatenate_raws(raw_list)
        del raw

        # ==================================================
        # define the channels corresponding to meg and eeg data
        # ==================================================

        picks_meg = mne.pick_types(raw_ref_runs.info, meg=True, eeg=False,
                                   eog=False, stim=False, exclude='bads')
        picks_eeg = mne.pick_types(raw_ref_runs.info, meg=False, eeg=True,
                                   eog=False, stim=False, exclude='bads')
        all_picks = {'meg': picks_meg, 'eeg': picks_eeg}

        print('Finding ICA components correlating with ECG and EOG epochs...')

        # ==================================================
        # define ch_types according to the recorded modalities
        # ==================================================

        ch_types = []
        if 'eeg' in config.ch_types:
            ch_types.append('eeg')
        if set(config.ch_types).intersection(('meg', 'grad', 'mag')):
            ch_types.append('meg')


        # ==================================================
        # main loop
        # ==================================================
        for ch_type in ch_types:
            print(ch_type)
            picks = all_picks[ch_type]

            # ==================================================
            # Load ICA
            # ==================================================

            fname_ica = op.join(meg_subject_dir,
                                '{0}_{1}_{2}-ica.fif'.format(subject,
                                                             config.study_name,
                                                             ch_type))
            print('Reading ICA: ' + fname_ica)
            ica[ch_type] = read_ica(fname=fname_ica)

            # ==================================================
            # report input and output names
            # ==================================================
            # Load previous report (.h5)
            report_fname = \
                '{0}_{1}_{2}-ica.h5'.format(subject,
                                                     config.study_name,
                                                     ch_type)
            report_fname = op.join(meg_subject_dir, report_fname)

            # set the name of the final report
            report_fname_html = \
                '{0}_{1}_{2}-ica.html'.format(subject,
                                                     config.study_name,
                                                     ch_type)
            report_fname_html = op.join(meg_subject_dir, report_fname_html)
            report= open_report(report_fname)


            # ==================================================
            # Correlation with ECG epochs
            # ==================================================

            pick_ecg = mne.pick_types(raw_ref_runs.info, meg=False, eeg=False,
                                      ecg=True, eog=False)

            # either needs an ecg channel, or avg of the mags (i.e. MEG data)
            if pick_ecg or ch_type == 'meg':

                picks_ecg = np.concatenate([picks, pick_ecg])

                # Create ecg epochs
                if ch_type == 'meg':
                    reject = {'mag': config.reject['mag'],
                              'grad': config.reject['grad']}
                elif ch_type == 'eeg':
                    reject = {'eeg': config.reject['eeg']}

                ecg_epochs = create_ecg_epochs(raw_ref_runs, picks=picks_ecg, reject=None,
                                               baseline=(None, 0), tmin=-0.5,
                                               tmax=0.5)
                ecg_average = ecg_epochs.average()

                ecg_inds, scores = \
                    ica[ch_type].find_bads_ecg(ecg_epochs, method='ctps',
                                      threshold=config.ica_ctps_ecg_threshold)
                del ecg_epochs
                # params = dict(exclude=ecg_inds, show=config.plot)
                params = dict(show=config.plot)  # The exclude parameter is deprecated and will be removed in version 0.20; specify excluded components using the ICA.exclude attribute instead.
                ica[ch_type].exclude = ecg_inds

                # == == == == == == == ==  plots appended to report = == == == == == == == == == == ==
                # Plot r score
                report.add_figs_to_section(ica[ch_type].plot_scores(scores,**params),
                                           captions=ch_type.upper() + ' - ECG - ' +
                                           'R scores')
                # Plot source time course
                report.add_figs_to_section(ica[ch_type].plot_sources(ecg_average,**params),
                                           captions=ch_type.upper() + ' - ECG - ' +
                                           'Sources time course')
                # Plot source time course
                report.add_figs_to_section(ica[ch_type].plot_overlay(ecg_average,**params),
                                           captions=ch_type.upper() + ' - ECG - ' +
                                           'Corrections')


            else:
                # XXX : to check when EEG only is processed
                print('no ECG channel is present. Cannot automate ICAs component '
                      'detection for EOG!')


            # ==================================================
            # Correlation with EOG epochs
            # ==================================================
            pick_eog = mne.pick_types(raw_ref_runs.info, meg=False, eeg=False,
                                      ecg=False, eog=True)

            if pick_eog.any():
                print('using EOG channel')
                picks_eog = np.concatenate([picks, pick_eog])
                # Create eog epochs
                eog_epochs = create_eog_epochs(raw_ref_runs, picks=picks_eog, reject=None,
                                               baseline=(None, 0), tmin=-0.5,
                                               tmax=0.5)

                eog_average = eog_epochs.average()
                eog_inds, scores = ica[ch_type].find_bads_eog(eog_epochs, threshold=3.0)
                del eog_epochs
                # params = dict(exclude=eog_inds, show=config.plot)
                params = dict(show=config.plot)  # The exclude parameter is deprecated and will be removed in version 0.20; specify excluded components using the ICA.exclude attribute instead.
                ica[ch_type].exclude = eog_inds

                # == == == == == == == ==  plots appended to report = == == == == == == == == == == ==
                # Plot r score
                report.add_figs_to_section(ica[ch_type].plot_scores(scores, **params),
                                           captions=ch_type.upper() + ' - EOG - ' +
                                           'R scores')

                # Plot source time course
                report.add_figs_to_section(ica[ch_type].plot_sources(eog_average, **params),
                                           captions=ch_type.upper() + ' - EOG - ' +
                                           'Sources time course')

                # Plot source time course
                report.add_figs_to_section(ica[ch_type].plot_overlay(eog_average, **params),
                                           captions=ch_type.upper() + ' - EOG - ' +
                                           'Corrections')

                report.save(report_fname, overwrite=True)
                report.save(report_fname_html, overwrite=True, open_browser=False)

            else:
                print('no EOG channel is present. Cannot automate ICAs component '
                      'detection for EOG!')

            ica_reject[ch_type] = list(ecg_inds) + list(eog_inds)

            # now reject the components
            print('Rejecting from %s: %s' % (ch_type, ica_reject))

        # ==================================================
        #  Visualize the data before and after cleaning
        # ==================================================
        fig = ica[ch_type].plot_overlay(raw_ref_runs, exclude=ica_reject[ch_type], show=config.plot)
        report.add_figs_to_section(fig, captions=ch_type.upper() +
                                                 ' - ALL(epochs) - Corrections')

        report.save(report_fname_html, overwrite=True, open_browser=False)
        report.save(report_fname, overwrite=True)


#_______________________________________________________________________________________________________________________
def apply_ica(subject):
    """
    ===============
    05. Apply ICA
    ===============

    This function loads the ica filter as well as the user confirmed/corrected set of components to reject (in order to remove the artefacts).
    It applies the filter and saves the ICA-filtered data

    """
    if subject not in config.rejcomps_man:
        raise Exception(
            'The EOG and ECG components were not saved in config.rejcomps_man\n You have to run 03-run-ica.py and 04-indentify_EOG_ECG_components.py first. Then manually enter the components in config.rejcomps_man.')


    print("Applying the ICA for subject: %s" % subject)
    meg_subject_dir = op.join(config.meg_dir, subject)

    # ==================================================
    # define ch_types according to the recorded modalities
    # ==================================================

    ch_types = []
    if 'eeg' in config.ch_types:
        ch_types.append('eeg')
    if set(config.ch_types).intersection(('meg', 'grad', 'mag')):
        ch_types.append('meg')


    ica = dict(meg=[], eeg=[])
    ica_reject = dict(meg=[], eeg=[])


    for ch_type in ch_types:
        print(ch_type)

        # ==================================================
        # Load ICA
        # ==================================================

        fname_ica = op.join(meg_subject_dir,
                            '{0}_{1}_{2}-ica.fif'.format(subject,
                                                         config.study_name,
                                                         ch_type))
        print('Reading ICA: ' + fname_ica)
        ica[ch_type] = read_ica(fname=fname_ica)

        ica_reject[ch_type] = list(config.rejcomps_man[subject][ch_type])
        print('Using user-defined bad ICA components')



    for run in config.runs_dict[subject]:
        if config.use_maxwell_filter:
            extension = run + '_sss_raw'
        else:
            extension = run + '_filt_raw'

        # = load the ICA =

        raw_fname_in = op.join(meg_subject_dir,
                               config.base_fname.format(**locals()))

        raw_before = mne.io.read_raw_fif(raw_fname_in, preload=True)

        raw_ica_eeg = ica['eeg'].apply(raw_before, exclude=ica_reject['eeg'])
        raw_ica_meg_eeg = ica['meg'].apply(raw_ica_eeg, exclude=ica_reject['meg'])
        extension = run + '_ica_raw'
        fname_out = op.join(meg_subject_dir,
                               config.base_fname.format(**locals()))

        print('Saving cleaned runs')
        raw_ica_meg_eeg.save(fname_out, overwrite=True)