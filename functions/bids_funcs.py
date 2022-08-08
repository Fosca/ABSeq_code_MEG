import sys
sys.path.append('/neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/ABSeq_scripts')
sys.path.append('/neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/')
import mne
import mne_bids
import config
from functions import epoching_funcs
import os
import numpy as np

def from_metadata_to_events_bids(subject,run):

    fields_of_interest = ['SequenceID','RunNumber','StimPosition','StimID','ViolationOrNot','Complexity','ViolationInSequence']

    epochs = epoching_funcs.load_epochs_items(subject, cleaned=False)
    epochs_AR = epoching_funcs.load_epochs_items(subject, cleaned=True)
    epochs.metadata['drop_log'] = epochs_AR.drop_log
    epochs_run = epochs['RunNumber == %i'%int(run)]
    events = epochs_run.events
    event_id = dict()

    for ii in range(len(events)):
        meta_ii = epochs_run[ii].metadata
        event_value = int(np.sum([meta_ii[fields_of_interest[k]].values[0]*10**(len(fields_of_interest)-k) for k in range(len(fields_of_interest))]))
        key_value = '/'.join([fields_of_interest[k]+'_'+str(int(meta_ii[fields_of_interest[k]].values[0])) for k in range(len(fields_of_interest))])
        events[ii,2] = event_value
        event_id[key_value] = event_value

    return events, event_id

def convert_to_bids(subject):

    input_dir = config.study_path
    session = ''
    task = 'BinaSeq'

    for run in config.runs_dict[subject]:
        run = run[3:5]
        bids_root = os.path.join(config.study_path, 'BIDS_good') + os.path.sep

        bids_path = mne_bids.BIDSPath(
            subject=subject[3:5],
            session=session,
            run=run,
            task=task,
            datatype='meg',
            suffix='meg',
            extension='.fif',
            root=bids_root
        )

        # Load raw files and find back the events from the epochs, it is easier

        raw = mne.io.read_raw(f'{input_dir}/MEG/{subject}/run{run}_raw.fif', allow_maxshield='yes')
        raw.info['bads'] = config.bads[subject][f'run{run}']

        events, event_id = from_metadata_to_events_bids(subject, run)

        mne_bids.write_raw_bids(
            raw=raw,
            bids_path=bids_path,
            events_data=events,
            event_id=event_id,
            overwrite=True,
            symlink=False,
            anonymize={'daysback': 40000},
            verbose=True
        )

        # write MEG calibration files
        cal_fname = config.study_path + '/system_calibration_files/sss_cal_nspn.dat'
        ct_fname = config.study_path + '/system_calibration_files/ct_sparse_nspn.fif'

        mne_bids.write_meg_calibration(
            calibration=cal_fname,
            bids_path=bids_path
        )
        mne_bids.write_meg_crosstalk(
            fname=ct_fname,
            bids_path=bids_path
        )
        # write anatomical fMRI data

        fname_t1 = f'{input_dir}/MRI/fs_converted/{subject}/mri/T1.mgz'
        if os.path.exists(fname_t1):
            bids_path_t1 = bids_path.copy().update(
                task=task,
                run=run,
                suffix='T1w',
                extension='.nii.gz',
                datatype='anat'
            )
            mne_bids.write_anat(
                image=fname_t1,
                bids_path=bids_path_t1,
                overwrite=True
            )








