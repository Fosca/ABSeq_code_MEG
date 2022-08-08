#!/bin/bash
        #PBS -N BIDS_sub01-pa_190002.py
        #PBS -q Nspin_run32 
        #PBS -l walltime=48:00:00
        #PBS -l nodes=1:ppn=1
        #PBS -o /neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/ABSeq_PUBLICATION/cluster//results_qsub/BIDS/std_BIDS_sub01-pa_190002.py
        #PBS -e /neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/ABSeq_PUBLICATION/cluster//results_qsub/BIDS/err_BIDS_sub01-pa_190002.py 
        cd /neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/ABSeq_PUBLICATION/cluster//results_qsub/
        python /neurospin/meg/meg_tmp/ABSeq_Samuel_Fosca2019/scripts/ABSeq_PUBLICATION/cluster//generated_jobs/BIDS/BIDS_sub01-pa_190002.py