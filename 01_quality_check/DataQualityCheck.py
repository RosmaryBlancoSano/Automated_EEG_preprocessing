"""DataQualityCheck.py -- automatic quality check of the raw Stroke recordings (EDF, 8 channels).

For every recording: FIR notch filters at 25, 50 and 100 Hz, 1 Hz high-pass and 35 Hz low-pass filters; detection of the high-amplitude segments
with the Global Variance of the Temporal Derivative (GVTD; threshold = mean + 10 left-sided standard deviations of the GVTD trace, +-0.5 s around
every sample above the threshold) and rejection of those segments; bad channels with the PREP deviation criterion (deviation threshold 2) and
flat channels. A recording is classified as BAD when at least 25% of the channels are bad or at least 50% of the data were rejected, GOOD otherwise.
Bad channels are only reported here: they are interpolated afterwards (06_real_datasets/Dataset_Stroke_new.m).

Usage
    python DataQualityCheck.py RAW_DIR [--subjects-csv FILE] [--plot]      RAW_DIR/sub-*/eeg/*.edf
    --subjects-csv: table of the participants (column participant_id) to select the recordings and add its columns to the summary
    --plot: histogram of the GVTD with its threshold for every recording

Output
    RAW_DIR/GoodQuality/<subject>/ and RAW_DIR/BadQuality/<subject>/  (EEGLAB .set files)
    RAW_DIR/EEG_quality_summary.csv  (channels and share of data rejected, quality of every recording)

Local modules in this folder: Functions.py (NoisyChannels, GVTD_EEG, eeg_reject), GVTD.py, ReadEDF32/ (EDF reader adapted from MNE).
The recordings are patient data and cannot be shared.
"""
import argparse
import os
import re

import mne
import numpy as np
import pandas as pd

import GVTD
from Functions import NoisyChannels, GVTD_EEG, logger, eeg_reject

NOTCH = [25, 50, 100]      # Hz
HIGHPASS = 1               # Hz
LOWPASS = 35               # Hz
GVTD_NSTD = 10             # threshold of the GVTD: mean + GVTD_NSTD left-sided standard deviations
DEVIATION_THRESHOLD = 2    # PREP deviation criterion
FLAT_THRESHOLD = 1e-6      # V
MAX_BAD_CHANNELS = 0.25    # a recording is BAD with at least this share of bad channels ...
MAX_BAD_DURATION = 50      # ... or at least this percentage of rejected data


def quality_check(data, sfreq, ch_names, plot=False):
    """data: channels x samples in volts. Returns the preprocessed (filtered, segments rejected) mne Raw and the quality measures."""
    info = mne.create_info(ch_names, sfreq, 'eeg')
    raw = mne.io.RawArray(data, info)
    filtered = mne.filter.notch_filter(raw.get_data(), Fs=sfreq, phase='zero', freqs=NOTCH, method='fir', mt_bandwidth=4, p_value=0.01,
                                       filter_length='20s')
    raw_hf = mne.io.RawArray(filtered, raw.info).filter(l_freq=HIGHPASS, h_freq=None)
    raw_lf = raw_hf.copy().filter(l_freq=None, h_freq=LOWPASS, method='fir', phase='zero', fir_design='firwin', skip_by_annotation='edge')

    # high-amplitude segments: GVTD above the threshold, +-0.5 s around each sample
    signal = raw_lf.get_data()
    gvtd_trace = GVTD.gvtd_processing(signal)
    if plot:
        GVTD.make_gvtd_hist(gvtd_trace, plotThresh=True, statType='Mean', nStd=GVTD_NSTD)
    thresh = GVTD.find_gvtd_thresh(gvtd_trace, statType='Mean', nStd=GVTD_NSTD)
    raw_gvtd, info, bad_segments = GVTD_EEG(signal, gvtd_trace, thresh, sfreq, raw_lf)
    raw_gvtd.info['bads'] = []
    rejected = eeg_reject(raw_gvtd.get_data(), bad_segments)
    cleaned = mne.io.RawArray(rejected, raw_gvtd.info)

    # bad channels on the data left after the rejection
    noisy = NoisyChannels(cleaned, random_state=1337)
    noisy.find_bad_ch_by_deviation(deviation_threshold=DEVIATION_THRESHOLD)
    noisy.find_bad_by_nan_flat(flat_threshold=FLAT_THRESHOLD)
    bads = noisy.get_bads(verbose=True)

    total = len(raw_gvtd.times) / sfreq
    kept = rejected.shape[1] / sfreq
    percent_bad_duration = (total - kept) / total * 100
    percent_bad_channels = len(bads) / len(ch_names) * 100
    bad_quality = len(bads) >= MAX_BAD_CHANNELS * len(ch_names) or percent_bad_duration >= MAX_BAD_DURATION
    logger.info(f"\nNoisyChannels REPORT\n------------------------\n{np.round(percent_bad_channels, 2)}% of the channels were detected as noisy."
                f"\n({len(bads)}) channels: {', '.join(bads)} Percentage of rejected data: {percent_bad_duration:.3f}%")
    quality = {'eeg_quality': 'bad' if bad_quality else 'good', 'n_bad_channels': len(bads), 'bad_channels': ', '.join(bads),
               'flat_channels': ', '.join(noisy.bad_by_flat), 'noisy_channels': ', '.join(noisy.bad_by_deviation),
               'total_duration': total, 'duration_after_rejection': kept, 'percent_bad_duration': percent_bad_duration}
    return cleaned, quality


# ------------------------------------------------------------------------------------------------------------------- recordings (EDF)
def stroke_recordings(root, subjects_csv=None):
    """Yields (subject, edf file, metadata) for RAW_DIR/sub-*/eeg/*.edf, in numerical order of the subjects."""
    meta = None
    if subjects_csv:
        info = pd.read_csv(subjects_csv)
        info['participant_id'] = info['participant_id'].str.strip().str.lower()
        meta = info.set_index('participant_id').to_dict('index')
    key = lambda s: int(re.search(r'\d+', s).group()) if re.search(r'\d+', s) else float('inf')
    for subject in sorted(os.listdir(root), key=key):
        if not (os.path.isdir(os.path.join(root, subject)) and subject.startswith('sub-')):
            continue
        if meta is not None and subject not in meta:
            continue
        for folder, _, files in os.walk(os.path.join(root, subject, 'eeg')):
            for f in files:
                if f.endswith('.edf'):
                    yield subject, os.path.join(folder, f), (meta or {}).get(subject, {})


def load_stroke(edf):
    from ReadEDF32 import read_raw_edf32      # EDF reader adapted from MNE, needed only for the Stroke recordings
    raw = read_raw_edf32(edf, preload=True)
    raw.rename_channels({ch: ch.replace('-Ref1', '') for ch in raw.ch_names})
    data = raw.get_data()
    if np.max(np.abs(data)) > 100:            # the data are in microvolts
        data = data * 1e-6
    return data, raw.info['sfreq'], raw.ch_names


def save_stroke(root, subject, edf, cleaned, quality):
    folder = os.path.join(root, 'GoodQuality' if quality['eeg_quality'] == 'good' else 'BadQuality', subject)
    os.makedirs(folder, exist_ok=True)
    out = os.path.join(folder, f"{subject}_{os.path.splitext(os.path.basename(edf))[0]}.set")
    mne.export.export_raw(out, cleaned, fmt='eeglab', overwrite=True)
    return out


# ------------------------------------------------------------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description='Automatic quality check of the raw Stroke recordings.')
    parser.add_argument('folder', help='raw data folder: folder/sub-*/eeg/*.edf')
    parser.add_argument('--subjects-csv', help='table of the participants (column participant_id) to select the recordings and add its columns to the summary')
    parser.add_argument('--plot', action='store_true', help='plot the histogram of the GVTD of every recording')
    args = parser.parse_args()

    results = []
    for subject, path, meta in stroke_recordings(args.folder, args.subjects_csv):
        try:
            data, sfreq, ch_names = load_stroke(path)
            cleaned, quality = quality_check(data, sfreq, ch_names, plot=args.plot)
            out = save_stroke(args.folder, subject, path, cleaned, quality)
            print(f"{subject}: {quality['eeg_quality']} -> {out}")
            results.append({'participant_id': subject, **meta, **quality})
        except Exception as e:
            print(f'Error processing {path}: {e}')
    summary = os.path.join(args.folder, 'EEG_quality_summary.csv')
    pd.DataFrame(results).to_csv(summary, index=False)
    print(f'\nEEG quality summary saved to {summary}')


if __name__ == '__main__':
    main()
