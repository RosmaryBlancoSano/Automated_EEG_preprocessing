import mne
import numpy as np
from pyprep.utils import _mat_iqr, _mad
from pyprep.removeTrend import removeTrend
from mne.utils import check_random_state, logger

class NoisyChannels:
    """Detect bad channels in an EEG recording using a range of methods.

    This class provides a number of methods for detecting bad channels across a
    full-session EEG recording. Specifically, this class implements all of the
    noisy channel detection methods used in the PREP pipeline, as described in [1]_.
    The detection methods in this class can be run independently, or can be run
    all at once using the :meth:`~.find_all_bads` method.

    At present, only EEG channels are supported and any non-EEG channels in the
    provided data will be ignored.

    Parameters
    ----------
    raw : mne.io.Raw
        An MNE Raw object to check for bad EEG channels.
    do_detrend : bool, optional
        Whether or not low-frequency (<1.0 Hz) trends should be removed from the
        EEG signal prior to bad channel detection. This should always be set to
        ``True`` unless the signal has already had low-frequency trends removed.
        Defaults to ``True``.
    random_state : {int, None, np.random.RandomState}, optional
        The seed to use for random number generation within RANSAC. This can be
        ``None``, an integer, or a :class:`~numpy.random.RandomState` object.
        If ``None``, a random seed will be obtained from the operating system.
        Defaults to ``None``.
    matlab_strict : bool, optional
        Whether or not PyPREP should strictly follow MATLAB PREP's internal
        math, ignoring any improvements made in PyPREP over the original code
        (see :ref:`matlab-diffs` for more details). Defaults to ``False``.

    References
    ----------
    .. [1] Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., Robbins, K. A.
       (2015). The PREP pipeline: standardized preprocessing for large-scale
       EEG analysis. Frontiers in Neuroinformatics, 9, 16.

    """

    def __init__(self, raw, do_detrend=True, random_state=None, matlab_strict=False):
        # Make sure that we got an MNE object
        assert isinstance(raw, mne.io.BaseRaw)

        raw.load_data()
        self.raw_mne = raw.copy()
        self.raw_mne.pick_types(eeg=True)
        self.sample_rate = raw.info["sfreq"]
        if do_detrend:
            self.raw_mne._data = removeTrend(
                self.raw_mne.get_data(), self.sample_rate, matlab_strict=matlab_strict
            )
        self.matlab_strict = matlab_strict

        # Extra data for debugging
        self._extra_info = {
            "bad_by_deviation": {},
            "bad_by_hf_noise": {},
            "bad_by_correlation": {},
            "bad_by_dropout": {},
            "bad_by_ransac": {},
        }

        # random_state
        self.random_state = check_random_state(random_state)

        # The identified bad channels
        self.bad_by_nan = []
        self.bad_by_flat = []
        self.bad_by_deviation = []
        self.bad_by_hf_noise = []
        self.bad_by_correlation = []
        self.bad_by_SNR = []
        self.bad_by_dropout = []
        self.bad_by_ransac = []

        # Get original EEG channel names, channel count & samples
        ch_names = np.asarray(self.raw_mne.info["ch_names"])
        self.ch_names_original = ch_names
        self.n_chans_original = len(ch_names)
        self.n_samples = raw.get_data().shape[1]

        # Before anything else, flag bad-by-NaNs and bad-by-flats
        self.find_bad_by_nan_flat(flat_threshold=1e-6)
        bads_by_nan_flat = self.bad_by_nan + self.bad_by_flat

        # Make a subset of the data containing only usable EEG channels
        self.usable_idx = np.isin(ch_names, bads_by_nan_flat, invert=True)
        self.EEGData = self.raw_mne.get_data(picks=ch_names[self.usable_idx])
        self.EEGFiltered = None

        # Get usable EEG channel names & channel counts
        self.ch_names_new = np.asarray(ch_names[self.usable_idx])
        self.n_chans_new = len(self.ch_names_new)

    def find_bad_by_nan_flat(self, flat_threshold=1e-15):
        """Detect channels than contain NaN values or have near-flat signals.

        A channel is considered flat if its standard deviation or its median
        absolute deviation from the median (MAD) are below the provided flat
        threshold (default: ``1e-15`` volts).

        This method is run automatically when a ``NoisyChannels`` object is
        initialized, preventing flat or NaN-containing channels from interfering
        with the detection of other types of bad channels.

        Parameters
        ----------
        flat_threshold : float, optional
            The lowest standard deviation or MAD value for a channel to be
            considered bad-by-flat. Defaults to ``1e-15`` volts (corresponds to
            10e-10 µV in MATLAB PREP).
        """
        # print('function flat_threshold =', flat_threshold)
        # Get all EEG channels from original copy of data
        EEGData = self.raw_mne.get_data()

        # Detect channels containing any NaN values
        nan_channel_mask = np.isnan(np.sum(EEGData, axis=1))
        nan_channels = self.ch_names_original[nan_channel_mask]

        # Detect channels with flat or extremely weak signals
        flat_by_mad = _mad(EEGData, axis=1) < flat_threshold
        flat_by_stdev = np.std(EEGData, axis=1) < flat_threshold
        flat_channel_mask = flat_by_mad | flat_by_stdev
        flat_channels = self.ch_names_original[flat_channel_mask]

        # Update names of bad channels by NaN or flat signal
        self.bad_by_nan = nan_channels.tolist()
        self.bad_by_flat = flat_channels.tolist()

    def calculate_empirical_IQR_to_SD(self, data):
        """
        Compute an empirical IQR_TO_SD value for the data, using the ratio between
        IQR and standard deviation to estimate the scaling factor.

        Parameters
        ----------
        data : np.ndarray
            The EEG data (for example, a 2D array in which every row is a channel).

        Returns
        -------
        empirical_IQR_TO_SD : float
            The empirical IQR_TO_SD value.
        """
        # IQR of every channel (row)
        Q1 = np.percentile(data, 25, axis=1)  # 25° percentile (1° quartile)
        Q3 = np.percentile(data, 75, axis=1)  # 75° percentile (3° quartile)
        IQR = Q3 - Q1  # Interquartile Range

        # standard deviation of every channel (row)
        std_dev = np.std(data, axis=1)

        # empirical ratio between IQR and standard deviation
        empirical_IQR_TO_SD = np.median(IQR / std_dev)  # the median reduces the influence of outliers

        return empirical_IQR_TO_SD

    def get_bads(self, verbose=False, as_dict=False):
        """Get the names of all channels currently flagged as bad.

        Note that this method does not perform any bad channel detection itself,
        and only reports channels already detected as bad by other methods.

        Parameters
        ----------
        verbose : bool, optional
            If ``True``, a summary of the channels currently flagged as by bad per
            category is printed. Defaults to ``False``.
        as_dict: bool, optional
            If ``True``, this method will return a dict of the channels currently
            flagged as bad by each individual bad channel type. If ``False``, this
            method will return a list of all unique bad channels detected so far.
            Defaults to ``False``.

        Returns
        -------
        bads : list or dict
            The names of all bad channels detected so far, either as a combined
            list or a dict indicating the channels flagged bad by each type.

        """
        bads = {
            "bad_by_nan": self.bad_by_nan,
            "bad_by_flat": self.bad_by_flat,
            "bad_by_deviation": self.bad_by_deviation,
            "bad_by_hf_noise": self.bad_by_hf_noise,
            "bad_by_correlation": self.bad_by_correlation,
            "bad_by_SNR": self.bad_by_SNR,
            "bad_by_dropout": self.bad_by_dropout,
            "bad_by_ransac": self.bad_by_ransac,
        }

        all_bads = set()
        for bad_chs in bads.values():
            all_bads.update(bad_chs)

        name_map = {"nan": "NaN", "hf_noise": "HF noise", "ransac": "RANSAC"}
        if verbose:
            out = f"Found {len(all_bads)} uniquely bad channels:\n"
            for bad_type, bad_chs in bads.items():
                bad_type = bad_type.replace("bad_by_", "")
                if bad_type in name_map.keys():
                    bad_type = name_map[bad_type]
                out += f"\n{len(bad_chs)} by {bad_type}: {bad_chs}\n"
            logger.info(out)

        if as_dict:
            bads["bad_all"] = list(all_bads)
        else:
            bads = list(all_bads)

        return bads

    def find_bad_ch_by_deviation(self, deviation_threshold=5.0):
        """Detect channels with abnormally high or low overall amplitudes.

        A channel is considered "bad-by-deviation" if its amplitude deviates
        considerably from the median channel amplitude, as calculated using a
        robust Z-scoring method and the given deviation threshold.

        Amplitude Z-scores are calculated using the formula
        ``(channel_amplitude - median_amplitude) / amplitude_sd``, where
        channel amplitudes are calculated using a robust outlier-resistant estimate
        of the signals' standard deviations (IQR scaled to units of SD), and the
        amplitude SD is the IQR-based SD of those amplitudes.

        Parameters
        ----------
        deviation_threshold : float, optional
            The minimum absolute z-score of a channel for it to be considered
            bad-by-deviation. Defaults to ``5.0``.

        """
        # IQR_TO_SD = 0.7413  # Scales units of IQR to units of SD, assuming normality
        # Reference: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/IQR.html

        IQR_TO_SD = self.calculate_empirical_IQR_to_SD(self.EEGData)

        # Get channel amplitudes and the median / robust SD of those amplitudes
        chan_amplitudes = _mat_iqr(self.EEGData, axis=1) * IQR_TO_SD
        amp_sd = _mat_iqr(chan_amplitudes) * IQR_TO_SD
        amp_median = np.nanmedian(chan_amplitudes)

        # Calculate robust Z-scores for the channel amplitudes
        amplitude_zscore = np.zeros(self.n_chans_original)
        amplitude_zscore[self.usable_idx] = (chan_amplitudes - amp_median) / amp_sd

        # Flag channels with amplitudes that deviate excessively from the median
        abnormal_amplitude = np.abs(amplitude_zscore) > deviation_threshold
        deviation_channel_mask = np.isnan(amplitude_zscore) | abnormal_amplitude

        # Update names of bad channels by excessive deviation & save additional info
        deviation_channels = self.ch_names_original[deviation_channel_mask]
        self.bad_by_deviation = deviation_channels.tolist()
        self._extra_info["bad_by_deviation"].update(
            {
                "median_channel_amplitude": amp_median,
                "channel_amplitude_sd": amp_sd,
                "robust_channel_deviations": amplitude_zscore,
            }
        )

def GVTD_EEG(data_np, gvtdTimeTrace, thresh, sfreq, data):
    # Global Variance in Temporal Derivative --------------------------------------------------
    # (GVTD) is a motion detection metric for multi-channel imaging
    # techniques inspired by the DVARS metric in fMRI. Introduced by Sherafati et al. 2020

    # dataMatrix = transpose(OD_MI_reconstructed)  # ._data  # from OD filtered data
    # np.savetxt('dataMatrix.txt', dataMatrix)

    bad_GVTD = []
    for ii in range(np.shape(data_np)[1]):
        if gvtdTimeTrace[ii] > thresh:
            a = data_np[0, ii]
            bad_GVTD.append(ii)
    # print('bad----------', len(bad_GVTD), bad_GVTD)
    # print('mean of matrix =', np.mean(np.mean(OD_MI_filter_np, axis=1)))
    # xs = [1, 0.1]
    # plt.plot(OD_MI_filter_np.T)
    # plt.vlines(x=[bad_GVTD], ymin=0, ymax=max(xs), color='black', alpha=0.8, lw=1.5)
    # plt.show()
    info = data.info
    raw_eeg_GVTD = mne.io.RawArray(data_np, info)

    # Automated annotation of bad segmants --------------------------------------------------
    bad_GVTD = np.array(bad_GVTD)
    # Calculate bad segment indices
    start = bad_GVTD
    bad_start_indices = (start - int(0.5 * sfreq)).astype(int)  # 0.5 second before
    bad_end_indices = (start + int(0.5 * sfreq)).astype(int)  # 0.5 second after
    bad_segment_indices = list(zip(bad_start_indices, bad_end_indices))
    # print('bad_segment_indices', bad_segment_indices)

    return raw_eeg_GVTD, info, bad_segment_indices

def eeg_reject(data, bad_segments):
    """
    Removes bad segments from EEG data.

    Parameters:
    - data: np.array (n_channels, n_samples), EEG data
    - bad_segments: list of tuples [(start_idx, end_idx)], indices of bad segments

    Returns:
    - cleaned_data: np.array (n_channels, new_n_samples), EEG data with bad segments removed
    """
    # Creiamo una maschera booleana per indicare i segmenti validi
    reject_mask = np.ones(data.shape[1], dtype=bool)

    # Segniamo i segmenti da rimuovere
    for segment in bad_segments:
        reject_mask[max(segment[0], 0):max(segment[1] + 1, 0)] = False

    # Applichiamo la maschera per rimuovere i segmenti cattivi
    cleaned_data = data[:, reject_mask]

    return cleaned_data


