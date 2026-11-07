"""
Create object containing (meta)data from EDF file.

Adapted from the MNE package (0.22.1). 
Changes are marked with ARvNA/MNvS

Changes include:

    - Add 32 bits reading
    - Changed DC prefilter reading 
    - Add impedance reading
    - Remove logging
    - Remove reading annotations

Classes:

    RawEDF32

Functions:

    read_raw_edf32
    _get_info32
    _read_header32
    _read_edf_header32

"""

# pylint: skip-file

# built-in packages
from datetime import datetime, timezone
import math
import os
import re
import warnings

# third party packages
import mne.io.edf.edf as edf
from mne.io.base import BaseRaw
from mne.io.constants import FIFF
import numpy as np


class RawEDF32(edf.RawEDF):
    '''
    RawEDF32 is a subclass of RawEDF, which is a subclass of BaseRaw.
    It enables the reading of 32 bits .edf files.
    Also, it reads data stored in the reserved space of the header
    and stores it in self.info['impedance'] (channels);
                    self.info['ref_impedance'] (reference);
                    self.info['gnd_impedance'] (ground).
    '''

    def __init__(self, input_fname, eog=None, misc=None,
                stim_channel='auto', exclude=(), preload=False, verbose=None):

        input_fname = os.path.abspath(input_fname)
        info, edf_info, orig_units = _get_info32(input_fname,
                                                stim_channel, eog, misc,
                                                exclude, preload)

        # Raw attributes
        last_samps = [edf_info['nsamples'] - 1]
        BaseRaw.__init__(self, info, preload, filenames=[input_fname],
                            raw_extras=[edf_info], last_samps=last_samps,
                            orig_format='int', orig_units=orig_units,
                            verbose=verbose)

        print('RawEDF32 initialized...', info)
        # BaseRaw.plot(self)


def read_raw_edf32(input_fname, eog=None, misc=None,
                 stim_channel='auto', exclude=(), preload=False, verbose=None):
    """Create object containing info and data from EDF"""
    input_fname = os.path.abspath(input_fname)
    ext = os.path.splitext(input_fname)[1][1:].lower()
    if ext != 'edf':
        raise NotImplementedError(
            'Only EDF files are supported by read_raw_edf, got %s' % (ext,))
    return RawEDF32(input_fname=input_fname, eog=eog, misc=misc,
                  stim_channel=stim_channel, exclude=exclude, preload=preload,
                  verbose=verbose)


def _get_info32(fname, stim_channel, eog, misc, exclude, preload):
    """Extract all the information from the EDF+, BDF or GDF file."""
    eog = eog if eog is not None else []
    misc = misc if misc is not None else []

    edf_info, orig_units = _read_header32(fname, exclude)

    # XXX: `tal_ch_names` to pass to `edf._check_stim_channel` should be computed
    #      from `edf_info['ch_names']` and `edf_info['tal_idx']` but 'tal_idx'
    #      contains stim channels that are not TAL.
    stim_channel_idxs, _ = edf._check_stim_channel(
        stim_channel, edf_info['ch_names'])

    sel = edf_info['sel']  # selection of channels not excluded
    ch_names = edf_info['ch_names']  # of length len(sel)
    n_samps = edf_info['n_samps'][sel]
    nchan = edf_info['nchan']
    physical_ranges = edf_info['physical_max'] - edf_info['physical_min']
    cals = edf_info['digital_max'] - edf_info['digital_min']
    bad_idx = np.where((~np.isfinite(cals)) | (cals == 0))[0]
    if bad_idx.size > 0:
        warnings.warn('Scaling factor is not defined in following channels:\n' +
             ', '.join(ch_names[i] for i in bad_idx))
        cals[bad_idx] = 1
    bad_idx = np.where(physical_ranges == 0)[0]
    if len(bad_idx) > 0:
        warnings.warn('Physical range is not defined in following channels:\n' +
             ', '.join(ch_names[i] for i in bad_idx))
        physical_ranges[bad_idx] = 1

    # Creates a list of dicts of eeg channels for raw.info
    #logger.info('Setting channel info structure...')
    chs = list()
    pick_mask = np.ones(len(ch_names))

    for idx, ch_name in enumerate(ch_names):
        chan_info = {}
        chan_info['cal'] = 1.
        chan_info['logno'] = idx + 1
        chan_info['scanno'] = idx + 1
        chan_info['range'] = 1.
        chan_info['unit_mul'] = FIFF.FIFF_UNITM_NONE
        chan_info['ch_name'] = ch_name
        chan_info['unit'] = FIFF.FIFF_UNIT_V
        chan_info['coord_frame'] = FIFF.FIFFV_COORD_HEAD
        chan_info['coil_type'] = FIFF.FIFFV_COIL_EEG
        chan_info['kind'] = FIFF.FIFFV_EEG_CH
        chan_info['loc'] = np.zeros(12)
        if ch_name in eog or idx in eog or idx - nchan in eog:
            chan_info['coil_type'] = FIFF.FIFFV_COIL_NONE
            chan_info['kind'] = FIFF.FIFFV_EOG_CH
            pick_mask[idx] = False
        elif ch_name in misc or idx in misc or idx - nchan in misc:
            chan_info['coil_type'] = FIFF.FIFFV_COIL_NONE
            chan_info['kind'] = FIFF.FIFFV_MISC_CH
            pick_mask[idx] = False
        elif idx in stim_channel_idxs:
            chan_info['coil_type'] = FIFF.FIFFV_COIL_NONE
            chan_info['unit'] = FIFF.FIFF_UNIT_NONE
            chan_info['kind'] = FIFF.FIFFV_STIM_CH
            pick_mask[idx] = False
            chan_info['ch_name'] = ch_name
            ch_names[idx] = chan_info['ch_name']
            edf_info['units'][idx] = 1
        chs.append(chan_info)

    edf_info['stim_channel_idxs'] = stim_channel_idxs

    if any(pick_mask):
        picks = [item for item, mask in zip(range(nchan), pick_mask) if mask]
        edf_info['max_samp'] = max_samp = n_samps[picks].max()
    else:
        edf_info['max_samp'] = max_samp = n_samps.max()

    # Info structure
    # -------------------------------------------------------------------------

    not_stim_ch = [x for x in range(n_samps.shape[0])
                   if x not in stim_channel_idxs]
    sfreq = np.take(n_samps, not_stim_ch).max() * \
        edf_info['record_length'][1] / edf_info['record_length'][0]
    info = edf._empty_info(sfreq)
    info['meas_date'] = edf_info['meas_date']
    info['chs'] = chs
    info['ch_names'] = ch_names

    '''# Filter settings
    highpass = edf_info['highpass']
    lowpass = edf_info['lowpass']
    if highpass.size == 0.:
        pass
    elif all(highpass):
        if highpass[0] == 'NaN':
            pass  # Placeholder for future use. Highpass set in edf._empty_info.
        elif highpass[0] == 'DC':
            info['highpass'] = 0.
        else:
            hp = highpass[0]
            try:
                hp = float(hp)
            except Exception:
                hp = 0.
            info['highpass'] = hp
    else:
        numeric_highpass = []
        for hp_str in highpass:
            try:
                numeric_highpass.append(float(hp_str))
            except ValueError:
                numeric_highpass.append(0.0)  # or another default value

        if numeric_highpass:
            info['highpass'] = float(np.max(numeric_highpass))
        else:
            info['highpass'] = 0.0

        warnings.warn('Channels contain different highpass filters. Highest filter '
                      'setting will be stored.')

    if np.isnan(info['highpass']):
        info['highpass'] = 0.

    if lowpass.size == 0:
        pass  # Placeholder for future use. Lowpass set in edf._empty_info.
    elif all(lowpass):
        if lowpass[0] in ('NaN', '0.', '0.0'):
            pass  # Placeholder for future use. Lowpass set in edf._empty_info.
        else:
            info['lowpass'] = float(lowpass[0])
    else:
        info['lowpass'] = float(np.min(lowpass))
        warnings.warn('Channels contain different lowpass filters. Lowest filter '
             'setting will be stored.')
    if np.isnan(info['lowpass']):
        info['lowpass'] = info['sfreq'] / 2.

    if info['highpass'] > info['lowpass']:
        warnings.warn(f'Highpass cutoff frequency {info["highpass"]} is greater than '
             f'lowpass cutoff frequency {info["lowpass"]}, '
             'setting values to 0 and Nyquist.')
        info['highpass'] = 0.
        info['lowpass'] = info['sfreq'] / 2.
    '''

    def process_filter(filter_values, is_highpass=True):
        """Processes filter values (highpass or lowpass)."""

        if filter_values is None or not filter_values.size:  # Check for None and empty array
            return None

        numeric_filter = []
        for val in filter_values:
            try:
                numeric_val = float(val)
                if np.isfinite(numeric_val):
                    numeric_filter.append(numeric_val)
                else:
                    warnings.warn(f"Invalid filter value: {val}. Ignored.")
            except ValueError:
                if str(val).strip().lower() in ('nan', '0.', '0.0'):  # convert val in string
                    pass  # Ignore 'NaN', '0.' and '0.0'
                elif str(val).strip().lower() == 'dc' and is_highpass:  # convert val in string
                    numeric_filter.append(0.0)  # Set 'DC' to 0.0 for highpass
                else:
                    warnings.warn(f"Invalid filter value: {val}. Ignored.")

        if numeric_filter:
            return float(np.max(numeric_filter) if is_highpass else np.min(numeric_filter))
        else:
            return None


    info['highpass'] = process_filter(edf_info['highpass'], is_highpass=True)
    info['lowpass'] = process_filter(edf_info['lowpass'], is_highpass=False)

    # default values if highpass or lowpass are None
    info['highpass'] = info['highpass'] if info['highpass'] is not None else 0.0
    info['lowpass'] = info['lowpass'] if info['lowpass'] is not None else info['sfreq'] / 2.0

    # check whether highpass is greater than lowpass
    if info['highpass'] > info['lowpass']:
        warnings.warn(f'Highpass cutoff frequency {info["highpass"]} is greater than '
                      f'lowpass cutoff frequency {info["lowpass"]}, '
                      'setting values to 0 and Nyquist.')
        info['highpass'] = 0.0
        info['lowpass'] = info['sfreq'] / 2.0


    # DC prefiltering (MNvS)
    # info['dc_gain'] = edf_info['dc_gain']
    # info['dc_bandwidth'] = edf_info['dc_bandwidth']

    # Some keys to be consistent with FIF measurement info
    info['description'] = None
    edf_info['nsamples'] = int(edf_info['n_records'] * max_samp)

    info._update_redundant()

    # Later used for reading
    edf_info['cal'] = physical_ranges / cals

    # Physical dimension in µV
    edf_info['offsets'] = (
        edf_info['physical_min'] - edf_info['digital_min'] * edf_info['cal'])
    del edf_info['physical_min']
    del edf_info['digital_min']

    if edf_info['subtype'] == 'bdf':
        edf_info['cal'][stim_channel_idxs] = 1
        edf_info['offsets'][stim_channel_idxs] = 0
        edf_info['units'][stim_channel_idxs] = 1

    # Impedances (ARvNA)
    # info['impedance'] = edf_info['impedance']
    # info['ref_impedance'] = edf_info['ref_impedance']
    # info['gnd_impedance'] = edf_info['gnd_impedance']

    return info, edf_info, orig_units


def _read_header32(fname, exclude):
    """Unify edf, bdf and gdf _read_header call.
    Parameters
    ----------
    fname : str
        Path to the EDF+, BDF, or GDF file.
    exclude : list of str
        Channel names to exclude. This can help when reading data with
        different sampling rates to avoid unnecessary resampling.
    Returns
    -------
    (edf_info, orig_units) : tuple
    """
    ext = os.path.splitext(fname)[1][1:].lower()
    if ext in ('bdf', 'edf'):
        return _read_edf_header32(fname, exclude)
    else:
        raise NotImplementedError(
            f'Only GDF, EDF, and BDF files are supported, got {ext}.')


def _read_edf_header32(fname, exclude):
    """Read header information from EDF+ or BDF file."""
    edf_info = {'events': []}

    with open(fname, 'rb') as fid:

        # Bits amount # (MNvS)
        nbit = fid.read(8).decode('latin-1')
        edf_info['version'] = nbit

        # Patient ID
        patient = {}
        id_info = fid.read(80).decode('latin-1').rstrip()
        id_info = id_info.split(' ')
        if len(id_info):
            patient['id'] = id_info[0]
            if len(id_info) >= 2:
                try:
                    birthdate = datetime.strptime(id_info[2], "%d-%b-%Y")
                except ValueError:
                    birthdate = "X"
                patient['sex'] = id_info[1]
                patient['birthday'] = birthdate
                patient['name'] = id_info[3]

        # Recording ID
        meas_id = {}
        rec_info = fid.read(80).decode('latin-1').rstrip().split(' ')
        valid_startdate = False
        if len(rec_info) == 5:
            try:
                startdate = datetime.strptime(rec_info[1], "%d-%b-%Y")
            except ValueError:
                startdate = "X"
            else:
                valid_startdate = True
            meas_id['startdate'] = startdate
            meas_id['study_id'] = rec_info[2]
            meas_id['technician'] = rec_info[3]
            meas_id['equipment'] = rec_info[4]

        # If startdate available in recording info, use it instead of the
        # file's meas_date since it contains all 4 digits of the year
        if valid_startdate:
            day = meas_id['startdate'].day
            month = meas_id['startdate'].month
            year = meas_id['startdate'].year
            fid.read(8)  # skip file's meas_date
        else:
            meas_date = fid.read(8).decode('latin-1')
            day, month, year = [int(x) for x in meas_date.split('.')]
            year = year + 2000 if year < 85 else year + 1900

        meas_time = fid.read(8).decode('latin-1')
        hour, minute, sec = [int(x) for x in meas_time.split('.')]
        try:
            meas_date = datetime(year, month, day, hour, minute, sec,
                                 tzinfo=timezone.utc)
        except ValueError:
            warnings.warn(f'Invalid date encountered ({year:04d}-{month:02d}-'
                 f'{day:02d} {hour:02d}:{minute:02d}:{sec:02d}).')
            meas_date = None

        header_nbytes = int(edf._edf_str(fid.read(8)))

        # The following 44 bytes sometimes identify the file type, but this is
        # not guaranteed. Therefore, we skip this field and use the file
        # extension to determine the subtype (EDF or BDF, which differ in the
        # number of bytes they use for the data records; EDF uses 2 bytes
        # whereas BDF uses 3 bytes).
        fid.read(44)
        subtype = os.path.splitext(fname)[1][1:].lower()

        n_records = int(edf._edf_str(fid.read(8)))
        record_length = float(edf._edf_str(fid.read(8)))
        record_length = np.array([record_length, 1.])  # in seconds
        if record_length[0] == 0:
            record_length = record_length[0] = 1.
            warnings.warn('Header information is incorrect for record length. Default '
                 'record length set to 1.')

        nchan = int(edf._edf_str(fid.read(4)))
        channels = list(range(nchan))
        ch_names = [fid.read(16).strip().decode('latin-1') for ch in channels]
        exclude = edf._find_exclude_idx(ch_names, exclude)
        tal_idx = edf._find_tal_idx(ch_names)
        exclude = np.concatenate([exclude, tal_idx])
        sel = np.setdiff1d(np.arange(len(ch_names)), exclude)
        for _ in channels:
            fid.read(80)  # transducer
        units = [fid.read(8).strip().decode('latin-1') for ch in channels]
        edf_info['units'] = list()
        for i, unit in enumerate(units):
            if i in exclude:
                continue
            if unit == 'uV':
                edf_info['units'].append(1e-6)
            elif unit == 'mV':
                edf_info['units'].append(1e-3)
            else:
                edf_info['units'].append(1)
        edf_info['units'] = np.array(edf_info['units'], float)

        ch_names = [ch_names[idx] for idx in sel]
        units = [units[idx] for idx in sel]

        # make sure channel names are unique
        ch_names = edf._unique_channel_names(ch_names)
        orig_units = dict(zip(ch_names, units))

        physical_min = np.array(
            [float(edf._edf_str(fid.read(8))) for ch in channels])[sel]
        physical_max = np.array(
            [float(edf._edf_str(fid.read(8))) for ch in channels])[sel]

        if nbit == '32      ': # Check if bits amount is 32 (MNvS)
            digital_min = -np.array([float(int(edf._edf_str(fid.read(8)),16))
                                    for ch in channels])[sel]
            digital_max = np.array([float(int(edf._edf_str(fid.read(8)),16))
                                    for ch in channels])[sel]
        else:
            digital_min = np.array(
                [float(edf._edf_str(fid.read(8))) for ch in channels])[sel]
            digital_max = np.array(
                [float(edf._edf_str(fid.read(8))) for ch in channels])[sel]

        prefiltering = [edf._edf_str(fid.read(80)).strip() for ch in channels][:-1]
        highpass, lowpass = edf._parse_prefilter_string(prefiltering)

        # read gain and bandwidth of DC attenuation filter (MNvS)
        if prefiltering[0][0:5] == 'DCatt':
            dc_filter = re.findall(r'\d+\.\d+', prefiltering[0])
            dc_gain = dc_filter[0]
            dc_bandwidth = dc_filter[1]
        else:
            dc_gain = 'NaN'
            dc_bandwidth = 'NaN'

        # number of samples per record
        n_samps = np.array([int(edf._edf_str(fid.read(8))) for ch in channels])

        # read impedances from reserved field (ARvNA)
        impedance = [fid.read(32).strip().decode('latin-1')[4:-2] for ch in channels[:-1]]
        try:
            impedance = [int(item) for item in impedance]
        except:
            impedance = [math.nan for ch in channels[:-1]]
        edf_info['impedance'] = impedance
        ref_gnd_imp = fid.read(32).decode('latin-1')[4:].strip().split('|')
        try:
            ref_impedance = int(ref_gnd_imp[0][5:-2])
        except:
            ref_impedance = math.nan
        edf_info['ref_impedance'] = ref_impedance
        try:
            gnd_impedance = int(ref_gnd_imp[1][5:-2])
        except:
            gnd_impedance = math.nan
        edf_info['gnd_impedance'] = gnd_impedance

        # Populate edf_info
        edf_info.update(
            ch_names=ch_names, data_offset=header_nbytes,
            digital_max=digital_max, digital_min=digital_min,
            highpass=highpass, dc_gain=dc_gain, dc_bandwidth=dc_bandwidth, # (MNvS)
            sel=sel, lowpass=lowpass, meas_date=meas_date,
            n_records=n_records, n_samps=n_samps, nchan=nchan,
            subject_info=patient, physical_max=physical_max,
            physical_min=physical_min, record_length=record_length,
            subtype=subtype, tal_idx=tal_idx,
            impedance=impedance, ref_impedance=ref_impedance, # (ARvNA)
            gnd_impedance=gnd_impedance,
            blob=None)                      # (ARvNA)

        assert fid.tell() == header_nbytes

        fid.seek(0, 2)
        n_bytes = fid.tell()
        n_data_bytes = n_bytes - header_nbytes

        if subtype == 'bdf':
            total_samps = n_data_bytes // 3
        elif nbit == '32      ':
            total_samps = n_data_bytes // 4  # (MNvS)
        else:
            total_samps = n_data_bytes // 2

        read_records = total_samps // np.sum(n_samps)
        if n_records != read_records:
            warnings.warn('Number of records from the header does not match the file '
                 'size (perhaps the recording was not stopped before exiting).'
                 ' Inferring from the file size.')
            edf_info['n_records'] = read_records
        del n_records

        if subtype == 'bdf':
            edf_info['dtype_byte'] = 3  # 24-bit (3 byte) integers
            edf_info['dtype_np'] = UINT8
        elif nbit == '32      ':
            edf_info['dtype_byte'] = 4  # 32-bit (4 byte) integers (MNvS)
            edf_info['dtype_np'] = INT32
        else:
            edf_info['dtype_byte'] = 2  # 16-bit (2 byte) integers
            edf_info['dtype_np'] = INT16

    return edf_info, orig_units




INT8 = '<i1'
UINT8 = '<u1'
INT16 = '<i2'
UINT16 = '<u2'
INT32 = '<i4'
UINT32 = '<u4'
INT64 = '<i8'
UINT64 = '<u8'
FLOAT32 = '<f4'
FLOAT64 = '<f8'
GDFTYPE_NP = (None, INT8, UINT8, INT16, UINT16, INT32, UINT32,
              INT64, UINT64, None, None, None, None,
              None, None, None, FLOAT32, FLOAT64)
GDFTYPE_BYTE = tuple(np.dtype(x).itemsize if x is not None else 0
                     for x in GDFTYPE_NP)


if __name__ == "__main__":
    edf_root_dir = os.environ.get('EEGDN_EDF_EXAMPLE_DIR', '.')
    edf_file = os.path.join(edf_root_dir, 'sub-001_task-rest_acq-1_eeg.edf')
    myInstance = RawEDF32(edf_file, GDFTYPE_NP, GDFTYPE_BYTE)

