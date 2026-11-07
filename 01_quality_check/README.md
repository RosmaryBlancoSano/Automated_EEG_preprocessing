# 01_quality_check

Automatic quality check of the raw Stroke recordings (EDF, 8 channels).

```
python DataQualityCheck.py RAW_DIR [--subjects-csv FILE]
```

It filters (FIR notch at 25, 50 and 100 Hz, 1 Hz high-pass, 35 Hz low-pass), rejects high-amplitude segments with the Global Variance of the Temporal Derivative (`GVTD.py`, threshold = mean + 10 left-sided standard deviations, after Sherafati et al.), finds bad channels with the PREP deviation criterion and flat channels (`Functions.py`, from PyPREP), and classifies a recording as BAD if at least 25% of the channels are bad or at least 50% of the data were rejected. Bad channels are not interpolated here. Output: `GoodQuality/` and `BadQuality/` (read by `06_real_datasets/Dataset_Stroke_new.m`) and `EEG_quality_summary.csv`. The recordings cannot be shared.

| File | Content |
|---|---|
| `DataQualityCheck.py` | the quality check |
| `Functions.py` | `NoisyChannels` (adapted from PyPREP, MIT licence), `GVTD_EEG`, `eeg_reject` |
| `GVTD.py` | GVTD time trace and threshold |
| `ReadEDF32/` | EDF reader with 32-bit reading, adapted from MNE 0.22.1 (BSD-3) |
| `requirements.txt` | Python packages |

A first visual screening of the raw traces (recordings with technical acquisition failure, e.g. isoelectric signals or amplifier saturation, or with less than one minute of usable signal) was done manually before this script and is not part of the code.
