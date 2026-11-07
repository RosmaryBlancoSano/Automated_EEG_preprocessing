"""Statistical analysis of the paper.

For every dataset, metric (and number of channels in the semi-simulated data) the methods, including the signal without denoising, are
compared with the Friedman test, with the simulated subjects, participants or recordings as blocks. The pairwise comparisons are paired
two-sided Wilcoxon signed-rank tests (exact up to 20 pairs, normal approximation otherwise), each with Bonferroni correction for the number
of comparisons within the same column (same metric, dataset and number of channels). Three families of comparisons:
  (i)   each method or setting against the signal without denoising;
  (ii)  the 10 pairs among the five default methods;
  (iii) the best value of the column against each of the others (Tables II and III), the ASR cutoffs and the artefact amplitudes.

Usage:  python statistical_analysis.py [simulated | real | amplitude | all]      (default: all)
Input:  the result files of 05_simulated_analysis (ALL_sim_rows.csv, sim_wt_coif4_variants.csv) and of 06_real_datasets
        (rec2_stroke/, rec2_dementia/, rec_auto/) in the results folder (environment variable EEGDN_RESULTS, or results/ inside or next to
        the repository).
Output: results/statistics/*.csv  (p-values after Bonferroni correction; p_adj >= 0.05 means no significant difference).
"""
import os, sys, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stat_tests import RD, friedman, wilcoxon, test_vs_base, load_dir

OUT = os.path.join(RD, 'statistics') + os.sep
os.makedirs(OUT, exist_ok=True)
ALPHA = 0.05

# ------------------------------------------------------------------------------------------------ semi-simulated data
def _load_sim():
    R = list(csv.DictReader(open(RD + 'ALL_sim_rows.csv')))
    V = list(csv.DictReader(open(RD + 'sim_wt_coif4_variants.csv')))
    return R, V

def _get(R, V, m, ch, met, scale=1.0):
    """Values of metric `met` of method `m` with `ch` channels, one per simulated subject (in the order of the subjects)."""
    rows = [r for r in R if r['Method'] == m and int(r['Nch']) == ch and float(r['Scale']) == scale]
    if not rows and scale == 1.0:   # the variants of the Coiflet 4 wavelet are in a separate file
        rows = [r for r in V if r['Method'] == m and int(r['Nch']) == ch]
    return [float(r[met]) for r in sorted(rows, key=lambda r: int(r['Subject']))]

DEFAULT = ['WT_repo', 'ASR_k5', 'EMD_2IMF', 'ASR5_WT', 'EMD2_WT']                                  # Table II
SETTINGS = ['ASR_k5', 'ASR_k10', 'ASR_k20', 'ASR_k30', 'WT_coif4_L5', 'WT_repo', 'WT_coif4_L9', 'WT_coif4_soft', 'WT_coif4_univ',
            'WT_den_hard', 'WT_den_soft', 'EMD_2IMF', 'EMD_3IMF', 'EMD_adapt', 'ASR5_WT', 'EMD2_WT']   # Table III
METRICS = {'SNR_aff': 1, 'CC_aff': 1, 'RRMSE_t_aff': -1, 'RRMSE_f_aff': -1}                      # +1: higher is better
NCH = (4, 8, 16, 31)

def simulated():
    R, V = _load_sim(); g = lambda m, ch, met: _get(R, V, m, ch, met)
    # (i) against the signal without denoising: Table II (four metrics) and Table III (SNR)
    rows = []
    for table, methods, metrics in (('II', DEFAULT, list(METRICS)), ('III', SETTINGS, ['SNR_aff'])):
        for ch in NCH:
            for met in metrics:
                fp, pp = test_vs_base(g('None', ch, met), {m: g(m, ch, met) for m in methods}, methods)
                rows += [[table, ch, met, m, fp, pp[m], len(methods)] for m in methods]
    _write('vs_no_denoising_simulated.csv', ['Table', 'Nch', 'Metric', 'Method', 'Friedman_p', 'p_adj', 'n_comparisons'], rows)
    # (ii) the 10 pairs among the five default methods
    rows = []; pairs = [(a, b) for i, a in enumerate(DEFAULT) for b in DEFAULT[i + 1:]]
    for ch in NCH:
        for met, sg in METRICS.items():
            D = {m: g(m, ch, met) for m in DEFAULT}; mean = {m: sum(v) / len(v) for m, v in D.items()}
            best = max(DEFAULT, key=lambda m: sg * mean[m])
            rows += [[ch, met, a, b, min(1, wilcoxon(D[a], D[b]) * len(pairs)), best] for a, b in pairs]
    _write('pairwise_default_methods_simulated.csv', ['Nch', 'Metric', 'Method_A', 'Method_B', 'p_adj', 'Best_method'], rows)
    # (iii) the best value of the column against each of the others: Tables II and III
    rows = []
    for table, methods, metrics in (('II', DEFAULT, list(METRICS)), ('III', SETTINGS, ['SNR_aff'])):
        k = len(methods) - 1
        for ch in NCH:
            for met in metrics:
                D = {m: g(m, ch, met) for m in methods}; mean = {m: sum(v) / len(v) for m, v in D.items()}
                best = max(methods, key=lambda m: METRICS[met] * mean[m])
                for m in methods:
                    p = float('nan') if m == best else min(1, wilcoxon(D[best], D[m]) * k)
                    rows.append([table, ch, met, best, m, mean[m], p, m == best or p >= ALPHA])
    _write('best_vs_others_simulated.csv', ['Table', 'Nch', 'Metric', 'Best_method', 'Method', 'Mean', 'p_adj', 'Not_different_from_best'], rows)
    # (iii) the ASR cutoffs: the best cutoff against the other three (3 comparisons) and against WT
    rows = []; ks = ['ASR_k5', 'ASR_k10', 'ASR_k20', 'ASR_k30']
    for ch in (31, 16, 8, 4):
        best = max(ks, key=lambda m: sum(g(m, ch, 'SNR_aff')))
        rows += [[ch, best, m, min(1, wilcoxon(g(best, ch, 'SNR_aff'), g(m, ch, 'SNR_aff')) * 3), 3] for m in ks if m != best]
        rows.append([ch, best, 'WT_repo', min(1, wilcoxon(g(best, ch, 'SNR_aff'), g('WT_repo', ch, 'SNR_aff'))), 1])
    _write('asr_cutoffs_simulated.csv', ['Nch', 'Best_cutoff', 'Compared_with', 'p_adj', 'n_comparisons'], rows)

# ------------------------------------------------------------------------------------------------ artefact amplitude
def amplitude():
    R, V = _load_sim(); M = ['None', 'WT_repo', 'ASR_k5', 'ASR_k20', 'EMD_2IMF', 'EMD_adapt', 'ASR5_WT', 'EMD2_WT']; rows = []
    for sc in (0.25, 0.5, 2, 4):
        for ch in (31, 16, 8, 4):
            D = {m: _get(R, V, m, ch, 'SNR_aff', sc) for m in M}
            if any(len(v) != 10 for v in D.values()): print('amplitude', sc, ch, 'missing', {m: len(v) for m, v in D.items()}); continue
            mean = {m: sum(v) / 10 for m, v in D.items()}; best = max(M, key=lambda m: mean[m]); k = len(M) - 1
            for m in M:
                vs_best = float('nan') if m == best else min(1, wilcoxon(D[best], D[m]) * k)
                vs_none = float('nan') if m == 'None' else min(1, wilcoxon(D[m], D['None']) * k)
                rows.append([sc, ch, m, mean[m], best, vs_best, vs_none])
    _write('artefact_amplitude_simulated.csv', ['Amplitude', 'Nch', 'Method', 'Mean_SNR', 'Best_method', 'p_adj_vs_best', 'p_adj_vs_no_denoising'], rows)

# ------------------------------------------------------------------------------------------------ real datasets
def real():
    DS = {'Stroke': load_dir('rec2_stroke'), 'Dementia': load_dir('rec2_dementia'), 'AutomaticArtefactRemoval': load_dir('rec_auto', 27)}
    default = ['WT', 'ASR_k5', 'EMD', 'EMD_WT', 'ASR_WT']
    cols = ['RMS_ratio', 'BPR_delta', 'BPR_theta', 'BPR_alpha', 'BPR_beta', 'Kurtosis', 'Entropy', 'CovPreserve', 'ACF_corr', 'GA50', 'GA100', 'GA150']
    sens = ['ASR_k5', 'ASR_k10', 'ASR_k20', 'ASR_k30', 'WT_L5', 'WT', 'WT_L9', 'WT_soft', 'EMD', 'EMD_3IMF', 'EMD_adapt']
    sens_cols = ['RMS_ratio', 'BPR_delta', 'GA100']; rows = []
    for ds, D in DS.items():
        recs = sorted(D['None'].keys())
        for table, methods, mcols in (('IV', default, cols), ('V', sens, sens_cols)):
            for c in mcols:
                base = [float(D['None'][r][c]) for r in recs]; oth = {m: [float(D[m][r][c]) for r in recs] for m in methods}
                fp, pp = test_vs_base(base, oth, methods)
                rows += [[table, ds, c, m, fp, pp[m], len(methods)] for m in methods]
    _write('vs_no_denoising_real.csv', ['Table', 'Dataset', 'Metric', 'Method', 'Friedman_p', 'p_adj', 'n_comparisons'], rows)

def _write(name, header, rows):
    with open(OUT + name, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
    print('written', OUT + name, len(rows), 'rows')

if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('simulated', 'all'): simulated()
    if what in ('amplitude', 'all'): amplitude()
    if what in ('real', 'all'): real()
