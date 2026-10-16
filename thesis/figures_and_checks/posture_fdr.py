"""Section 5.4: correct the session-trend p-values of the 40 foot measures for multiple testing.

Reads the summary written by analysis/posture_across_session.m (one row per measure:
4 trial types x 5 time windows x 2 size metrics) and prints Benjamini-Hochberg (false discovery
rate) and Bonferroni values for the measures reported in the thesis.
Usage: python posture_fdr.py [path to trend_summary_40features.csv]
"""
import csv, os, sys

SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.join('output', 'trial_evolution_40features', 'trend_summary_40features.csv')
rows = list(csv.DictReader(open(SRC)))
ps = {r['feature']: float(r['p_value']) for r in rows}
m = len(ps)
order = sorted(ps, key=ps.get)
q, prev = {}, 1.0
for i in range(m - 1, -1, -1):                       # Benjamini-Hochberg step-up
    prev = min(prev, ps[order[i]] * m / (i + 1))
    q[order[i]] = prev
print(f'{m} measures; {sum(p < 0.05 for p in ps.values())} with p < 0.05')
for f in ['Expected_Reflex(0to0.5s)_RMS', 'Unexpected_Reflex(0to0.5s)_RMS',
          'Expected_LateComp(1to1.5s)_RMS', 'Unexpected_LateComp(1to1.5s)_RMS']:
    print(f'{f:36s} p = {ps[f]:.4f}  FDR q = {q[f]:.3f}  Bonferroni = {min(1, ps[f] * m):.3f}')
