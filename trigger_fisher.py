import numpy as np, os

compiler = 'python3'
pgmname = 'fisher.py'
use_ilc_nl = 0
fwhm_arcmins = 1.4 #arcmins
rms_map_T_arr = np.arange(1., 10., 0.1)
round_results_arr = [1]#0, 1]

for round_results in round_results_arr:
    for rms_map_T in rms_map_T_arr:
        cmd = '%s %s -use_ilc_nl %s -rms_map_T %.3f -fwhm_arcmins %.2f -round_results %s' %(compiler, pgmname, use_ilc_nl, rms_map_T, fwhm_arcmins, round_results)
        print('\n#### %s ####' %(cmd))
        os.system(cmd)

