"""
#CMB fisher forecasting
"""
############################################################################################################

import numpy as np, scipy as sc, sys, argparse, os
sys.path.append('modules')
import tools


############################################################################################################
#get the necessary arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('-paramfile', dest='paramfile', action='store', help='paramfile', type=str, default='params/params_planck_r_0.0_2015_cosmo_lensed_LSS_JM.txt')
parser.add_argument('-which_spectra', dest='which_spectra', action='store', help='which_spectra', type=str, default='lensed_scalar', choices=['lensed_scalar', 'unlensed_scalar'])
##parser.add_argument('-use_thetastar', dest='use_thetastar', action='store', type=int, help='use_thetastar', default= 1)
##parser.add_argument('-use_cosmomc_theta', dest='use_cosmomc_theta', action='store', type=int, help='use_cosmomc_theta', default= 0)

args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
    param_value = args_keys[kargs]
    if isinstance(param_value, str):
        cmd = '%s = "%s"' %(kargs, param_value)
    else:
        cmd = '%s = %s' %(kargs, param_value)
    exec(cmd)

############################################################################################################
#folders containing inputs
camb_folder = 'data/CMB_spectra_derivatives_for_code_comparison/'
draft_results_folder = 'data/DRAFT_results_20200601/s4like_mask/TT-EE-TE/baseline/'
############################################################################################################

#get fiducual LCDM power spectra computed using CAMB
print('\n\tget/read fiducual LCDM power spectra computed using CAMB')
camb_fname = '%s/cmb_spectra_%s_Srini.txt' %(camb_folder, which_spectra.replace('_scalar',''))

cl_camb = np.loadtxt(camb_fname)
els = cl_camb[:,0]
cl_dic = {}
cl_dic['TT'] = cl_camb[:,1]
cl_dic['EE'] = cl_camb[:,2]
cl_dic['TE'] = cl_camb[:,3]

############################################################################################################
#read derivatives
print('\n\tget/read derivatives')
camb_deriv_fname = '%s/cmb_spectra_derivs_%s_Srini.npy' %(camb_folder, which_spectra.replace('_scalar',''))
cl_deriv_dic_tmp = np.load(camb_deriv_fname, allow_pickle = 1).item()
cl_deriv_dic = {}
param_names = []
for p in sorted( cl_deriv_dic_tmp ):
    if p == 'ell': continue
    cl_deriv_dic[p]={}
    cl_deriv_dic[p]['TT'] = cl_deriv_dic_tmp[p][0]
    cl_deriv_dic[p]['EE'] = cl_deriv_dic_tmp[p][1]
    cl_deriv_dic[p]['TE'] = cl_deriv_dic_tmp[p][2]
    param_names.append( p )
############################################################################################################
#get experiment specs
print('\n\tget experiment specs and nl')
pspectra_to_use = ['TT', 'EE', 'TE']
min_l_temp, max_l_temp = 30, 5000
min_l_pol, max_l_pol = 30, 5000
fix_params = [] #['mnu', 'ws'] but curretnly nothing to fix as we only have a 6+1(neff) LCDM model
prior_dic = {'tau':0.007}#02}#02} #Planck tau prior
desired_param = 'neff' #desired parameter for which we are computing the constraints. set to None if you want to analyse the full fisher matrix
include_gal = 0
gal_mask = 3 #only valid if galaxy is included
fsky = 0.77
############################################################################################################
#get nl
print('\n\tget nl')
if not include_gal:
    nlfile = '%s/S4_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_AZ.npy' %(draft_results_folder)
else:
    nlfile = '%s/S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_galmask%s_AZ.npy' %(draft_results_folder, gal_mask)
nl_dic, fsky = tools.get_nldic(nlfile, els)
############################################################################################################
#get delta_cl
print('\n\tget delta Cl')
delta_cl_dic = tools.fn_delta_Cl(els, cl_dic, nl_dic, fsky)
############################################################################################################
#get Fisher / COV matrices
print('\n\tget fisher')
F_mat = tools.get_fisher_mat(els, cl_deriv_dic, delta_cl_dic, param_names, pspectra_to_use = pspectra_to_use,\
            min_l_temp = min_l_temp, max_l_temp = max_l_temp, min_l_pol = min_l_pol, max_l_pol = max_l_pol)
#print(F_mat)
############################################################################################################
#fix params
print('\n\tfixing paramaters, if need be')
F_mat, param_names = tools.fn_fix_params(F_mat, param_names, fix_params)
param_names = np.asarray(param_names)
#print(F_mat)
############################################################################################################
#add prior
print('\n\tadding prior')
F_mat = tools.fn_add_prior(F_mat, param_names, prior_dic)
############################################################################################################
#get cov matrix now
print('\n\tget covariance matrix')
Cov_mat = sc.linalg.pinv2(F_mat) #made sure that COV_mat_l * Cinv_l ~= I
############################################################################################################
#extract sigma(neff)
if desired_param is not None:
    print('\n\textract sigma(neff)')
    pind = np.where(param_names == desired_param)[0][0]
    pcntr1, pcntr2 = pind, pind
    cov_inds_to_extract = [(pcntr1, pcntr1), (pcntr1, pcntr2), (pcntr2, pcntr1), (pcntr2, pcntr2)]
    cov_extract = np.asarray( [Cov_mat[ii] for ii in cov_inds_to_extract] ).reshape((2,2))
    sigma = cov_extract[0,0]**0.5
    sigma = round(sigma, 4)
    opline = '\n\t\t\simga(neff) = %.4f using observables = %s; fsky = %s; power spectra = %s\n' %(sigma, str(pspectra_to_use), fsky, which_spectra)
    print(opline)
sys.exit()
############################################################################################################







