"""
#CMB fisher forecasting
"""
############################################################################################################

import numpy as np, scipy as sc, sys, argparse, os
sys.path.append('modules')
import exp_specs, tools


############################################################################################################
#get the necessary arguments
parser = argparse.ArgumentParser(description='')
parser.add_argument('-paramfile', dest='paramfile', action='store', help='paramfile', type=str, default='params/params_planck_r_0.0_2015_cosmo_lensed_LSS_JM.txt')
parser.add_argument('-delta_l', dest='delta_l', action='store', help='delta_l', type=int, default=1)
parser.add_argument('-which_spectra', dest='which_spectra', action='store', help='which_spectra', type=str, default='lensed_scalar')
parser.add_argument('-use_thetastar', dest='use_thetastar', action='store', type=int, help='use_thetastar', default= 1)
parser.add_argument('-use_cosmomc_theta', dest='use_cosmomc_theta', action='store', type=int, help='use_cosmomc_theta', default= 0)
parser.add_argument('-Alens', dest='Alens', action='store', help='Alens', type=float, default=None)
parser.add_argument('-add_lensing', dest='add_lensing', action='store', help='add_lensing', type=int, default=0)
parser.add_argument('-also_te', dest='also_te', action='store', help='also_te', type=int, default=1)
parser.add_argument('-expname', dest='expname', action='store', help='expname', type=str, default = 's4')
parser.add_argument('-include_planck_on_large_scales', dest='include_planck_on_large_scales', action='store', help='include_planck_on_large_scales', type=int, default=0)##0)
parser.add_argument('-set_YHe_using_BBN', dest='set_YHe_using_BBN', action='store', type=int, default= 1, help='set_YHe_using_BBN')

## only for S4
parser.add_argument('-totlf_s4wide', dest='totlf_s4wide', action='store', help='totlf_s4wide', type=float, default=2)
parser.add_argument('-totmf_s4wide', dest='totmf_s4wide', action='store', help='totmf_s4wide', type=float, default=12)
parser.add_argument('-tothf_s4wide', dest='tothf_s4wide', action='store', help='tothf_s4wide', type=float, default=5)
parser.add_argument('-include_gal', dest='include_gal', action='store', type=int, default= 1, help='include_gal')
parser.add_argument('-galmask', dest='galmask', action='store', type=int, default= -1, help='galmask')
parser.add_argument('-force_fsky_val', dest='force_fsky_val', action='store', type=float, default= -1., help='force_fsky_val')

#lensing stuffs
parser.add_argument('-just_Aphiphi', dest='just_Aphiphi', action='store', help='just_Aphiphi', type=int, default = 0)
parser.add_argument('-lensing_est', dest='lensing_est', action='store', help='lensing_est', type=str, default = 'MV')
parser.add_argument('-include_atm_noise', dest='include_atm_noise', action='store', help='include_atm_noise', type=int, default=1)##0)
parser.add_argument('-use_ilc_N0', dest='use_ilc_N0', action='store', help='use_ilc_N0', type=int, default=1)##0)
parser.add_argument('-mod_lensing_N0_fac', dest='mod_lensing_N0_fac', action='store', type=float, default= 1., help='mod_lensing_N0_fac')

#code checks
#parser.add_argument('-old_stepsizes', dest='old_stepsizes', action='store', type=int, default= 0, help='old_stepsizes')
#parser.add_argument('-new_h_stepsize', dest='new_h_stepsize', action='store', type=int, default= 1, help='new_h_stepsize')
#parser.add_argument('-old_code', dest='old_code', action='store', type=float, default= 0, help='old_code')


args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
    param_value = args_keys[kargs]
    if isinstance(param_value, str):
        cmd = '%s = "%s"' %(kargs, param_value)
    else:
        cmd = '%s = %s' %(kargs, param_value)
    exec(cmd)
delta_l = int(delta_l)
if delta_l != 1: 
    print('\n\n\t\t delta_l must be 1 currently. aborting code\n\n')
    sys.exit()
############################################################################################################
#initialise and set cosmology
print('\n\n\n\tinitialise and set cosmology')
param_dict = tools.fn_ini_param_dic(paramfile)
if set_YHe_using_BBN:
    param_dict['YHe'] = None
else:
    param_dict['YHe'] = 0.2454006
#print(param_dict)
#cosmo = astropy.cosmology.FlatLambdaCDM(H0 = param_dict['h']*100., Om0 = param_dict['omega_m'])

min_l_limit, max_l_limit = param_dict['min_l_limit'], param_dict['max_l_limit']

if Alens is not None:
    param_dict['Alens'] = Alens
Alens = param_dict['Alens']

#get fiducual LCDM power spectra
print('\n\tget fiducual LCDM power spectra')
pars, els, Cl_dic = tools.fn_set_CAMB_como(param_dict, which_spectra, add_lensing = add_lensing, use_thetastar = use_thetastar, use_cosmomc_theta = use_cosmomc_theta)
if which_spectra == 'lens_potential' and just_Aphiphi:# add_lensing and which_spectra != 'lens_potential':
    pars, els, Cl_lens_dic = tools.fn_set_CAMB_como(param_dict, 'lens_potential')
    for k in Cl_lens_dic:
        Cl_dic[k] = Cl_lens_dic[k]

if expname.find('s4like_simple')>-1:
    #store spectra derivatives
    dropbox_folder = '/Users/sraghunathan/Dropbox/S4_Neff/results/'

    #save CMB spectra
    cmb_spectra_arr = np.asarray( [els, Cl_dic['TT'], Cl_dic['EE'], Cl_dic['TE'], Cl_dic['BB']] ).T

    #dropbox_folder_updated = '%s/CMB_spectra/' %(dropbox_folder)
    #dropbox_folder_updated = '%s/CMB_spectra_v2/' %(dropbox_folder)
    #dropbox_folder_updated = '%s/CMB_spectra_v3/' %(dropbox_folder)
    #dropbox_folder_updated = '%s/CMB_spectra_v4/' %(dropbox_folder)
    #dropbox_folder_updated = '%s/CMB_spectra_v5/' %(dropbox_folder)

    dropbox_folder_updated = '/Users/sraghunathan/Dropbox/S4_Neff/data/CMB_spectra_derivatives_for_code_comparison/'

    os.system('mkdir -p %s' %(dropbox_folder_updated))

    cmb_spectra_opf = '%s/cmb_spectra_%s_Srini.txt' %(dropbox_folder_updated, which_spectra.split('_')[0])
    np.savetxt(cmb_spectra_opf, cmb_spectra_arr, header = 'ell tt ee te bb')
    #sys.exit()

############################################################################################################
#get experiment specs
print('\n\tget experiment specs')
Nlfile = None
if expname == 's4':
    if galmask == -1:
        galmask = 1
elif expname.find('spt')>-1:
    if which_spectra == 'lens_potential':
        
        if not use_ilc_N0:
            if include_atm_noise:
                sptlensingnoisecurvesfolder = 'data/spt_summer_fields/lensing_noise_curves/with_oneoverfnoise/'
            else:
                sptlensingnoisecurvesfolder = 'data/spt_summer_fields/lensing_noise_curves/white_noise/'
        else:
            sptlensingnoisecurvesfolder = 'data/spt_summer_fields/lensing_noise_curves/ilc/'
        lmin, lmax = 100, 3000
        Nlfile = '%s/%s_lmin%s_lmax%s.npy' %(sptlensingnoisecurvesfolder, expname, lmin, lmax)
        if (0):
            Nlfile = '%s/%s_lmin400_lmax3000.npy' %(sptlensingnoisecurvesfolder, expname)
            #Nlfile = '%s/%s_lmin400_lmax3000_withplanck.npy' %(sptlensingnoisecurvesfolder, expname)
        print('\n\tNlfile is: %s\n' %(Nlfile))

if include_planck_on_large_scales:
    exp_dic_planck, expnamearr_planck = exp_specs.define_exp_specs('planck', which_spectra = which_spectra, els = els, \
        include_gal = include_gal, galmask = galmask, totlf = totlf_s4wide, totmf = totmf_s4wide, tothf = tothf_s4wide, also_te = 1, Nlfile = Nlfile,\
        lensing_est = lensing_est, mod_lensing_N0_fac = mod_lensing_N0_fac, add_lensing = add_lensing)
    Nldic_planck = exp_dic_planck[expnamearr_planck[0]]['Nldic']
else:
    Nldic_planck = None

exp_dic, expnamearr = exp_specs.define_exp_specs(expname, which_spectra = which_spectra, els = els, \
    include_gal = include_gal, galmask = galmask, totlf = totlf_s4wide, totmf = totmf_s4wide, tothf = tothf_s4wide, also_te = 1, Nlfile = Nlfile,\
    lensing_est = lensing_est, Nldic_for_inv_var = Nldic_planck, include_planck_on_large_scales = include_planck_on_large_scales, force_fsky_val = force_fsky_val, mod_lensing_N0_fac = mod_lensing_N0_fac, add_lensing = add_lensing)

totexp = len(expnamearr)
############################################################################################################

print('\n\tget delta Cl')
exp_delta_Cl_dic = {}
for expcnt, exp in enumerate(expnamearr):
    #if expcnt%500 == 0:
    #    print('\t\t%s of %s' %(expcnt+1, len(expnamearr)))

    #delta_Cl = fn_delta_Cl(els, Cl_TT, exp_dic[exp])

    binned_Cl_dic = {}
    for spec in Cl_dic:
        Cl = Cl_dic[spec]        
        binned_el, binned_Cl = tools.fn_el_binning(els, Cl, delta_l = delta_l)
        binned_Cl_dic[spec] = binned_Cl

    delta_Cl_dic = tools.fn_delta_Cl(binned_el, binned_Cl_dic, exp_dic[exp], Nlfile = exp_dic[exp]['Nlfile'], Nldic = exp_dic[exp]['Nldic'], delta_l = delta_l)
    exp_delta_Cl_dic[exp] = delta_Cl_dic

els = np.copy(binned_el)
Cl_dic = binned_Cl_dic
#################################################################################
#################################################################################

if (1):##not old_stepsizes: #based on Joel Meyers' step sizes
    cosmo_param_dict = {\
    'ombh2' : (0.0008,(0.0005, 0.000001)),\
    'omch2' : (0.0030,(0.005, 0.00001)),
    'tau' : (0.020, (0.1, 0.0001)),\
    'As' : (0.1e-9,(0.5,1e-3)),\
    'ns' : (0.010, (0.05, 0.0001)),\
    'ws' : (-1e-2,(0.5, 0.0001)),\
    'neff': (0.080,(0.5,0.01)),\
    'mnu': (0.02,(1.,0.0001)),\
    ##############'YHe': (0.005,(1.,0.0001)),\
    ###'Alens': (1e-2, (1., 0.3)),\
    ###'Aphiphi': (1e-2, (1., 0.3)),\
    }
    if use_thetastar:
        cosmo_param_dict['thetastar'] = (0.000050, (0.05,0.001))
        #cosmo_param_dict['thetastar'] = (0.000010, (0.05,0.001))
    elif use_cosmomc_theta:
        cosmo_param_dict['cosmomc_theta'] = (0.000050, (0.05,0.001))
    else:
        #cosmo_param_dict['h']= (.00001, (0.05,0.001))
        cosmo_param_dict['h']= (.005, (0.05,0.001))
    if not set_YHe_using_BBN:
        cosmo_param_dict['YHe']= (0.005,(1.,0.0001))


fix_params = ['ws', 'neff', 'mnu', 'Alens', 'Aphiphi']
#fix_params = ['As']
#fix_params = ['ws', 'mnu']

cosmo_param_pl_chars_dict = {\
'ombh2' : [1, r'$\Omega_{b}h^{2}$'], 
'omch2' : [2, r'$\Omega_{c}h^{2}$'], 
'tau' : [4, r'$\tau$'], 
'As' : [5, r'$A_{s}\ [10^{9}]$'], 
'ns' : [6, r'$n_{s}$'],
'ws' : [7, r'$w_{0}$'], 
'neff': [8, r'$N_{\rm eff}$'], 
'mnu': [9, r'$\sum m_{\nu}$'],
'YHe': [10, r'$Y_{p}$'],
'Alens': (11, r'$A_{\rm lens}$'),\
'Aphiphi': (11, r'$A_{\phi \phi}$'),\
}
if use_thetastar:
    cosmo_param_pl_chars_dict['thetastar'] = [3, r'$theta_{s}$']
elif use_cosmomc_theta:
    cosmo_param_dict['cosmomc_theta'] = (0.000050, (0.05,0.001))
else:
    cosmo_param_pl_chars_dict['h'] = [3, r'$h$']

if expname.find('spt')>-1: #SPT CMB+lensing
    cosmo_param_dict = {\
    'ombh2' : (0.0008,(0.0005, 0.000001)),\
    'omch2' : (0.0030,(0.005, 0.00001)),
    'tau' : (0.020, (0.1, 0.0001)),\
    'As' : (0.1e-9,(0.5,1e-3)),\
    'ns' : (0.010, (0.05, 0.0001)),\
    'ws' : (-1e-2,(0.5, 0.0001)),\
    'neff': (0.080,(0.5,0.01)),\
    'mnu': (0.02,(1.,0.0001)),\
    'Alens': (1e-2, (1., 0.3)),\
    'Aphiphi': (1e-2, (1., 0.3)),\
    'omk': (0.001,(0.,0.0001)),\
    ####'YHe': (0.005,(1.,0.0001)),\
    }

    if use_thetastar:
        cosmo_param_dict['thetastar'] = (0.000050, (0.05,0.001))
        #cosmo_param_dict['thetastar'] = (0.000010, (0.05,0.001))
    elif use_cosmomc_theta:
        cosmo_param_dict['cosmomc_theta'] = (0.000050, (0.05,0.001))
    else:
        #cosmo_param_dict['h']= (.00001, (0.05,0.001))
        cosmo_param_dict['h']= (.005, (0.05,0.001))

if which_spectra == 'lens_potential' and just_Aphiphi:
    cosmo_param_dict = {\
    'Aphiphi': (1e-3, (1., 0.3)),\
    }
    fix_params = []
    cosmo_param_pl_chars_dict = {\
    'Aphiphi': (11, r'$A_{\phi \phi}$'),\
    }

if expname == 's4like_simple':
    cosmo_param_dict = {\
    'ombh2' : (0.0008,(0.0005, 0.000001)),\
    'omch2' : (0.0030,(0.005, 0.00001)),
    'tau' : (0.020, (0.1, 0.0001)),\
    'As' : (0.1e-9,(0.5,1e-3)),\
    'ns' : (0.010, (0.05, 0.0001)),\
    ###'ws' : (-1e-2,(0.5, 0.0001)),\
    'neff': (0.080,(0.5,0.01)),\
    'mnu': (0.02,(1.,0.0001)),\
    ##############'YHe': (0.005,(1.,0.0001)),\
    ###'Alens': (1e-2, (1., 0.3)),\
    ###'Aphiphi': (1e-2, (1., 0.3)),\
    }
    if use_thetastar:
        cosmo_param_dict['thetastar'] = (0.000050, (0.05,0.001))
        #cosmo_param_dict['thetastar'] = (0.000010, (0.05,0.001))
    elif use_cosmomc_theta:
        cosmo_param_dict['cosmomc_theta'] = (0.000050, (0.05,0.001))
    else:
        #cosmo_param_dict['h']= (.00001, (0.05,0.001))
        cosmo_param_dict['h']= (.005, (0.05,0.001))
    if not set_YHe_using_BBN:
        cosmo_param_dict['YHe']= (0.005,(1.,0.0001))

param_names = sorted(cosmo_param_dict.keys())
#################################
#get derivatives
print('\n\tget derivatives')

'''
derived_param_names_dic_to_store = None
if (0):#expname.find('spt')>-1:
    #derived_param = ['sigma8', 'results.get_sigmaR(R=8., hubble_units=False, return_R_z=False)[0]']
    derived_param_names_dic = {\
        'As_e-2tau': '(1e9 * pars.InitPower.As) * np.exp( -2. * pars.Reion.optical_depth )',\
        'omegam': 'pars.omegam',\
        #'omegab': 'pars.omegab',\
        #'omegac': 'pars.omegac',\
        'sigma8': 'results.get_sigma8_0()',\
        #'sigma8': 'results.get_sigmaR(R=8., hubble_units=False, return_R_z=False)[0]', \
        }
    derived_param = []
    for der_p in derived_param_names_dic:
        derived_param.append([der_p, derived_param_names_dic[der_p]])
    Cl_deriv_dic, derived_param_deriv = tools.fn_get_dervatives(Cl_dic, param_dict, cosmo_param_dict, which_spectra, delta_l = delta_l, add_lensing = add_lensing, use_thetastar = use_thetastar, \
    use_cosmomc_theta = use_cosmomc_theta, \
    derived_param = derived_param, both_Cl_param = 1)#, fix_params = fix_params)

    derived_param_names_dic_to_store = {}
    derived_param_names_dic_to_store['derived_params'] = derived_param_names_dic
    derived_param_names_dic_to_store['derivatives'] = derived_param_deriv
'''


if which_spectra == 'lens_potential' and just_Aphiphi:
    Aphiphi, delta_Aphiphi = param_dict['Aphiphi'], cosmo_param_dict['Aphiphi'][0]
    Cl_deriv_dic = tools.fn_get_dervatives_Aphiphi(Cl_dic, Aphiphi, delta_Aphiphi)
else:
    Cl_deriv_dic = tools.fn_get_dervatives(Cl_dic, param_dict, cosmo_param_dict, which_spectra, delta_l = delta_l, add_lensing = add_lensing, use_thetastar = use_thetastar, use_cosmomc_theta = use_cosmomc_theta)#, fix_params = fix_params)

if expname.find('s4like_simple')>-1:

    #store spectra derivatives
    #save CMB spectra derivatives
    Cl_deriv_dic_to_store = {}
    for tmppname in Cl_deriv_dic:
        Cl_deriv_dic_to_store[tmppname] = np.asarray( [Cl_deriv_dic[tmppname]['TT'], Cl_deriv_dic[tmppname]['EE'], Cl_deriv_dic[tmppname]['TE']] )
        print(tmppname)
    Cl_deriv_dic_to_store['ell'] = els

    #dropbox_folder_updated = '%s/CMB_spectra_derivs/' %(dropbox_folder)
    #dropbox_folder_updated = '%s/CMB_spectra_derivs_v2/' %(dropbox_folder)

    dropbox_folder_updated = '/Users/sraghunathan/Dropbox/S4_Neff/data/CMB_spectra_derivatives_for_code_comparison/'

    os.system('mkdir -p %s' %(dropbox_folder_updated))
    
    cmb_spectra_deriv_opf = '%s/cmb_spectra_derivs_%s_Srini.npy' %(dropbox_folder_updated, which_spectra.split('_')[0])
    np.save(cmb_spectra_deriv_opf, Cl_deriv_dic_to_store)
    sys.exit()


###from IPython import embed; embed()
#################################
#get Fisher / COV matrices
print('\n\tget fisher')
F_dic = {}
Cov_dic = {}

#for exp in sorted( exp_dic ):
for cntr, keyname in enumerate( expnamearr ):

    ##print keyname
    #print('\t\t%s' %(keyname))
    #expname_split = expname.split('_')

    if expname == 's4' or expname == 'planck':
        param_names_str = '-'.join(param_names)
        opfolder = 'results/20200601/'
        #opfolder = 'results/S4_ilc_20203030/'
        #opfolder = 'results/S4_ilc_20205xx_with_gal/'

        '''
        if old_stepsizes:
            opfolder = '%s/old_stepsizes' %(opfolder)
        else:
            opfolder = '%s/new_stepsizes' %(opfolder)

        try:
            opfolder = '%s/new_h_stepsize/' %(opfolder)
        except:
            pass
        '''

        if set_YHe_using_BBN:
            opfolder = '%s/Yp_using_BBN/' %(opfolder)

        if use_thetastar:
            opfolder = '%s/with_thetastar' %(opfolder)
        elif use_cosmomc_theta:
            opfolder = '%s/with_cosmomc_theta' %(opfolder)

        if param_dict['Alens'] != 1.0:
            opfolder = '%s/%s_Alens%.2f/s4like_mask/TT-EE-TE/' %(opfolder, which_spectra, param_dict['Alens'])
        else:
            opfolder = '%s/%s/s4like_mask/TT-EE-TE/' %(opfolder, which_spectra)

        if totlf_s4wide == 2 and totmf_s4wide == 12 and tothf_s4wide == 5:
            opfolder = '%s/baseline/' %(opfolder)
        else:
            opfolder = '%s/tubes_mod/' %(opfolder)

        if include_planck_on_large_scales:
            opfolder = '%s/with_planck_on_large_scales' %(opfolder)

        if expname == 'planck':
            opfolder = '%s/%s' %(opfolder, expname)

        if add_lensing:
            opfolder = '%s/with_lensing' %(opfolder)

        if not os.path.exists(opfolder): os.system('mkdir -p %s' %(opfolder))
        opfname = '%s/%s_%s-%scosmo.npy' %(opfolder, keyname, param_names_str,len(param_names))

        print(opfname)
        ##remove me later
        ###if os.path.exists(opfname): continue
        ##remove me later

        ##from IPython import embed;embed()

    else:##if expname.find('spt')>-1:

        param_names_str = '-'.join(param_names)
        if expname.find('spt')>-1:
            opfolder = 'results/spt/'
        else:
            opfolder = 'results/%s/' %(expname)

        opfolder = '%s/%s' %(opfolder, which_spectra)

        if include_planck_on_large_scales:
            opfolder = '%s/with_planck_on_large_scales' %(opfolder)

        if add_lensing:
            opfolder = '%s/with_lensing' %(opfolder)

        if not os.path.exists(opfolder): os.system('mkdir -p %s' %(opfolder))

        opfname = '%s/%s_%s-%scosmo.npy' %(opfolder, keyname, param_names_str,len(param_names))
        #if os.path.exists(opfname): continue

        ##from IPython import embed;embed()

    print('\t\t%s of %s: %s' %(cntr+1, totexp, keyname))

    if (which_spectra == 'lens_potential' and just_Aphiphi):

        F_mat = tools.fn_fisher_forecast_Aphiphi(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'], \
            min_l = exp_dic[keyname]['min_l_pol'], max_l = exp_dic[keyname]['max_l_pol'])
    else:

        '''
        if old_code:
            F_mat = tools.fn_fisher_forecast_old(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
                min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
        else:
            F_mat_v1 = tools.fn_fisher_forecast_updated(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
                min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
            F_mat_v2 = tools.fn_fisher_forecast_with_trace(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
                min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
        '''
        '''
        F_mat = tools.fn_fisher_forecast(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
            min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
        F_mat_old = tools.fn_fisher_forecast_old(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
            min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])

        #F_mat_v2 = tools.fn_fisher_forecast_new(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
        #    min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
        '''

        #F_mat_v1 = tools.fn_fisher_forecast(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
        #    min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
        if expname.find('spt')>-1: #SPT CMB+lensing
            F_mat = tools.fn_fisher_forecast_with_trace_spt(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
                min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])
        else:
            F_mat = tools.fn_fisher_forecast_with_trace(els, Cl_deriv_dic, exp_delta_Cl_dic[keyname], param_names, pspectra_to_use = exp_dic[keyname]['pspectra_to_use'],\
                min_l_temp = exp_dic[keyname]['min_l_temp'], max_l_temp = exp_dic[keyname]['max_l_temp'], min_l_pol = exp_dic[keyname]['min_l_pol'], max_l_pol = exp_dic[keyname]['max_l_pol'])

        #from IPython import embed; embed()
        #sys.exit()



    if expname == 's4like_simple':

        opdic = {}
        opdic['exp_chars'] = exp_dic[keyname]
        opdic['param_names'] = param_names
        opdic['F_mat'] = F_mat
        #opdic['derived_params_derivatives'] = derived_param_names_dic_to_store
        ##opdic['Cov_mat'] = Cov_mat
        ##from IPython import embed; embed(); sys.exit()
        print('\n\t\t\t%s\n' %(opfname))
        np.save(opfname, opdic)
        sys.exit()

        from IPython import embed; embed()
        using_BW = 0
        using_JM = 1
        if using_BW:
            other_fisher_file = 'checks/s4_simple/fisher_matrices_basic.npy'
            other_fisher_dic = np.load(other_fisher_file, allow_pickle = 1).item()
            other_fisher_mat = other_fisher_dic['TTTEEE']['CAMB']['more accurate']
            other_fisher_params_ori = other_fisher_dic['parameter list']
            other_fisher_params_to_params_mapper_dic = {'ln10^{10}A_s':'As', 'N_ur':'neff', 'n_s':'ns', 'omega_b':'ombh2', 'omega_cdm':'omch2', 'tau_reio':'tau', '100*theta_s': 'thetastar'}#, 'YHe':}
            other_fisher_params = []
            for tmp in other_fisher_params_ori:
                other_fisher_params.append(other_fisher_params_to_params_mapper_dic[tmp])
            other_key = 'BW'

        if using_JM:
            if which_spectra == 'lensed_scalar':
                other_fisher_file = 'checks/s4_simple/fisher_matrices_basic_JW.npy'
            elif which_spectra == 'unlensed_scalar':
                other_fisher_file = 'checks/s4_simple/fisher_matrices_basic_JW_unlensed.npy'
            other_fisher_dic = np.load(other_fisher_file, allow_pickle = 1).item()
            other_fisher_mat = other_fisher_dic['TTEETE']
            other_fisher_params_ori = other_fisher_dic['params']
            other_fisher_params_to_params_mapper_dic = {'A_s':'As', 'tau':'tau', 'N_eff': 'neff', 'n_s':'ns', 'omega_b_h2':'ombh2', 'omega_c_h2':'omch2', 'theta_s': 'thetastar', 'mnu': 'mnu'}#, 'YHe':}
            other_fisher_params = []
            for tmp in other_fisher_params_ori:
                other_fisher_params.append(other_fisher_params_to_params_mapper_dic[tmp])
            other_key = 'JM'

            other_fisher_mat[:,2] = other_fisher_mat[:,2] / 1e9
            other_fisher_mat[2,:] = other_fisher_mat[2,:] / 1e9

        which_param = 'neff' ##'thetastar'#neff'##, 'omch2', 'tau' ##'As' #neff'

        F_mat_arr = [F_mat, other_fisher_mat]
        param_names_arr = [sorted(cosmo_param_dict.keys()), other_fisher_params]

        #F_mat_arr = [other_fisher_mat]
        #param_names_arr = [other_fisher_params]

        for fcntr in range(len(F_mat_arr)):
                F_mat_current = np.copy(F_mat_arr[fcntr])
                param_names = np.copy( param_names_arr[fcntr] )
                fix_Yhe = iter
                if fcntr == 0:
                    fix_params = ['ws', 'mnu', 'YHe', 'Alens', 'Aphiphi']
                else:
                    fix_params = []

                #if using_JM: fix_params.append('As')
                prior_dic = {'tau':0.007}

                #now fix parameters
                F_mat_current, param_names = tools.fn_fix_params(F_mat_current, param_names, fix_params)
                param_names = np.asarray(param_names)
                if fcntr == 0:
                    reqd_param_names = param_names


                #add prior
                F_mat_current = tools.fn_add_prior(F_mat_current, param_names, prior_dic)

                #get cov matrix now
                Cov_mat = sc.linalg.pinv2(F_mat_current) #made sure that COV_mat_l * Cinv_l ~= I

                param_names = np.asarray(param_names)
                pind = np.where(param_names == which_param)[0][0]
                pcntr1, pcntr2 = pind, pind
                cov_inds_to_extract = [(pcntr1, pcntr1), (pcntr1, pcntr2), (pcntr2, pcntr1), (pcntr2, pcntr2)]

                #cov_extract = np.asarray( [Cov_mat[ii] for ii in cov_inds_to_extract] ).reshape((2,2))
                cov_extract = []
                for ii in cov_inds_to_extract:
                    cov_extract.append( Cov_mat[ii] )
                cov_extract = np.asarray(cov_extract).reshape((2,2))

                sigma = cov_extract[0,0]**0.5
                sigma = round(sigma, 4)
                print(sigma, np.sqrt(np.diag(Cov_mat)))

        from IPython import embed; embed(); 

        #reqd_param_names = sorted(cosmo_param_dict.keys())
        F_mat_dic = {'SR': F_mat, other_key: other_fisher_mat}
        param_names_arr = [sorted(cosmo_param_dict.keys()), other_fisher_params]
        full_fisher_mat = {}
        for keycntr, keyname in enumerate( F_mat_dic ):
            F_mat_current = np.copy(F_mat_dic[keyname])
            param_names = np.copy( param_names_arr[keycntr] )
            F_mat_refined = []
            for pcntr1, p1 in enumerate( reqd_param_names ):
                pind1 = np.where(param_names == p1)[0][0]
                print(p1, pind1, keyname)
                for pcntr2, p2 in enumerate( reqd_param_names ):
                    if p1 not in other_fisher_params or p2 not in other_fisher_params: continue
                    pind2 = np.where(param_names == p2)[0][0]
                    F_mat_refined.append( (F_mat_current[pind2, pind1]) )

            totparams = len( reqd_param_names )
            F_mat_refined = np.asarray( F_mat_refined ).reshape( (totparams, totparams) )
            full_fisher_mat[keyname] = F_mat_refined

        full_fisher_dic_to_save = {}
        full_fisher_dic_to_save['fisher_mat'] = {}
        full_fisher_dic_to_save['fisher_mat'][other_key] = full_fisher_mat[other_key]
        if (0):
            full_fisher_dic_to_save['fisher_mat']['SR_with_Yp'] = full_fisher_mat['SR']
            full_fisher_dic_to_save['fisher_mat']['SR_no_Yp'] = full_fisher_mat_v2['SR']
        full_fisher_dic_to_save['fisher_mat']['SR'] = full_fisher_mat['SR']
        full_fisher_dic_to_save['params'] = reqd_param_names
        full_fisher_dic_to_save['exp_dic'] = exp_dic

        if using_JM:
            np.save('checks/s4_simple/fisher_comparison_JM_SR_lensed.npy', full_fisher_dic_to_save)
        elif using_BW:
            np.save('checks/s4_simple/fisher_comparison_BW_SR.npy', full_fisher_dic_to_save)

        sys.exit()



    #from IPython import embed; embed()
    #sign, logdetval = np.linalg.slogdet(F_mat)
    #logdetval = logdetval * sign
    #from IPython import embed; embed()
    if (1):##expname == 's4' or expname == 'planck' or expname.find('spt')>-1:
        opdic = {}
        opdic['exp_chars'] = exp_dic[keyname]
        opdic['param_names'] = param_names
        opdic['F_mat'] = F_mat
        opdic['param_dict'] = param_dict
        opdic['param_dict_derivatives'] = cosmo_param_dict
        ##opdic['Cov_mat'] = Cov_mat
        #from IPython import embed; embed(); #sys.exit()
        #print(F_mat)
        print('\n\t\t\t%s\n' %(opfname))
        np.save(opfname, opdic)        
        #sys.exit()
    else:

        #remove parameters that must be fixed    
        F_mat, param_names = tools.fn_fix_params(F_mat, param_names, fix_params)
        F_mat = np.mat( F_mat )# + np.eye(F_mat.shape[0])

        #adding prior
        if 'prior_dic' in exp_dic[keyname]:
            print('\t\t\tadding prior')
            F_mat = tools.fn_add_prior(F_mat, param_names, exp_dic[keyname]['prior_dic'])

        ##from IPython import embed; embed()
        #Cov_mat = np.linalg.pinv(F_mat)
        Cov_mat = sc.linalg.pinv2(F_mat) #made sure that COV_mat_l * Cinv_l ~= I

        F_dic[keyname] = F_mat
        Cov_dic[keyname] = Cov_mat

if (0):
    from IPython import embed; embed()
    from pylab import *
    F1 = F_dic['planck']
    F2 = F_dic['planck2']

    subplot(131);imshow(F1); colorbar();subplot(132);imshow(F2); colorbar();subplot(133);imshow(F1-F2); colorbar(); show();sys.exit()

    C1 = sc.linalg.pinv2(F1)
    C2 = sc.linalg.pinv2(F2)
#print Cov_dic
'''
opdic = {}
opdic['exp_dic'] = exp_dic
opdic['F_dic'] = F_dic
opdic['Cov_dic'] = Cov_dic
opfolder = 'results/'
if not os.path.exists('results'): os.system('mkdir results')

param_names_str = '-'.join(param_names)
opfname = '%s/%s-%scosmo.npy' %(opfolder, param_names_str,len(param_names))

np.save(opfname, opdic)
'''
if expname == 's4' or expname == 'planck':
    print('\n\n\n\t\t\tDOne\n\n\n\n')
    sys.exit()

if expname.find('spt')>-1 and which_spectra == 'lens_potential' and just_Aphiphi:
    Cov_mat = sc.linalg.pinv2(F_mat)
    sigma_Alens = np.sqrt(np.diag(Cov_mat))[0]

    if which_spectra == 'lens_potential':
        print('\n\tLensing SNR for %s = %.2f; Sigma Alens = %.3f; include_atm_noise = %s; lensing_est = %s; mod_lensing_N0_fac = %s\n' %(expname, 1./sigma_Alens, sigma_Alens, include_atm_noise, lensing_est, mod_lensing_N0_fac))
    else:
        from IPython import embed; embed()
        #pass
    sys.exit()

    from pylab import *
    from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
    from matplotlib.patches import Ellipse
    rcParams['font.family'] = 'serif'
    rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

    Nldic = exp_dic[expname]['Nldic']

    clf()
    phi_kappa_fac = 1e7 ##( (els * (els+1))**2. /4. )
    Cl_phiphi = Cl_dic['PP']
    deltaCl_phiphi = delta_Cl_dic['PP']

    nl_phiphi = Nldic['PP']
    #adjust lengths
    reclen_padding_zeros = max(els) + 1 - len(nl_phiphi)
    nl_phiphi = np.concatenate( (nl_phiphi, np.zeros(reclen_padding_zeros)) )
    nl_phiphi[nl_phiphi == 0] = max(nl_phiphi) * 1e6 #some large number
    nl_phiphi = nl_phiphi[els]

    ax = subplot(111, xscale = 'log')
    plot(els, Cl_phiphi * phi_kappa_fac)
    errorbar(els, Cl_phiphi * phi_kappa_fac, yerr = deltaCl_phiphi * phi_kappa_fac, color = 'black', ecolor = 'orangered', elinewidth = 1.)
    plot(els, nl_phiphi * phi_kappa_fac, color = 'darkred')
    ylabel(r'$C_{L}^{\kappa \kappa}$', fontsize = 14)
    xlabel(r'Multipole $L$', fontsize = 14)
    #ylim(1e-10,1e-4);
    ylim(-0.1, 2.)
    xlim(0, 6000)
    title(r'%s' %(expname))
    show()

    sys.exit()

sys.exit()
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

#make traingle plots
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
from matplotlib.patches import Ellipse
rcParams['font.family'] = 'serif'
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

tr, tc = len(param_names), len(param_names)
fig = figure(figsize=(tr + 3, tc + 3))
subplots_adjust(hspace = 0.2, wspace = 0.2)

totparams = len(param_names)
diag_matrix = np.arange( totparams**2 ).reshape((totparams, totparams)) + 1

sbpl_locs_dic = {}
for p1 in param_names:
    for p2 in param_names:
        sbpl_locs_dic[(p1,p2)] = cosmo_param_pl_chars_dict[p1][0] + ((cosmo_param_pl_chars_dict[p2][0]-1) * len(param_names))

for pcntr1, p1 in enumerate( param_names ):
    for pcntr2, p2 in enumerate( param_names ):        

        sbpl = sbpl_locs_dic[(p1,p2)]
        if sbpl not in np.tril(diag_matrix): continue

        if p1 in fix_params or p2 in fix_params: continue

        cov_inds_to_extract = [(pcntr1, pcntr1), (pcntr1, pcntr2), (pcntr2, pcntr1), (pcntr2, pcntr2)]
        x = param_dict[p1]
        y = param_dict[p2]

        deltax, epsilon_x = cosmo_param_dict[p1][1]
        deltay, epsilon_y = cosmo_param_dict[p2][1]

        if p1 == 'As': x*=1e9
        if p2 == 'As':  y*=1e9

        x1, x2 = x - deltax, x + deltax
        y1, y2 = y - deltay, y + deltay

        ax = subplot(tr, tc, sbpl)#, aspect = 'equal')

        p1str = cosmo_param_pl_chars_dict[p1][1]
        p2str = cosmo_param_pl_chars_dict[p2][1]

        if sbpl<=(tr*(tc-1)):
            setp(ax.get_xticklabels(), visible=False)
        else:
            xlabel(p1str, fontsize = 10);

        if ((sbpl-1)%tc == 0):
            ylabel(p2str, fontsize = 10);
        else:
            setp(ax.get_yticklabels(), visible=False)

        for label in ax.get_xticklabels(): label.set_fontsize(7)
        for label in ax.get_xticklabels(): label.set_fontsize(7)
        print(p1, p2, sbpl)
        for exp in expnamearr:

            exp_COV = Cov_dic[exp]

            cov_extract = np.asarray( [exp_COV[ii] for ii in cov_inds_to_extract] ).reshape((2,2))
            #marginalised_F = sc.linalg.pinv2(cov_extract)
            #cov_extract = sc.linalg.pinv2(marginalised_F)

            #from IPython import embed; embed()

            colorarr = [exp_dic[exp]['color']]
            alphaarr = [1., 0.5]
            for ss in range(0,1):

                lwval = 0.75
                if p1 == p2:

                    widthval = cov_extract[0,0]**0.5##/2.35
                    hor, ver = tools.fn_get_Gaussian(x, widthval, x1, x2, epsilon_x)
                    plot(hor, ver, color = colorarr[ss], lw = lwval, label = r'%g' %(widthval))
                    legend(loc = 3, framealpha = 1, fontsize = 5, edgecolor = 'None')

                    xlim(x1, x2)
                    ylim(0., 1.)
                    setp(ax.get_yticklabels(), visible=False); tick_params(axis='y',left='off')
                    title(p1str, fontsize = 10);

                else:

                    Ep = tools.fn_get_ellipse_specs(cov_extract, howmanysigma = ss + 1)
                    widthval, heightval = Ep[0], Ep[1]
                    ellipse = Ellipse(xy=[x,y], width=2.*widthval, height=2.*heightval, angle=np.degrees(Ep[2]))

                    ax.add_artist(ellipse)
                    ellipse.set_clip_box(ax.bbox)
                    ellipse.set_facecolor('None')#colorarr[ss])
                    ellipse.set_edgecolor(colorarr[ss])
                    ellipse.set_linewidth(lwval)
                    #ellipse.set_alpha(alphaarr[ss])

                    xlim(x1, x2)
                    ylim(y1, y2)


#make legends
ax = subplot(tr, tc, 5)
for exp in expnamearr:
    #colorval = exp_pl_chars_dic[exp][0]
    #expnamespl = exp_pl_chars_dic[exp][1]
    colorval = exp_dic[exp]['color']
    expnamespl = exp_dic[exp]['label']
    plot([],[], color = colorval, label = r'%s' %(expnamespl))
axis('off')
legend(loc = 1, fontsize = 10, framealpha = 0)
#show();quit()
#savefig('plots/Neff_fisher_constraints_planck_S4.pdf')

#plname = 'plots/fisher_constraints_planck_smica_S4_%scosmoparams_lmax%s.pdf' %(len(param_names), max_l)
expnames_str = '-'.join(expnamearr)
param_names_str = '-'.join(param_names)
#pspectra_to_use_str = '-'.join(pspectra_to_use)
#plname = 'plots/planck_S4_%s-%scosmo_%sspec_lmax%s.pdf' %(param_names_str, len(param_names), pspectra_to_use_str, max_l)
plname = 'plots/%s-%scosmo_%sdeltaelbinning_%s.pdf' %(param_names_str, len(param_names), delta_l, expnames_str)


savefig(plname)
sys.exit()










