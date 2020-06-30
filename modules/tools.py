import numpy as np, camb, sys, scipy as sc, os
#from camb import model, initialpower
#from pylab import *
from copy import deepcopy


########################################################################################################################

def fn_ini_param_dic(fpath = 'params/params_planck_r_0.0_2015_cosmo_lensed_LSS.txt'):
    """
    read params file and initialise cosmology
    """
    try:
        params = np.recfromtxt(fpath, delimiter = '=', encoding = 'utf-8')
    except:
        params = np.recfromtxt(fpath, delimiter = '=')
    param_dict = {}
    for rec in params:
        val = rec[1].strip()##.decode("utf-8")
        try:
            if val.find('.')>-1:
                val = float(val)
            else:
                val = int(val)
        except:
            val = str(val)

        if val == 'None':
            val = None
        paramname = rec[0].strip()#.decode("utf-8")
        param_dict[paramname] = val

    return param_dict

########################################################################################################################

def fn_get_Nl(els, fwhm, rms_map_T, rms_map_P = None, CMB = 1):
    """
    compute Nl - white noise + beam
    """

    if rms_map_P == None:
        rms_map_P = rms_map_T * 1.414

    fwhm_radians = np.radians(fwhm/60.)
    #Bl = np.exp((-fwhm_radians**2.) * els * (els+1) /2.35)
    sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
    sigma2 = sigma ** 2
    Bl = np.exp(els * (els+1) * sigma2)

    rms_map_T_radians = rms_map_T * np.radians(1/60.)
    rms_map_P_radians = rms_map_P * np.radians(1/60.)

    Nl_TT = (rms_map_T_radians)**2. * Bl
    Nl_PP = (rms_map_P_radians)**2. * Bl

    return Bl, Nl_TT, Nl_PP

########################################################################################################################

def fn_pad(el, cl):
    reclen_padding_zeros = max(el) - len(cl) + 1
    cl = np.concatenate( (cl, np.zeros(reclen_padding_zeros)) )
    cl[cl == 0] = max(cl) * 1e6 #some large number        
    return cl

########################################################################################################################

def fn_delta_Cl(els, Cl_dic, exp_dic, delta_l = 1., Nlfile = None, Nldic = None):
    if Nldic is not None:
        if 'T' in Nldic:
            Nl_TT = Nldic['T'][els]
        if 'P' in Nldic:
            Nl_PP = Nldic['P'][els]
        if 'TT' in Nldic:
            Nl_TT = fn_pad(els, Nldic['TT'])
            Nl_TT = Nl_TT[els]

            Nl_TT = Nldic['TT']

        if 'EE' in Nldic:
            Nl_EE = fn_pad(els, Nldic['EE'])
            Nl_EE = Nl_EE[els]

            Nl_EE = Nldic['EE']

        if 'TE' in Nldic:
            Nl_TE = Nldic['TE']
            if Nldic['TE'] is not None:
                Nl_TE = fn_pad(els, Nl_TE)
                Nl_TE = Nl_TE[els]

            Nl_TE = Nldic['TE']

        if 'TB' in Nldic:
            Nl_TB = fn_pad(els, Nldic['TB'])
            Nl_TB = Nl_TB[els] 

            Nl_TB = Nldic['TB']           

        if 'BB' in Nldic:
            Nl_BB = fn_pad(els, Nldic['BB'])
            Nl_BB = Nl_BB[els]            

            Nl_BB = Nldic['BB']

        if 'EB' in Nldic:
            Nl_EB = fn_pad(els, Nldic['EB'])
            Nl_EB = Nl_EB[els]            

            Nl_EB = Nldic['EB']

        if 'PP' in Nldic: #lensing phiphi N0
            Nl_PP = fn_pad(els, Nldic['PP'])            
            Nl_PP = Nl_PP[els] 
            ##Nl_PP = Nldic['PP']


        if (0):
            loglog(Nl_TT)
            loglog(Nl_EE)
            ylim(1e-6, 1e3)
            show();sys.exit()

    else:
        if Nlfile is not None:
            #from IPython import embed; embed()        
            Nl = np.load(Nlfile)

            if isinstance(Nl, dict):
                Nl_TT = Nl['TT']
                Nl_PP = Nl['PP']
            else:
                Nl_TT = Nl
                Nl_PP = Nl_TT * np.sqrt(2.)

            #adjust lengths
            reclen_padding_zeros = max(els) + 1 - len(Nl_TT)
            Nl_TT = np.concatenate( (Nl_TT, np.zeros(reclen_padding_zeros)) )
            Nl_PP = np.concatenate( (Nl_PP, np.zeros(reclen_padding_zeros)) )

            #pick the desired ell
            Nl_TT = Nl_TT[els]
            Nl_PP = Nl_PP[els]
        else:
            exp = [exp_dic['beam'], exp_dic['deltaT']] #beam, deltaT, deltaP = None --> 1.414 * deltaT
            Bl, Nl_TT, Nl_EE = fn_get_Nl(els, *exp)
            Nl_TE = np.zeros( len(els) )

    #exp = [exp_dic['beam'], exp_dic['deltaT']] #beam, deltaT, deltaP = None --> 1.414 * deltaT
    #Bl, Nl_TT_white, Nl_EE_white = fn_get_Nl(els, *exp)

    delta_Cl_dic = {}
    for XX in Cl_dic:
        #print(XX)       
        if XX == 'TT':
            Nl = Nl_TT
        elif XX == 'EE' or XX == 'BB':
            Nl = Nl_EE ##Nl_PP
        elif XX == 'TE':
            Nl = Nl_TE

            if (1):
                print('\n\n\n\t\t\t\tnulling galaxy TE\n\n\n')
                Nl = np.copy(Nl) * 0.

        elif XX == 'PP':
            Nl = Nl_PP
        else:
            Nl = np.zeros( len( els ) )

        #print(XX, Nl)
        if Nl is None:
            Nl = np.zeros( len( els ) )
        #print(XX, Nl)

        Cl = Cl_dic[XX]
        if (0): #20200602 - this is wrong - being done twice. Nl was already picked at binned el centres.
            els_2, Nl = fn_el_binning(np.arange(min(els),max(els)+1), Nl, delta_l = delta_l)
            #els_2, Nl = fn_el_binning(np.arange(min(els),len(Nl)+2), Nl, delta_l = delta_l)

        ##delta_Cl_dic[XX] = np.sqrt(2./ (2*els + 1) / exp_dic['fsky'] / (1.*delta_l) ) * (Cl + Nl)
        delta_Cl_dic[XX] = np.sqrt(2./ (2.*els + 1.) / exp_dic['fsky'] ) * (Cl + Nl)

        #from IPython import embed; embed()

        if (0):
            ax = subplot(111, yscale = 'log')
            plot(Nl);
            errorbar(els, Cl, yerr = delta_Cl_dic[XX]);
            title(XX);
            ylim(1e-8, 1e3)
            show(); #sys.exit()

    if (0):
        ax = subplot(111, yscale = 'log')
        dls_fac = 1. ##(els * (els+1))/2/np.pi
        plot(Cl_dic['TT'] * dls_fac, 'k-'); plot(Cl_dic['EE'] * dls_fac, 'r-'); #plot(Cl_dic['TE'] * dls_fac, 'g-')
        plot(Nl_TT * dls_fac, 'k--'); plot(Nl_EE * dls_fac, 'r--'); #plot(Nl_TE * dls_fac, 'g--')
        plot(Nl_TT_white * dls_fac, 'k:'); plot(Nl_EE_white * dls_fac, 'r:'); #plot(Nl_TE * dls_fac, 'g:')
        #plot(delta_Cl_dic['TT'] * dls_fac, 'k--'); plot(delta_Cl_dic['EE'] * dls_fac, 'r--'); #plot(delta_Cl_dic['TE'] * dls_fac, 'g--')
        show()
        sys.exit()

    #print(delta_Cl_dic); sys.exit()
    return delta_Cl_dic

########################################################################################################################

def fn_fix_params(F_mat, param_names, fix_params):

    #remove parameters that must be fixed    
    F_mat_refined = []
    for pcntr1, p1 in enumerate( param_names ):
        for pcntr2, p2 in enumerate( param_names ):
            if p1 in fix_params or p2 in fix_params: continue
            F_mat_refined.append( (F_mat[pcntr2, pcntr1]) )

    totparamsafterfixing = int( np.sqrt( len(F_mat_refined) ) )
    F_mat_refined = np.asarray( F_mat_refined ).reshape( (totparamsafterfixing, totparamsafterfixing) )

    param_names_refined = []
    for p in param_names:
        if p in fix_params: continue
        param_names_refined.append(p)


    return F_mat_refined, param_names_refined

########################################################################################################################

def get_jacobian(param_dict, param_dict_derivatives, param_names, derived_param_names_dic, derived_param_names = None, der_params_derv_fname = None):

    if (0):
        print('\n\n\tremoving non linear lensing\n')
        param_dict['lens_potential_accuracy'] = 0

    if (1):
        #print('\n\n\tusing cosmomc_theta now\n')
        param_dict['cosmomc_theta'] = 0.010406810462970357

    param_dict_derivatives_refined = {}
    for p in param_names:
        param_dict_derivatives_refined[p] = param_dict_derivatives[p]

    if os.path.exists(der_params_derv_fname):
        der_params_derv_dic = np.load(der_params_derv_fname, allow_pickle = 1).item()
    else:
        der_params_derv_dic = {}

    if derived_param_names is None:
        derived_param_names = sorted( derived_param_names_dic.keys() )
        derived_param_names = np.concatenate( (param_names, derived_param_names) )

    Cl_dic = None
    nparams = len( param_names )
    nderparams = len( derived_param_names )
    J_mat = np.zeros( (nparams, nderparams) )
    for pcntr, der_p in enumerate( derived_param_names ):

        if 'thetastar' in param_dict_derivatives:
            use_thetastar = 1
        else:
            use_thetastar = 0

        if 'cosmomc_theta' in param_dict_derivatives:
            use_cosmomc_theta = 1
        else:
            use_cosmomc_theta = 0

        #check if this is fine for all parameters
        if (0):##der_p.find('omega')>-1 or der_p.find('sigma8')>-1 or der_p == 'cosmomc_theta':
            use_thetastar = 0 ###1
            #use_cosmomc_theta = 1

        print(der_p, use_thetastar, use_cosmomc_theta)

        if der_p not in der_params_derv_dic:
            if der_p in param_names:
                derived_param_deriv = {}
                for tmp_p in param_names:
                    if tmp_p == der_p:
                        derived_param_deriv[tmp_p] = 1.
                    else:
                        derived_param_deriv[tmp_p] = 0.
            else:
                der_p_cmd = derived_param_names_dic[der_p]
                derived_param_deriv = fn_get_dervatives(Cl_dic, param_dict, param_dict_derivatives_refined, use_thetastar = use_thetastar, use_cosmomc_theta = use_cosmomc_theta, derived_param = [der_p, der_p_cmd])
            der_params_derv_dic[der_p] = [derived_param_deriv, use_thetastar]
            np.save(der_params_derv_fname, der_params_derv_dic)
        else:
            derived_param_deriv = der_params_derv_dic[der_p][0]

        for pcntr2, p2 in enumerate(param_names):
            if derived_param_deriv[p2] == 0.:
                J_mat[pcntr2, pcntr] = derived_param_deriv[p2]
            else:
                J_mat[pcntr2, pcntr] = 1./derived_param_deriv[p2]

        '''
        for pcntr2, p2 in enumerate(param_names):
            J_mat[pcntr2, pcntr] = derived_param_deriv[p2]
        J_mat = sc.linalg.pinv2(J_mat)
        '''

    from IPython import embed; embed()
    return J_mat

########################################################################################################################

def fisher_transform(F_mat, param_dict, param_dict_derivatives, param_names, derived_param_names_dic, J_mat = None):

    F_mat = np.mat( F_mat )
    if J_mat is None:
        J_mat = get_jacobian(param_dict, param_dict_derivatives, param_names, derived_param_names_dic)
    J_mat = np.mat( J_mat )

    F_mat_updated = np.dot(J_mat.T, np.dot( F_mat, J_mat ) )

    #from IPython import embed; embed()
    
    return F_mat_updated

########################################################################################################################

def fn_add_prior(F_mat, param_names, prior_dic):

    for pcntr1, p1 in enumerate( param_names ):
        for pcntr2, p2 in enumerate( param_names ):
            if p1 == p2 and p1 in prior_dic:
                prior_val = prior_dic[p1]
                F_mat[pcntr2, pcntr1] += 1./prior_val**2.

    return F_mat

########################################################################################################################

def get_planck_prior(param_names, fix_params = [], add_BAO = 1):

    if add_BAO:
        planck_prior_dic_stored = { 
                'ombh2' : 0.00014, 'omch2' : 0.00091, 
                'h': 0.0042, 'tau' : 0.0071, 'ns' : 0.0038,
                'As' : 0.2917e-10 * 1e9, 
                #'ws' : (-1e-1,(0.5, 0.0001)), 'neff': (1e-3,(0.5,0.01)), 'mnu': (1e-3,(1.,0.0001)), 'YHe': (1e-3,(1.,0.0001)),\
                #'Alens': (1e-3, (1., 0.3)), 
                }
    else:
        planck_prior_dic_stored = { 
                'ombh2' : 0.00015, 'omch2' : 0.0012, 
                'h': 0.0054, 'tau' : 0.0073, 'ns' : 0.0042,
                'As' : 0.2917e-10 * 1e9, 
                }

    planck_prior_dic = {}
    for p1 in param_names:
        if p1 in fix_params: continue
        if p1 not in planck_prior_dic_stored: continue
        planck_prior_dic[p1] = planck_prior_dic_stored[p1]

    F_mat_planck = np.zeros( ( len(param_names), len(param_names) ) )
    F_mat_planck = fn_add_prior(F_mat_planck, param_names, planck_prior_dic)

    return planck_prior_dic, F_mat_planck

########################################################################################################################

def fn_get_ellipse_specs(COV, howmanysigma = 1):
    """
    Refer https://arxiv.org/pdf/0906.4123.pdf
    """
    assert COV.shape == (2,2)
    confsigma_dic = {1:2.3, 2:6.17, 3: 11.8}

    sig_x2, sig_y2 = COV[0,0], COV[1,1]
    sig_xy = COV[0,1]
    
    t1 = (sig_x2 + sig_y2)/2.
    t2 = np.sqrt( (sig_x2 - sig_y2)**2. /4. + sig_xy**2. )
    
    a2 = t1 + t2
    b2 = t1 - t2

    a = np.sqrt(a2)
    b = np.sqrt(b2)

    t1 = 2 * sig_xy
    t2 = sig_x2 - sig_y2
    theta = np.arctan2(t1,t2) / 2.
    
    alpha = np.sqrt(confsigma_dic[howmanysigma])
    
    #return (a*alpha, b*alpha, theta)
    return (a*alpha, b*alpha, theta, alpha*(sig_x2**0.5), alpha*(sig_y2**0.5))

########################################################################################################################

def fn_get_Gaussian(mean, sigma, minx, maxx, delx):

        x = np.arange(minx, maxx, delx)

        #return x, 1./(2*np.pi*sigma)**0.5 * np.exp( -(x - mean)**2. / (2 * sigma**2.)  )
        return x, np.exp( -(x - mean)**2. / (2 * sigma**2.)  )

########################################################################################################################

def fn_fisher_forecast_Aphiphi(els, Cl_deriv_dic, delta_Cl_dic, params, pspectra_to_use, min_l = 0, max_l = 6000):

    npar = len(params)
    F = np.zeros([npar,npar])
    #els = np.arange( len( delta_Cl_dic.values()[0] ) )

    pspectra_to_use_full = np.asarray( ['PP'] )

    for lcntr, l in enumerate( els ):

        if l<min_l or l>max_l:
            continue

        PP = delta_Cl_dic['PP'][lcntr]
        COV_mat_l = PP**2.
        COV_mat_l = np.mat( COV_mat_l )
        Cinv_l = sc.linalg.pinv2(COV_mat_l) #made sure that COV_mat_l * Cinv_l ~= I
        #print l, p, p2, fprime1_l_vec, fprime2_l_vec, COV_mat_l

        pspec_combinations = []
        for X in pspectra_to_use:
            for Y in pspectra_to_use:
                xind = np.where(pspectra_to_use_full == X)[0][0]
                yind = np.where(pspectra_to_use_full == Y)[0][0]
                if [Y,X, yind, xind] in pspec_combinations: continue
                pspec_combinations.append([X, Y, xind, yind])

        param_combinations = []
        for pcnt,p in enumerate(params):
            for pcnt2,p2 in enumerate(params):
                ##if [p2,p,pcnt2,pcnt] in param_combinations: continue
                param_combinations.append([p,p2, pcnt, pcnt2])

        for (p,p2, pcnt, pcnt2) in param_combinations:
            for (X,Y, xind, yind) in pspec_combinations:

                der1 = np.asarray( [Cl_deriv_dic[p]['PP'][lcntr]] )
                der2 = np.asarray( [Cl_deriv_dic[p2]['PP'][lcntr]] )

                fprime1_l_vec = np.zeros(len(der1))
                fprime2_l_vec = np.zeros(len(der2))

                fprime1_l_vec[xind] = der1[xind]
                fprime2_l_vec[yind] = der2[yind]

                #if l > 100:
                #    from IPython import embed; embed()

                curr_val = np.dot(fprime1_l_vec, np.dot( Cinv_l, fprime2_l_vec ))

                F[pcnt2,pcnt] += curr_val

    return F    

########################################################################################################################
########################################################################################################################
########################################################################################################################

#def fn_fisher_forecast(els, Cl_dic, Cl_deriv_dic, delta_Cl_dic, params, pspectra_to_use, min_l_temp = None, max_l_temp = None, min_l_pol = None, max_l_pol = None):
def fn_fisher_forecast_with_trace(els, Cl_deriv_dic, delta_Cl_dic, params, pspectra_to_use, min_l_temp = None, max_l_temp = None, min_l_pol = None, max_l_pol = None):

    if min_l_temp is None: min_l_temp = 0
    if max_l_temp is None: max_l_temp = 10000

    if min_l_pol is None: min_l_pol = 0
    if max_l_pol is None: max_l_pol = 10000

    npar = len(params)
    F = np.zeros([npar,npar])
    #els = np.arange( len( delta_Cl_dic.values()[0] ) )

    with_lensing = 0
    if 'PP' in pspectra_to_use:
        with_lensing = 1

    all_pspectra_to_use = []
    for tmp in pspectra_to_use:
        if isinstance(tmp, list):      
            all_pspectra_to_use.extend(tmp)
        else:
            all_pspectra_to_use.append(tmp)

    print('!!Using trace now!!')

    for lcntr, l in enumerate( els ):

        '''
        #creating the covariance matrix for this multipole
        TT, EE, TE = delta_Cl_dic['TT'][lcntr], delta_Cl_dic['EE'][lcntr], delta_Cl_dic['TE'][lcntr]
        if with_lensing:
            Tphi, Ephi, PP = delta_Cl_dic['Tphi'][lcntr], delta_Cl_dic['Ephi'][lcntr], delta_Cl_dic['PP'][lcntr]
        else:
            Tphi = Ephi = PP = 0.
        '''

        TT, EE, TE = 0., 0., 0.
        Tphi = Ephi = PP = 0.
        if 'TT' in delta_Cl_dic:
            TT = delta_Cl_dic['TT'][lcntr]
        if 'EE' in delta_Cl_dic:
            EE = delta_Cl_dic['EE'][lcntr]
        if 'TE' in delta_Cl_dic:
            TE = delta_Cl_dic['TE'][lcntr]
        if with_lensing:
            Tphi, Ephi, PP = delta_Cl_dic['Tphi'][lcntr], delta_Cl_dic['Ephi'][lcntr], delta_Cl_dic['PP'][lcntr]

        null_TT, null_EE, null_TE = 0, 0, 0
        if l<min_l_temp or l>max_l_temp:
            null_TT = 1
        if l<min_l_pol or l>max_l_pol: 
            null_EE = 1
            null_TE = 1
        null_PP = 0 #Lensing noise curves already have pretty large noise outside desired L range
        #if l<min_l_TE or l>max_l_TE:  
        #    null_TE = 1

        #20200611
        if 'TT' not in all_pspectra_to_use:
            null_TT = 1
        if 'EE' not in all_pspectra_to_use:
            null_EE = 1
        if 'TE' not in all_pspectra_to_use:
            if 'TT' not in pspectra_to_use and 'EE' not in pspectra_to_use:
                null_TE = 1
        if ['TT', 'EE', 'TE'] in pspectra_to_use:
            null_TT = 0
            null_EE = 0
            null_TE = 0
        #20200611

        ##if (null_TT and null_EE and null_TE): continue# and null_PP): continue
        if (null_TT and null_EE and null_TE and null_PP): continue

        param_combinations = []
        for pcnt,p in enumerate(params):
            for pcnt2,p2 in enumerate(params):
                ##if [p2,p,pcnt2,pcnt] in param_combinations: continue
                param_combinations.append([p,p2, pcnt, pcnt2])

        def get_cov(TT, EE, TE, PP, TP, EP):

            C = np.zeros( (3,3) ) #TT, EE, PP
            C[0,0] = TT
            C[1,1] = EE
            C[0,1] = C[1,0] = TE

            C[2,2] = PP
            C[0,2] = C[2,0] = TP
            C[1,2] = C[2,1] = 0. ##EP

            return np.mat( C )

        #20200621 - nulling unwanted fields
        if null_TT and null_TE: TT = 0
        if null_EE and null_TE: EE = 0
        if null_TE and (null_TT and null_EE): TE = 0
        if null_PP: PP = Tphi = EPhi = 0
        if null_TT: Tphi = 0
        if null_EE: Ephi = 0
        #20200621 - nulling unwanted fields

        COV_mat_l = get_cov(TT, EE, TE, PP, Tphi, Ephi)
        inv_COV_mat_l = sc.linalg.pinv2(COV_mat_l)

        #print(l, null_TT, null_EE, null_EE)

        if (0):##l%500 == 0: 
            from IPython import embed; embed()
            print(l, null_TT, null_EE, null_TE, null_PP)
            print(COV_mat_l)

        for (p,p2, pcnt, pcnt2) in param_combinations:

            '''
            TT_der1, EE_der1, TE_der1 = Cl_deriv_dic[p]['TT'][lcntr], Cl_deriv_dic[p]['EE'][lcntr], Cl_deriv_dic[p]['TE'][lcntr]
            TT_der2, EE_der2, TE_der2 = Cl_deriv_dic[p2]['TT'][lcntr], Cl_deriv_dic[p2]['EE'][lcntr], Cl_deriv_dic[p2]['TE'][lcntr]
            '''

            TT_der1, EE_der1, TE_der1 = 0., 0., 0.
            TT_der2, EE_der2, TE_der2 = 0., 0., 0.

            if 'TT' in Cl_deriv_dic[p]:
                TT_der1 = Cl_deriv_dic[p]['TT'][lcntr]
                TT_der2 = Cl_deriv_dic[p2]['TT'][lcntr]
            if 'EE' in Cl_deriv_dic[p]:
                EE_der1 = Cl_deriv_dic[p]['EE'][lcntr]
                EE_der2 = Cl_deriv_dic[p2]['EE'][lcntr]
            if 'TE' in Cl_deriv_dic[p]:
                TE_der1 = Cl_deriv_dic[p]['TE'][lcntr]
                TE_der2 = Cl_deriv_dic[p2]['TE'][lcntr]


            if with_lensing:
                PP_der1, TPhi_der1, EPhi_der1 = Cl_deriv_dic[p]['PP'][lcntr], Cl_deriv_dic[p]['Tphi'][lcntr], Cl_deriv_dic[p]['Ephi'][lcntr]
                PP_der2, TPhi_der2, EPhi_der2 = Cl_deriv_dic[p2]['PP'][lcntr], Cl_deriv_dic[p2]['Tphi'][lcntr], Cl_deriv_dic[p2]['Ephi'][lcntr]
            else:
                PP_der1 = PP_der2 = 0.
                TPhi_der1 = TPhi_der2 = 0. 
                EPhi_der1 = EPhi_der2 = 0.

            '''
            if null_TT: TT_der1 = TT_der2 = 0.
            if null_EE: EE_der1 = EE_der2 = 0.
            if null_TE: TE_der1 = TE_der2 = 0
            if null_PP: PP_der1 = PP_der2 = 0
            '''

            if null_TT: TT_der1 = TT_der2 = TPhi_der1 = TPhi_der2 = 0
            if null_EE: EE_der1 = EE_der2 = EPhi_der1 = EPhi_der2 = 0
            if null_TE: TE_der1 = TE_der2 = 0
            if null_PP: PP_der1 = PP_der2 = 0

            fprime1_l_vec = get_cov(TT_der1, EE_der1, TE_der1, PP_der1, TPhi_der1, EPhi_der1)
            fprime2_l_vec = get_cov(TT_der2, EE_der2, TE_der2, PP_der2, TPhi_der2, EPhi_der2)

            if (0):##l==5000:#l % 500 == 0:
                #from IPython import embed; embed()
                print(null_TT, null_EE, null_TE, null_PP, with_lensing)
                print(fprime1_l_vec)
                print(fprime2_l_vec)
                print(COV_mat_l) 

            #curr_val = np.trace( inv_COV_mat_l * fprime1_l_vec * inv_COV_mat_l * fprime2_l_vec)
            curr_val = np.trace( np.dot( np.dot(inv_COV_mat_l, fprime1_l_vec), np.dot(inv_COV_mat_l, fprime2_l_vec) ) )

            F[pcnt2,pcnt] += curr_val

    if (0):
        from IPython import embed; embed()
        F = np.mat(F)

    if (0):
        from IPython import embed; embed()
        F = np.mat(F)
        C = sc.linalg.pinv2(F) #made sure that COV_mat_l * Cinv_l ~= I

    return F   

########################################################################################################################
