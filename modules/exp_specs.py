import numpy as np
from pylab import *
#import py_ini


def get_lensing_noise(nlfile, mod_lensing_N0_fac = 1., lensing_est = 'MV', expname = None):
    result_dic = np.load(nlfile, allow_pickle = 1, encoding = 'latin1').item()
    lmin, lmax, lmax_tt = result_dic['lmin'], result_dic['lmax'], result_dic['lmax_tt']
    el, cl_kk = result_dic['els'], result_dic['cl_kk']

    min_L, max_L = 100, 5000#len(cl_kk)

    if expname == 'sptpol' or expname == 'sptsz':
        min_L, max_L = 100, 2000 ##2000
        min_l_temp, max_l_temp = min_L, max_L
        min_l_pol, max_l_pol = min_L, max_L
        print('\n\tpicking only %s to %s for lensing multipoles\n' %(min_L, max_L))

    '''
    nl_tt, nl_eb, nl_mv, nl_mvpol = result_dic['Nl_TT'].real, result_dic['Nl_EB'].real, result_dic['Nl_MV'].real, result_dic['Nl_MVpol'].real

    """
    Dls_fac = (el * (el+1))**2. /4.
    nl_tt /= Dls_fac
    nl_mvpol /= Dls_fac
    nl_mv /= Dls_fac
    nl_eb /= Dls_fac
    """

    nl_tt[np.isnan(nl_tt)] = 0.
    nl_mvpol[np.isnan(nl_mvpol)] = 0.
    nl_eb[np.isnan(nl_eb)] = 0.
    nl_mv[np.isnan(nl_mv)] = 0.

    Nldic = {}
    Nldic['TT'] = nl_tt
    Nldic['MV'] = nl_mv
    Nldic['MVpol'] = nl_mvpol
    Nldic['EB'] = nl_eb
    '''

    nl_name_str = 'Nl_%s' %(lensing_est)            
    nl_pp = result_dic[nl_name_str].real
    nl_pp[np.isnan(nl_pp)] = 1e10
    nl_pp[el>max_L] = 1e10
    nl_pp[el<min_L] = 1e10
    #semilogy(nl_pp); ylim(1e-9, 1e-4); xlim(0, 5000); show(); sys.exit()

    if (0):
        Dls_fac = (el * (el+1))**2. /2./np.pi
        nl_pp /= Dls_fac
    
    if mod_lensing_N0_fac != 1.:
        print('\n\tModifying lensing noise by x%s\n' %(mod_lensing_N0_fac))
        nl_pp *= mod_lensing_N0_fac
    return nl_pp, min_L, max_L

def get_Nldic(Nlfile, els):
    dic = np.load(Nlfile, allow_pickle = 1, encoding = 'latin1').item()
    el_nl, cl_residual = dic['el'], dic['cl_residual']

    if 'T' in cl_residual:
        Nl_TT, Nl_EE = cl_residual['T'], cl_residual['P']
        Nl_TT = np.interp(els, el_nl, Nl_TT)
        Nl_EE = np.interp(els, el_nl, Nl_EE)
        Nl_TE = None
    else:
        Nl_TT, Nl_EE = cl_residual['TT'], cl_residual['EE']
        if 'TE' in cl_residual:
            Nl_TE = cl_residual['TE']
        else:
            Nl_TE = None
        el_nl = np.arange(len(Nl_TT))
        Nl_TT = np.interp(els, el_nl, Nl_TT)
        Nl_EE = np.interp(els, el_nl, Nl_EE)
        if Nl_TE is not None:
            Nl_TE = np.interp(els, el_nl, Nl_TE)
        else:
            Nl_TE = np.zeros(len(els))

    Nldic = {}
    Nldic['TT'] = Nl_TT
    Nldic['EE'] = Nl_EE
    Nldic['TE'] = Nl_TE

    return Nldic

def add_planck_via_inv_variance(Nldic, Nldic_for_inv_var, include_planck_on_large_scales):

    extra_nl_dic_str = None
    if Nldic_for_inv_var is not None:
        Nldic_ref = {}
        for XX in Nldic_for_inv_var:
            nl1, nl2 = Nldic[XX], Nldic_for_inv_var[XX]

            if np.sum(nl2) == 0: #nothing to add here
                Nldic_ref[XX] = nl1
                continue

            padzeros_len = len(nl1) - len(nl2)
            nl2_pad = np.zeros(padzeros_len) + 1e10
            nl2 = np.concatenate( (nl2, nl2_pad) )
            nl_ref = 1./ ( (1./nl1) + (1./nl2) )

            Nldic_ref[XX] = nl_ref

            if (0):
                ax = subplot(111, yscale = 'log')
                plot(nl1,'red'); plot(nl2, 'green'); plot(nl_ref, 'black'); title(r'%s' %XX); ylim(1e-6, 1e1); show()

        extra_nl_dic_str = ''
        if include_planck_on_large_scales:
            min_l_temp_arr_ori = [2.]
            min_l_pol_arr_ori = [2.]
            extra_nl_dic_str = 'withplanckonlargescales'

    for XX in Nldic:
        if XX in Nldic_ref:
            Nldic[XX] = Nldic_ref[XX]

    return Nldic, min_l_temp_arr_ori, min_l_pol_arr_ori, extra_nl_dic_str
        
def define_exp_specs(expname, which_spectra = None, els = None, Nlfile = None, include_gal = 1, galmask = None, totlf = 2, totmf = 12, tothf = 5, also_te = 1, force_fsky_val = -1., \
    lensing_est = 'MV', Nldic_for_inv_var = None, include_planck_on_large_scales = 0, mod_lensing_N0_fac = 1., add_lensing = 0):
    if expname == 's4' or expname == 's4wide' or expname == 'planck': #clean this bit later
        exp_dic = {}
        #expnamearr = ['planck', 's4']
        #expnamearr = ['s4']
        expnamearr = [expname]
        for expname in expnamearr:
            '''
            if expname == 'planck':
                pspectra_to_use_arr = [['TT','EE','TE']]
                #Nlfile_arr = ['data/Planck/DR2/Nl_smica_halfmission_HDHM.npy']
                Nlfile = 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy'
            '''
            extrastr_tubesmoved = None
            Nldic = {}
            if expname == 's4':
                pspectra_to_use_arr = [['TT'], ['EE']]#, ['TE'], ['EE','TE']]#, ['TT','EE','TE']]
                #pspectra_to_use_arr = [['TT'], ['EE'], ['TE'], ['EE','TE'], ['TT','EE','TE']]
                #pspectra_to_use_arr = [['TT'], ['EE'], ['TE']]
                #pspectra_to_use_arr = [['EE','TE']]
                #pspectra_to_use_arr = [['EE'], ['TE']]
                #pspectra_to_use_arr = [['TE']]
                pspectra_to_use_arr = [['EE']]
                #Nlfile_arr = [None]
                
                #Nlfile_for_S4 = 'S4_ilc_20203030.npy'

                #s4_results after inlcuding gaalctic sims
                #s4_ilc_folder = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/scripts/notebooks/results/galactic_sims/s4like_mask/'
                #s4_ilc_folder = '/Users/sraghunathan/Research/SPTPol/analysis/git/DRAFT/scripts/notebooks/results/galactic_sims/s4like_mask/tt_ee_te/20200528/with_gal/'
                #Nlfile_for_S4 = '%s/S4_ilc_zonca_sims_20204028_galaxy0_93-145-225-278_TT-EE_AZ.npy' %(s4_ilc_folder)

                #galmask = 0
                #Nlfile_for_S4 = '%s/S4_ilc_zonca_sims_20204028_galaxy1_93-145-225-278_TT-EE_galmask%s_AZ.npy' %(s4_ilc_folder, galmask)

                '''
                #galmask = 1
                #totlf, totmf, tothf = 2, 12, 5
                extrastr_tubesmoved = 'lf%d-mf%d-hf%d' %(totlf, totmf, tothf)
                if (1):##also_te:
                    Nlfile_for_S4 = '%s/S4_ilc_zonca_sims_20204028_galaxy1_27-39-93-145-225-278_TT-EE-TE_%s_galmask%s_AZ.npy' %(s4_ilc_folder, extrastr_tubesmoved, galmask)
                    pspectra_to_use_arr = [['TT', 'EE'], ['TT','EE','TE']]            
                else:
                    Nlfile_for_S4 = '%s/S4_ilc_zonca_sims_20204028_galaxy1_27-39-93-145-225-278_TT-EE_%s_galmask%s_AZ.npy' %(s4_ilc_folder, extrastr_tubesmoved, galmask)
                    pspectra_to_use_arr = [['TT', 'EE'], ['TT','EE','TE']]
                '''


                s4_ilc_folder = 'data/DRAFT/results/20200601/s4like_mask/TT-EE-TE/'
                if totlf == 2 and totmf == 12 and tothf == 5:
                    s4_ilc_folder = '%s/baseline/' %(s4_ilc_folder)
                else:
                    s4_ilc_folder = '%s/tubes_mod/' %(s4_ilc_folder)
                extrastr_tubesmoved = 'lf%d-mf%d-hf%d' %(totlf, totmf, tothf)
                if include_gal:
                    Nlfile_for_S4 = '%s/S4_ilc_galaxy%s_27-39-93-145-225-278_TT-EE-TE_%s_galmask%s_AZ.npy' %(s4_ilc_folder, include_gal, extrastr_tubesmoved, galmask)
                else:
                    Nlfile_for_S4 = '%s/S4_ilc_galaxy%s_27-39-93-145-225-278_TT-EE-TE_%s_AZ.npy' %(s4_ilc_folder, include_gal, extrastr_tubesmoved)

                pspectra_to_use_arr = [['TT','EE','TE']]

                dic = np.load(Nlfile_for_S4, allow_pickle = 1, encoding = 'latin1').item()
                el_nl, cl_residual = dic['el'], dic['cl_residual']
                #force_fsky_val = dic['fsky_val']
                #force_fsky_dic = {0:0.53, 1:0.59, 2:0.65, 3:0.77}
                force_fsky_dic = {0:0.5298, 1:0.5885, 2:0.6504, 3:0.7732}
                force_fsky_val = force_fsky_dic[galmask]
                force_fsky_val = round(force_fsky_val, 4)

                if (0):
                    Nlfile_for_S4 = 'versions/20200529/S4_ilc_20203030.npy'
                    pspectra_to_use_arr = [['TT','EE','TE']]
                    dic = np.load(Nlfile_for_S4, allow_pickle = 1, encoding = 'latin1').item()
                    el_nl, cl_residual = dic['el'], dic['cl_residual']                
                    force_fsky_val = 0.6

                if 'T' in cl_residual:
                    Nl_TT, Nl_EE = cl_residual['T'], cl_residual['P']
                    Nl_TT = np.interp(els, el_nl, Nl_TT)
                    Nl_EE = np.interp(els, el_nl, Nl_EE)
                    Nl_TE = None
                else:
                    Nl_TT, Nl_EE = cl_residual['TT'], cl_residual['EE']
                    if 'TE' in cl_residual:
                        Nl_TE = cl_residual['TE']
                    else:
                        Nl_TE = None
                    el_nl = np.arange(len(Nl_TT))
                    Nl_TT = np.interp(els, el_nl, Nl_TT)
                    Nl_EE = np.interp(els, el_nl, Nl_EE)
                    if Nl_TE is not None:
                        Nl_TE = np.interp(els, el_nl, Nl_TE)

                Nldic = {}
                Nldic['TT'] = Nl_TT
                Nldic['EE'] = Nl_EE
                Nldic['TE'] = Nl_TE

                if (0):
                    ax = subplot(111, yscale = 'log')
                    plot(els, Nl_TT, ls = '-', color = 'k', label = r'TT'); 
                    plot(els, Nl_EE, ls = '-', color = 'orangered', label = r'EE'); 
                    if Nl_TE is not None:
                        plot(els, Nl_TE, ls = '-.', color = 'darkred', label = r'TE')
                    legend(loc = 1)
                    ylim(1e-9,1.);
                    show(); sys.exit()

                if add_lensing:
                    #nlfile = 'data/DRAFT/results/20200601/lensing_noise_curves/20200601/s4like_mask/TT-EE-TE/baseline/ilc/S4_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE_lf2-mf12-hf5_AZ_lmin100_lmax5000_lmaxtt3000.npy'
                    if extrastr_tubesmoved == 'lf2-mf12-hf5':
                        tmpstr = 'baseline/'
                    else:
                        tmpstr = 'tubes_mod/'
                    if include_gal == 0:
                        nlfile = 'data/DRAFT/results/20200601/lensing_noise_curves/20200601/s4like_mask/TT-EE-TE/%s/ilc/S4_ilc_galaxy0_27-39-93-145-225-278_TT-EE-TE_%s_AZ_lmin100_lmax5000_lmaxtt3000.npy' %(tmpstr, extrastr_tubesmoved)
                    else:
                        nlfile = 'data/DRAFT/results/20200601/lensing_noise_curves/20200601/s4like_mask/TT-EE-TE/%s/ilc/S4_ilc_galaxy1_27-39-93-145-225-278_TT-EE-TE_%s_galmask%s_AZ_lmin100_lmax5000_lmaxtt3000.npy' %(tmpstr, extrastr_tubesmoved, galmask)
                    print(nlfile)
                    Nl_PP, min_L, max_L = get_lensing_noise(nlfile, mod_lensing_N0_fac = mod_lensing_N0_fac, lensing_est = lensing_est)
                    Nldic['PP'] = Nl_PP
                    pspectra_to_use_arr[0].append('PP')

            if expname.lower().find('planck')>-1:
                beam = 5.
                deltaT = 30.
                deltaP = deltaT * np.sqrt(2.)
                fsky = 0.67
                #min_l_temp_arr_ori = [30]
                #max_l_temp_arr_ori = [2060]
                #min_l_pol_arr_ori = [2]
                #max_l_pol_arr_ori = [29]
                fsky_arr = [fsky]
                pspectra_to_use_arr = [['TT','EE','TE']]
                #Nlfile_arr = ['data/Planck/DR2/Nl_smica_halfmission_HDHM.npy']
                Nlfile = 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy'

                Nl_TT = np.load(Nlfile, allow_pickle=1)
                Nl_EE = np.sqrt(2.) * Nl_TT
                Nl_TE = np.copy(Nl_TT) * 0.

                Nldic = {}
                Nldic['TT'] = Nl_TT
                Nldic['EE'] = Nl_EE
                Nldic['TE'] = Nl_TE
                Nlfile = None

                if force_fsky_val != -1.:
                    fsky_arr = [force_fsky_val]

                min_l_temp_arr_ori = [2]
                max_l_temp_arr_ori = [2500]
                min_l_pol_arr_ori = [2]
                max_l_pol_arr_ori = [29]

            elif expname.lower().find('s4')>-1:
                beam = 1.4
                deltaT = 2.0
                deltaP = deltaT * np.sqrt(2.)
                #fsky = 0.4

                fsky_arr = np.arange( 0.1, 0.7, 0.1)
                #fsky_arr = [0.5]
                #fsky_arr = np.arange( 0.6, 0.8, 0.1)
                #min_l_temp_arr_ori = [30., 50., 100., 200., 250., 300.]
                max_l_temp_arr_ori = np.arange(2000., 6001., 500)
                #min_l_pol_arr_ori = [20., 50., 100., 200., 250., 300.]
                #min_l_pol_arr_ori = [30., 50., 100., 200., 250., 300.]
                max_l_pol_arr_ori = np.arange(2000., 6001., 500)

                max_l_temp_arr_ori = [2000., 3500., 5000.]#, 6000.]
                max_l_pol_arr_ori = [2000., 3500., 5000.]#, 6000.]


                if force_fsky_val != -1.:
                    fsky_arr = [force_fsky_val]

                min_l_temp_arr_ori = [30.]
                min_l_pol_arr_ori = [30.]

            extra_nl_dic_str = None
            if Nldic_for_inv_var is not None:
                Nldic, min_l_temp_arr_ori, min_l_pol_arr_ori, extra_nl_dic_str = add_planck_via_inv_variance(Nldic, Nldic_for_inv_var, include_planck_on_large_scales)

            if (0):##include_gal:
                pspectra_to_use_arr = [['TT'],['EE']]

            for pspectra_to_use in pspectra_to_use_arr:
                if pspectra_to_use == ['TT']:
                    min_l_temp_arr = np.copy( min_l_temp_arr_ori )
                    max_l_temp_arr = np.copy( max_l_temp_arr_ori )
                    min_l_pol_arr = [None]
                    max_l_pol_arr = [None]
                elif pspectra_to_use == ['EE'] or pspectra_to_use == ['TE'] or pspectra_to_use == ['EE','TE'] or pspectra_to_use == ['EE','TE']:
                    min_l_temp_arr = [None]
                    max_l_temp_arr = [None]
                    min_l_pol_arr = np.copy( min_l_pol_arr_ori )
                    max_l_pol_arr = np.copy( max_l_pol_arr_ori )
                else:
                    min_l_temp_arr = np.copy( min_l_temp_arr_ori )
                    max_l_temp_arr = np.copy( max_l_temp_arr_ori )
                    max_l_temp_arr = [3000.]
                    min_l_pol_arr = np.copy( min_l_pol_arr_ori )
                    max_l_pol_arr = np.copy( max_l_pol_arr_ori )

                max_l_temp_arr = [3000.] #3000.,]
                max_l_pol_arr = [5000.]

                for min_l_temp in min_l_temp_arr:
                    for max_l_temp in max_l_temp_arr:
                        for min_l_pol in min_l_pol_arr:
                            for max_l_pol in max_l_pol_arr:

                                if (max_l_temp is not None and max_l_pol is not None):
                                    if len(max_l_temp_arr) == 1:
                                        pass
                                    elif max_l_temp != max_l_pol: 
                                        continue

                                for fsky in fsky_arr:
                                #for Nlfile in Nlfile_arr:

                                    pspectra_to_use_str = '-'.join(pspectra_to_use)
                                    '''
                                    if Nlfile is not None:
                                        Nlfile_str = Nlfile.split('/')[-1].replace('_','#')
                                    else:
                                        Nlfile_str = None
                                    otherkeys = '%sminltemp_%smaxltemp_%sminlpol_%smaxlpol_%s_%sNlfile' %(min_l_temp, max_l_temp, min_l_pol, max_l_pol, pspectra_to_use_str, Nlfile_str)
                                    '''
                                    otherkeys = '%sfsky_%sminltemp_%smaxltemp_%sminlpol_%smaxlpol_%s' %(fsky, min_l_temp, max_l_temp, min_l_pol, max_l_pol, pspectra_to_use_str)
                                    keyname = '%s_%s_galaxy%s' %(expname, otherkeys, include_gal)
                                    if galmask is not None:
                                        keyname = '%s_galmask%s' %(keyname, galmask)
                                    if extrastr_tubesmoved is not None:
                                        keyname = '%s_%s' %(keyname, extrastr_tubesmoved)
                                    if extra_nl_dic_str is not None:
                                        keyname = '%s_%s' %(keyname, extra_nl_dic_str)
                                    if add_lensing:
                                        keyname = '%s_pluslensing' %(keyname)

                                    print(keyname)

                                    exp_dic[keyname] = {\
                                    'beam': beam, 'fsky':fsky, 'deltaT': deltaT, 'deltaP': deltaP, 'Nlfile': Nlfile, 'Nldic': Nldic, \
                                    'min_l_temp': min_l_temp, 'min_l_pol': min_l_pol, 'max_l_temp': max_l_temp, 'max_l_pol': max_l_pol,\
                                    'pspectra_to_use':pspectra_to_use,
                                    }
        expnamearr = sorted( exp_dic.keys() )
    elif expname.find('spt')>-1:
        beam = 1.2
        if expname == 'sptpol':
            deltaT = 6.
            fsky = (500./41253.)
        if expname == 'sptpollowl':
            deltaT = 6.
            fsky = (500./41253.)
        if expname == 'sptpolhighl':
            deltaT = 6.
            fsky = (500./41253.)
        elif expname == 'sptpolsummer':
            deltaT = 27.
            fsky = (2700./41253.)
        elif expname == 'sptpol3gsummer':
            deltaT = 12.
            fsky = (600./41253.)
        elif expname == 'sptsz':
            deltaT = 18.
            fsky = (2500./41253.)
        elif expname == 'spt3gsummer':
            deltaT = 15.
            #fsky = (1900./41253.)
            fsky = (1500./41253.)
        elif expname == 'spt3gy12':
            deltaT = 9.
            fsky = (1500./41253.)

        #fsky = 1.
        #print('\n\nsetting all fsky = 1\n')

        min_l_temp, max_l_temp = 100, 3000
        min_l_pol, max_l_pol = 100, 8000

        if expname == 'sptpollowl':
            min_l_temp, max_l_temp = 50, 1000
            min_l_pol, max_l_pol = 50, 1000
        elif expname == 'sptpolhighl':
            min_l_temp, max_l_temp = 1000, 8000
            min_l_pol, max_l_pol = 1000, 8000


        deltaP = None
        ##pspectra_to_use = ['TT', 'EE', 'TE']
        if expname.find('sptpol')>-1:
            pspectra_to_use = ['EE', 'TE']

        Nldic = None

        Nlfile_spt = 'data/DRAFT/results/spt/%s_ilc_90-150-220_TT-EE.npy' %(expname)
        Nlfile_spt = Nlfile_spt.replace('lowl','').replace('highl','')

        if which_spectra != 'lens_potential':
            Nldic = get_Nldic(Nlfile_spt, els)
            if Nldic_for_inv_var is not None:
                Nldic, min_l_temp_arr_ori, min_l_pol_arr_ori, extra_nl_dic_str = add_planck_via_inv_variance(Nldic, Nldic_for_inv_var, include_planck_on_large_scales)
        else:
            Nldic = {}

        if add_lensing or which_spectra == 'lens_potential':

            if Nlfile is None:
                sptlensingnoisecurvesfolder = 'data/spt_summer_fields/lensing_noise_curves/ilc/'
                lmin, lmax = 100, 3000
                Nlfile = '%s/%s_lmin%s_lmax%s.npy' %(sptlensingnoisecurvesfolder, expname, lmin, lmax)
                Nlfile = Nlfile.replace('lowl','').replace('highl','')
            if which_spectra == 'lens_potential':
                pspectra_to_use = ['PP']
            else:
                pspectra_to_use.append( 'PP' )

            Nl_PP, min_L, max_L = get_lensing_noise(Nlfile, mod_lensing_N0_fac = mod_lensing_N0_fac, lensing_est = lensing_est)

            min_l_temp, max_l_temp = min_L, max_L
            min_l_pol, max_l_pol = min_L, max_L

            Nldic['PP'] = Nl_PP

        spt_dic = \
        {'beam':beam, 'deltaT':deltaT, 'deltaP': deltaP, 'fsky':fsky, 'Nlfile': Nlfile, 'Nldic': Nldic,\
        'pspectra_to_use': pspectra_to_use,\
        'min_l_temp': min_l_temp, 'min_l_pol': min_l_pol, 'max_l_temp': max_l_temp, 'max_l_pol': max_l_pol, 'color': 'black', 'label': r'%s' %(expname),\
        }
        exp_dic = {}
        exp_dic = {expname:spt_dic}
        expnamearr = [expname]        
    elif expname == 's4like_simple':
        beam = 1.4
        deltaT = 2.
        deltaP = np.sqrt(2.) * deltaT
        fsky = 0.4
        Nlfile = None
        Nldic = None
        min_l_temp, max_l_temp = 30, 5000
        min_l_pol, max_l_pol = 30, 5000
        #pspectra_to_use = [['TT'], ['TT', 'EE', 'TE']]
        pspectra_to_use = [['TT', 'EE', 'TE']]

        s4like_simple_dic = \
        {'beam':beam, 'deltaT':deltaT, 'deltaP': deltaP, 'fsky':fsky, 'Nlfile': Nlfile, 'Nldic': Nldic,\
        'pspectra_to_use': pspectra_to_use,\
        'min_l_temp': min_l_temp, 'min_l_pol': min_l_pol, 'max_l_temp': max_l_temp, 'max_l_pol': max_l_pol, 'color': 'black', 'label': r'S4-like simple',\
        }
        exp_dic = {}
        exp_dic = {expname:s4like_simple_dic}
        expnamearr = [expname]        
    elif expname == 'planck_s4':
        cmap = py_ini.get_planck_cmap()
        cmap = cm.jet
        coloarr = [cmap(int(d)) for d in np.linspace(0,255,5)]
        planck_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': None,\
        #'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy', 
        #'pspectra_to_use': '[TT, EE, TE]',\
        'pspectra_to_use': ['TT', 'EE', 'TE'],\
        'min_l_temp': 30, 'min_l_pol': 2, 'max_l_temp': 2000, 'max_l_pol': 29, 'color': 'black', 'label': r'{\it Planck} TT$+$lowP',\
        }
        planck_smica_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy',\
        'pspectra_to_use': ['TT'],\
        'min_l_temp': 30, 'min_l_pol': 2, 'max_l_temp': 2500, 'max_l_pol': 29, 'color': coloarr[0], 'label': r'{\it Planck} \texttt{SMICA} TT',\
        }
        planck2_smica_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy',\
        'pspectra_to_use': ['TT', 'EE'],\
        #'min_l_temp': 30, 'min_l_pol': 2, 'max_l_temp': 2500, 'max_l_pol': 29, 'color': 'goldenrod', 'label': r'{\it Planck} \texttt{SMICA} TT+EE',\
        'min_l_temp': 30, 'min_l_pol': 30, 'max_l_temp': 2500, 'max_l_pol': 100, 'color': coloarr[1], 'label': r'{\it Planck} \texttt{SMICA} TT+EE (T:30 - 2500, P:30-100)',\
        }
        planck3_smica_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy',\
        'pspectra_to_use': ['TT', 'EE', 'TE'],\
        'min_l_temp': 30, 'min_l_pol': 30, 'max_l_temp': 2500, 'max_l_pol': 100, 'color': coloarr[2], 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE (T:30 - 2500, P:30-100)',\
        #'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 29, 'max_l_pol': 29, 'color': 'forestgreen', 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE',\
        }
        planck4_smica_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy',\
        'pspectra_to_use': ['TT', 'EE', 'TE'],\
        'min_l_temp': 30, 'min_l_pol': 10, 'max_l_temp': 2500, 'max_l_pol': 29, 'color': coloarr[3], 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE  (T:30 - 2500, P:10-29)',\
        }
        planck5_smica_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy',\
        'pspectra_to_use': ['TT', 'EE', 'TE'],\
        'min_l_temp': 30, 'min_l_pol': 2, 'max_l_temp': 2500, 'max_l_pol': 29, 'color': coloarr[4], 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE  (T:30 - 2500, P:2-29)',\
        #'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 29, 'max_l_pol': 29, 'color': 'darkred', 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE',\
        }
        planck6_smica_dic = \
        {'beam':5., 'deltaT':30., 'deltaP': None, 'fsky':.67, 'Nlfile': 'data/Planck/DR2/Nl_smica_halfmission_HDHM.npy',\
        'pspectra_to_use': ['TT', 'EE', 'TE'],\
        'min_l_temp': 30, 'min_l_pol': 30, 'max_l_temp': 2500, 'max_l_pol': 100, 'color': 'k', 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE (T:30 - 2500, P:30-100) + $\tau$ prior',\
        'prior_dic': {'tau':0.01},
        #'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 29, 'max_l_pol': 29, 'color': 'forestgreen', 'label': r'{\it Planck} \texttt{SMICA} TT+EE+TE',\
        }
        cmbs4_TT_dic  = {'beam':1.4, 'deltaT':2.0, 'deltaP': None, 'fsky':.4, 'Nlfile': None, 'pspectra_to_use': ['TT'], 'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 2000, 'max_l_pol': 3000,\
        'color': 'orangered', 'label': r'CMB-S4 TT $\ell_{\rm max} = 2000$',\
        }
        cmbs4_all_dic  = {'beam':1.4, 'deltaT':2.0, 'deltaP': None, 'fsky':.4, 'Nlfile': None, 'pspectra_to_use': ['TT', 'EE', 'TE'], 'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 2000, 'max_l_pol': 3000,\
        'color': 'goldenrod', 'label': r'CMB-S4 TTEETE $\ell_{\rm max}$ = (T:2000, P:3000)',\
        }
        cmbs4_EETE_lmax3000pol_dic  = {'beam':1.4, 'deltaT':2.0, 'deltaP': None, 'fsky':.4, 'Nlfile': None, 'pspectra_to_use': ['EE', 'TE'], 'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 2000, 'max_l_pol': 3000,\
        'color': 'forestgreen', 'label': r'CMB-S4 EETE $\ell_{\rm max} =3000$',\
        }
        cmbs4_EETE_lmax6000pol_dic  = {'beam':1.4, 'deltaT':2.0, 'deltaP': None, 'fsky':.4, 'Nlfile': None, 'pspectra_to_use': ['EE', 'TE'], 'min_l_temp': 2, 'min_l_pol': 2, 'max_l_temp': 2000, 'max_l_pol': 6000,\
        'color': 'darkred', 'label': r'CMB-S4 EETE $\ell_{\rm max} =6000$',\
        }

        exp_dic = {}
        exp_dic = {'planck':planck_smica_dic, 'cmbs4all': cmbs4_all_dic, 'cmbs4TT': cmbs4_TT_dic, 'cmbs4EETElmax3000pol': cmbs4_EETE_lmax3000pol_dic}#, 'cmbs4EETElmax6000pol':cmbs4_EETE_lmax6000pol_dic}
        #exp_dic = {'cmbs4_EETE_lmax3000pol': cmbs4_EETE_lmax3000pol_dic, 'cmbs4_EETE_lmax6000pol':cmbs4_EETE_lmax6000pol_dic}
        expnamearr = ['planck', 'cmbs4TT', 'cmbs4all']

        exp_dic = {'planck':planck_smica_dic, 'planck2':planck2_smica_dic, 'planck3':planck3_smica_dic, 'planck4':planck4_smica_dic, 'planck5':planck5_smica_dic, 'planck6':planck6_smica_dic}
        expnamearr = ['planck', 'planck2', 'planck3', 'planck4', 'planck5', 'planck6']
        #expnamearr = ['planck2', 'planck3']

    return exp_dic, expnamearr
