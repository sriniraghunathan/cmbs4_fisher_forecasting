def get_latex_param_str(param):
    params_str_dic= {\
    'norm_YszM': r'${\rm log}(Y_{\ast})$', 'alpha_YszM': r'$\alpha_{_{Y}}$',\
    'beta_YszM': r'$\beta_{_{Y}}$', 'gamma_YszM': r'$\gamma_{_{Y}}$', \
    'alpha': r'$\eta_{\rm v}$', 'sigma_8': r'$\sigma_{\rm 8}$', \
    'one_minus_hse_bias': r'$1-b_{\rm SZ}$', 'omega_m': r'$\Omega_{\rm m}$', \
    'h0':r'$h$', 'm_nu':r'$\sum m_{\nu}$', 'ombh2': r'$\Omega_{b}h^{2}$', 'omch2': r'$\Omega_{c}h^{2}$', 'w0': r'$w_{0}$', 'wa': r'$w_{a}$', \
    'tau': r'$\tau_{\rm re}$', 'As': r'$A_{\rm s}$', 'ns': r'$n_{\rm s}$', 'neff': r'N$_{\rm eff}$', \
    'mnu': r'$\sum m_{\nu}$', 'thetastar': r'$\theta_{\ast}$', \
    }

    return params_str_dic[param]

################################################

import numpy as np, glob
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = r'\boldmath'
from matplotlib.patches import Ellipse
rcParams['font.family'] = 'serif'
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

rounding = 0
if not rounding:
    fd = 'no_rounding/'
else:
    fd = 'rounding_0p4/'
flist = glob.glob('%s/*.txt' %(fd))

fsval = 14
colorarr = ['black', 'orangered']
for fcntr, f in enumerate( flist ):
    pname, fsky, fwhm = f.split('/')[-1].replace('.txt','').split('_')
    fsky, fwhm = float(fsky.replace('fsky','')), float(fwhm.replace('fwhm','').replace('am',''))
    print(pname, fsky, fwhm)
    pstr = get_latex_param_str(pname)
    pstr = '(%s)' %(pstr)
    noise, sigma = np.loadtxt(f, unpack = 1)
    if pname == 'ns': 
        sigma *= 10.
        pstr = r'%s $\times$10' %(pstr)

    #labval = r'$f_{\rm sky}$ = %.2f; Beam = $%.2f^{\prime}$' %(fsky, fwhm)
    labval = r'$\sigma$%s: f$_{\rm sky}$ = %.2f; Beam = $%.1f^{\prime}$' %(pstr, fsky, fwhm)
    plot(noise, sigma, label = labval, color = colorarr[fcntr])
    legend(loc = 2, fancybox = True)#, fontsize = fsval-2)

    if (1):
        ylabel(r'$\sigma(\theta)$', fontsize = fsval)
        xlabel(r'Noise $\Delta_{T}$ [$\mu$K-arcmin]', fontsize = fsval)
        ylim(0.001, 0.1)
    else:
        if pname == 'ns':
            ylim(0.001, 0.010)
        elif pname == 'neff':
            ylim(0.01, 0.1)
        #ylabel(r'$\sigma$(%s)' %(pstr))

    xlim(0., 10.)
    grid(True, which='major', axis = 'both', lw = 0.5, alpha = 0.1)
    grid(True, which='minor', axis = 'both', lw = 0.1, alpha = 0.1)
    titstr = r'Model: 6-LCDM + %s; Rounding = %s' %(get_latex_param_str('neff'), rounding)
    title(titstr, fontsize = fsval)
plname = 'parameter_constraints'
if not rounding:
    plname = '%s_norounding.png' %(plname)
else:
    plname = '%s_rounding0p4.png' %(plname)
savefig(plname, dpi = 200.)
#show()