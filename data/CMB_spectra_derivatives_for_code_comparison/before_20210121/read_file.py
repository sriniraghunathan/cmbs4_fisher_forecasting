import numpy as np, sys

"""
#reading CMB spectra from the text file
"""
fname = 'cmb_spectra_lensed_Srini.txt' ## 'cmb_spectra_unlensed_Srini.txt'
ell, tt, ee, te, bb = np.loadtxt(fname, unpack = 1)
#print(ell, tt, ee, te, bb)


"""
#reading CMB spectra derivatives from the python numpy file (dictionary)
"""
fname = 'cmb_spectra_derivs_lensed_Srini.npy' ##'cmb_spectra_derivs_unlensed_Srini.npy'
cl_deriv_dic = np.load(fname, allow_pickle = 1).item()
#it this does not work then try with enconding = 'latin1' option
#cl_deriv_dic = np.load(fname, allow_pickle = 1, enconding = 'latin1').item()
print( '\n\tCosmological parameters are: %s (ignore \ell)\n' %(sorted(cl_deriv_dic.keys())) )
ell = cl_deriv_dic['ell'] #\els 

#looping over all parameters
print('\n\tReading derivatives of Cl w.r.t cosmological parameters\n')
for param_name in sorted( cl_deriv_dic ):
    if param_name == 'ell': continue
    print( '\n\t\t%s' %(param_name) )
    cl_der_tt, cl_der_ee, cl_der_te = cl_deriv_dic[param_name]
print('\n')
sys.exit()
