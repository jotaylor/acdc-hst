from scipy.optimize import minimize
from astropy.io import fits
import numpy as np
from math import *
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})
import glob as glob

def read_dark (filename):
	dark = np.load(filename)

	return dark

def linear_combination (darks, coeffs):
	if len(darks) != len(coeffs):
		print ('Size of files is wrong')
		sys.exit(0)

	for i in range (0, len(coeffs)):
		if i == 0:
			final = np.copy(darks[0]) * coeffs [0]
		else:
			final = final + darks[i] * coeffs [i]
		#print ('i -- ', i)
		
	return final

def C_stat (combined_superdark, science_exposure, excluded_rows):
	Csum = 0.0
	for i in range (0, combined_superdark.shape[0]):
		for j in range (0, combined_superdark.shape[1]):
			if science_exposure[i][j] > 0 and (i not in excluded_rows) and combined_superdark[i][j] > 0:
				Csum = Csum + 2.0*(combined_superdark[i][j] - science_exposure[i][j] + science_exposure[i][j] * (log(science_exposure[i][j]) - log(combined_superdark[i][j])))	
				#Csum = Csum + 2.0 * (combined_superdark[i][j] - int(science_exposure[i][j]) * log(combined_superdark[i][j]) + log(factorial(int(science_exposure[i][j]))) )	
				#if science_exposure[i][j] > 10:
				#	print ('Warning: ', i, j, science_exposure[i][j])
			elif (i not in excluded_rows):
				Csum = Csum + abs(combined_superdark[i][j]) * 2
				#print (i, j, Csum)

	return Csum

def fun_opt (coeff, darks, science_exposure, excluded_rows):
	combined_superdark = linear_combination (darks, coeff)
	Cval = C_stat (combined_superdark, science_exposure, excluded_rows)
	print (coeff, Cval)
	
	return Cval


def read_and_convert_science (filename):
	hdul = fits.open(filename)
	data = hdul[1].data

	xstart = 1264
	xend   = 15112
	ystart = 376
	yend   = 660
	val1 = np.zeros ((12, 289)) 

	## Reading corrtag file and preparing binned image 
	for i in range (0, len(data['XCORR'][:])):
		if (data['XCORR'][i] > xstart) and (data['XCORR'][i] < xend-6):
			if (data['YCORR'][i] > ystart) and (data['YCORR'][i] < yend) and (data['PHA'][i] < 28) and (data['PHA'][i] > 3):
				jj = int ((data['XCORR'][i] - xstart) / 6.0 / 8.0)
				kk = int ((data['YCORR'][i] - ystart) / 12.0 / 2.0)
				if ((int(data['PHA'][i]) > 3) and (int(data['PHA'][i]) <28)): 
					val1 [kk, jj] = val1[kk, jj] + data['EPSILON'][i]	

	return val1	

exp_superdark = 259358.81600000005

sciences = glob.glob('/Users/sveash/Dropbox/COS_dark/low_snr/M83-1/*_corrtag_a.fits')
rootdir = '/Users/sveash/Dropbox/COS_dark/low_snr/M83-1/dark_corr/v2/'

for s in sciences:
    science = read_and_convert_science (s)  ## file with science exposure
    rootn = s.split('-1/')[1].split('_corr')[0]
    dark0 = read_dark ('sp_FUVA_167_300days.npy') ## quiescent state dark
    dark1 = read_dark ('superdark_ba_max.npy')    ## any number of peculiar superdarks
    plt.imshow(dark0, aspect='auto')
    plt.show()

    exp_obs = fits.getheader(s,1)['EXPTIME']
    tau_exposure = exp_superdark / exp_obs

    darks  = [dark0, dark1]
    coeffs = [0.5 / tau_exposure, 0.5 / tau_exposure]
    combined_superdark = linear_combination (darks, coeffs)

    plt.title(f'{rootn}')
    plt.imshow(science, aspect='auto')
    plt.show()

    bb = plt.imshow(combined_superdark, aspect='auto')
    plt.colorbar(bb)
    plt.show()

    excluded_rows = [2,3,4,7,8] #May need to adjust these depending on LP

    for i in range (0, science.shape[0]):
            if i not in excluded_rows:
                    plt.plot (science[i], label=i)
    plt.legend()
    plt.show()

    val_C = C_stat (combined_superdark, science, excluded_rows)
    print ('First ', val_C)	

    #coeffs = [0.5 / tau_exposure, 0.5]
    #val_C = fun_opt (coeffs, darks, science_exposure, excluded_rows)
    #print ('Second ', val_C)

    x0 = [0.007, 0.005]
    res = minimize(fun_opt, x0, method='Nelder-Mead', tol=1e-6, args=(darks, science, excluded_rows))
    print (res.x)

    combined_superdark = linear_combination (darks, res.x)
    bb = plt.imshow(combined_superdark, aspect='auto')
    plt.colorbar(bb)
    plt.show()

    NM = 9 ## row with science data
    plt.plot (science[NM], label='Outside of science extraction')
    plt.plot (combined_superdark[NM], label='Predicted dark level')
    plt.xlabel('x')
    plt.ylabel('Number of photons')
    plt.legend()
    plt.savefig(f'{rootdir}{rootn}_predicted_dark.pdf')

    plt.show()
    np.save (f'{rootdir}{rootn}_noise', combined_superdark[NM])
    np.save (f'{rootdir}{rootn}_signal', science[NM])

    print (np.mean(science[NM][180:]))
    print (np.mean(combined_superdark[NM][180:]))

    np.save (f'{rootdir}{rootn}_noise_complete', combined_superdark)

