from astropy.io import fits
import sys
from scipy.optimize import minimize
from math import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})


def read_and_convert_science(filename, fact):
	hdul = fits.open(filename)
	data = hdul[1].data

	xstart = 1264
	xend   = 15112
	ystart = 376
	yend   = 660

	#val1 = np.zeros ((12, 2312)) 
	#n_rec = np.zeros ((12, 2312))
	val1 = np.zeros ((12, 289*fact)) 
	n_rec = np.zeros ((12, 289*fact))

	## Reading corrtag file and preparing binned image 
	for i in range (0, len(data['XCORR'][:])):
		if (data['XCORR'][i] > xstart) and (data['XCORR'][i] < xend-6):
			if (data['YCORR'][i] > ystart) and (data['YCORR'][i] < yend) and (data['PHA'][i] < 28) and (data['PHA'][i] > 3):
				jj = int ((data['XCORR'][i] - xstart) / 6.0 / 8.0 * fact)
				kk = int ((data['YCORR'][i] - ystart) / 12.0 / 2.0)
				if ((int(data['PHA'][i]) > 3) and (int(data['PHA'][i]) <28)): 
					val1 [kk, jj]  = val1[kk, jj] + data['EPSILON'][i]	
					n_rec [kk, jj] = n_rec[kk, jj] + 1.0

	return [val1, n_rec]	


fact = 2


sciences = glob.glob('/astro/sveash/cos_dark/corr_code/testing_M83-1/ldaw01a7q*corrtag_a.fits')
rootdir = 'corr_code'

for s in sciences:
        ## get rootname of file
        rootn = s.split('-1/')[1].split('_corr')[0]
        
        ## Read predicted noise level
        noise = np.load(f'{rootdir}{rootn}_noise_complete.npy')
        
        ## Read science exposure
        signal, n_rec = read_and_convert_science (s, fact)
        print('Size signal: ', len(signal[0]))

        hdul = fits.open(s) 
        data = hdul[1].data

        xstart = 1264
        xend   = 15112
        ystart = 376
        yend   = 660

        wwf = open ('tst.txt', 'w')

        logic = np.zeros((signal.shape[0], signal.shape[1]))
        print(len(logic), np.shape(logic))

        #sys.exit(0)

        for i in range (0, len(data['XCORR'][:])):
                if (data['XCORR'][i] > xstart) and (data['XCORR'][i] < xend-6):
                        if (data['YCORR'][i] > ystart) and (data['YCORR'][i] < yend) and (data['PHA'][i] < 28) and (data['PHA'][i] > 3):
                                jj = int ((data['XCORR'][i] - xstart) / 6.0 / 8.0 * fact)
                                kk = int ((data['YCORR'][i] - ystart) / 12.0 / 2.0)

                                #ss = int(data['PHA'][i] / 2.0)
                                #print('jj, kk ', jj, kk )

                                delta_eps = noise[kk][int(jj / fact)] / float(fact)
                                logic [kk][jj] = 1
                                delta_eps_ind = delta_eps / float(n_rec[kk, jj])

                                wwf.write(str(data['EPSILON'][i]) + '\t' + str(delta_eps) + '\t' + str(delta_eps_ind) + '\n')
                                data['EPSILON'][i] = data['EPSILON'][i] - delta_eps_ind

        #hdul.writeto('lcm802inq_corrtag_a_after_correction.fits')

        plt.imshow(logic, aspect="auto", origin="lower")
        plt.show()

        print(logic.shape[0])
        #sys.exit(0)

        nmb_zero = logic.shape[0] * logic.shape[1] - np.sum(logic)

        if nmb_zero > 1:
                new_records = data[len(data)-int(nmb_zero):].copy()
        print('Null records: ', nmb_zero, ' copied: ', len(new_records))

        #sys.exit(0)

        cnt = 0
        for i in range (0, logic.shape[0]):
                for j in range (0, logic.shape[1]):
                        if logic[i,j] == 0:
                                y_coord = 12.0 * 2.0 * i       + ystart
                                x_coord = 6.0 * 8.0 * j / fact + xstart

                                new_records['RAWX'][cnt] = x_coord
                                new_records['RAWY'][cnt] = y_coord

                                new_records['XCORR'][cnt] = x_coord
                                new_records['YCORR'][cnt] = y_coord

                                new_records['XFULL'][cnt] = x_coord
                                new_records['YFULL'][cnt] = y_coord

                                new_records['EPSILON'][cnt] = -noise[i][int(j / fact)] / fact
                                new_records['PHA'][cnt] = 10

                                cnt = cnt + 1

        data = np.concatenate ([data, new_records])	
        hdul[1].data = data
        hdul.writeto(f'corrected_{rootn}_corrtag_a.fits', overwrite=True)
                                #val1[kk, jj, ss] = val1[kk, jj, ss] + data['EPSILON'][i]

        #plt.plot (cmb)
        #plt.plot (signal)
        #plt.show()

        #noise_0 = np.zeros(len(noise))

        #cmb, clean = combined (list_of_mu, list_of_gaussian, noise_0)

        #plt.plot (cmb)
        #plt.plot (signal)
        #plt.show()

        #np.save ('clean', cmb)

