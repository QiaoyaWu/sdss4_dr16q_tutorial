#!/usr/bin/python3

import os,sys,time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

if __name__ == "__main__":

    generate_global_fits_para = True
    generate_local_fits_para = True
    if generate_global_fits_para:
        ## generate fitting parameters for global fit
                               # lambda cpname min/maxwav lname          ngauss inisig  minsig      maxsig      voff  vindex windex findex  fvalue
        newdata = np.rec.array([(6564.61,'Ha',6400.,6800.,'Halpha_br',      3,  5e-3,   0.002,      0.05,       0.015,  0,  0,  0,  0.01),\
                                (6564.61,'Ha',6400.,6800.,'Halpha_na',      1,  1e-3,   2.3e-4,     0.002,      5e-3,   1,  1,  0,  0.002),\
                                (6549.85,'Ha',6400.,6800.,'NII6549',        1,  1e-3,   2.3e-4,     0.002,      5e-3,   1,  1,  1,  0.001),\
                                (6585.28,'Ha',6400.,6800.,'NII6585',        1,  1e-3,   2.3e-4,     0.002,      5e-3,   1,  1,  1,  0.003),\
                                (6718.29,'Ha',6400.,6800.,'SII6718',        1,  1e-3,   2.3e-4,     0.002,      5e-3,   1,  1,  2,  0.001),\
                                (6732.67,'Ha',6400.,6800.,'SII6732',        1,  1e-3,   2.3e-4,     0.002,      5e-3,   1,  1,  2,  0.001),\

                                (4862.68,'Hb',4640.,5100.,'Hbeta_br',       3,  5e-3,   0.002,      0.05,       0.01,   0,  0,  0,  0.01),\
                                (4862.68,'Hb',4640.,5100.,'Hbeta_na',       1,  1e-3,   2.3e-4,     0.002,      0.01,   1,  1,  0,  0.002),\
                                (4960.30,'Hb',4640.,5100.,'OIII4959c',      1,  1e-3,   2.3e-4,     0.002,      0.01,   1,  1,  0,  0.002),\
                                (5008.24,'Hb',4640.,5100.,'OIII5007c',      1,  1e-3,   2.3e-4,     0.002,      0.01,   1,  1,  0,  0.004),\
                                (4960.30,'Hb',4640.,5100.,'OIII4959w',      1,  3e-3,   2.3e-4,     0.004,      0.01,   2,  2,  0,  0.001),\
                                (5008.24,'Hb',4640.,5100.,'OIII5007w',      1,  3e-3,   2.3e-4,     0.004,      0.01,   2,  2,  0,  0.002),\
                                (4687.02,'Hb',4640.,5100.,'HeII4687_br',    1,  5e-3,   0.002,      0.05,       0.005,  0,  0,  0,  0.001),\
                                (4687.02,'Hb',4640.,5100.,'HeII4687_na',    1,  1e-3,   2.3e-4,     0.002,      0.005,  1,  1,  0,  0.001),\

                                #(3934.78,'CaII',3900.,3960.,'CaII3934',     2,  1e-3,   3.333e-4,   0.002,     0.01,   99, 0,  0,  -0.001),\

                                #(3728.48,'OII',3650.,3800.,'OII3728',       1,  1e-3,   3.333e-4,   0.002,     0.01,   1,  1,  0,  0.001),\

                                #(3426.84,'NeV',3380.,3480.,'NeV3426',       1,  1e-3,   3.333e-4,   0.002,     0.01,   0,  0,  0,  0.001),\
                                #(3426.84,'NeV',3380.,3480.,'NeV3426_br',    1,  5e-3,   0.0025,     0.02,       0.01,   0,  0,  0,  0.001),\

                                (2798.75,'MgII',2700.,2900.,'MgII_br',      2,  5e-3,   0.002,      0.05,       0.0015, 0,  0,  0,  0.05),\
                                (2798.75,'MgII',2700.,2900.,'MgII_na',      1,  1e-3,   5e-4,       0.002,      0.01,   1,  1,  0,  0.002),\

                                (1908.73,'CIII',1700.,1970.,'CIII_br1',     1,  5e-3,   0.002,      0.05,       0.015,  3, 0,  0,  0.01),\
                                (1908.73,'CIII',1700.,1970.,'CIII_br2',     1,  5e-3,   0.002,      0.05,       0.015,  3, 0,  0,  0.01),\
                                #(1908.73,'CIII',1700.,1970.,'CIII_na',      1,  1e-3,   5e-4,       0.002,      0.01,   4,  4,  0,  0.002),\
                                (1892.03,'CIII',1700.,1970.,'SiIII1892',    1,  2e-3,   0.001,      0.015,      0.003,  1,  1,  0,  0.005),\
                                (1857.40,'CIII',1700.,1970.,'AlIII1857',    1,  2e-3,   0.001,      0.015,      0.003,  1,  1,  0,  0.005),\
                                (1816.98,'CIII',1700.,1970.,'SiII1816',     1,  2e-3,   0.001,      0.015,      0.01,   2,  2,  0,  0.0002),\
                                (1750.26,'CIII',1700.,1970.,'NIII1750',     1,  2e-3,   0.001,      0.015,      0.01,   2,  2,  0,  0.001),\
                                (1718.55,'CIII',1700.,1900.,'NIV1718',      1,  2e-3,   0.001,      0.015,      0.01,   2,  2,  0,  0.001),\

                                (1549.06,'CIV',1500.,1700.,'CIV_br',        3,  5e-3,   0.001,      0.05,       0.015,  0,  0,  0,  0.05),\
                                #(1549.06,'CIV',1500.,1700.,'CIV_na',        1,  1e-3,   5e-4,       0.002,      0.01,   1,  1,  0,  0.002),\
                                (1640.42,'CIV',1500.,1700.,'HeII1640',      1,  1e-3,   5e-4,       0.002,      0.008,  1,  1,  0,  0.002),\
                                (1663.48,'CIV',1500.,1700.,'OIII1663',      1,  1e-3,   5e-4,       0.002,      0.008,  1,  1,  0,  0.002),\
                                (1640.42,'CIV',1500.,1700.,'HeII1640_br',   1,  5e-3,   0.0025,     0.02,       0.008,  2,  2,  0,  0.002),\
                                (1663.48,'CIV',1500.,1700.,'OIII1663_br',   1,  5e-3,   0.0025,     0.02,       0.008,  2,  2,  0,  0.002),\

                                (1402.06,'SiIV',1290.,1450.,'SiIV_OIV1',    1,  5e-3,   0.002,      0.05,       0.015,  1,  1,  0,  0.05),\
                                (1396.76,'SiIV',1290.,1450.,'SiIV_OIV2',    1,  5e-3,   0.002,      0.05,       0.015,  1,  1,  0,  0.05),\
                                (1335.30,'SiIV',1290.,1450.,'CII1335',      1,  2e-3,   0.001,      0.015,      0.01,   2,  2,  0,  0.001),\
                                (1304.35,'SiIV',1290.,1450.,'OI1304',       1,  2e-3,   0.001,      0.015,      0.01,   2,  2,  0,  0.001),\

                                (1215.67,'Lya',1150.,1290.,'Lya_br',        3,  5e-3,   0.002,      0.05,       0.02,   0,  0,  0,  0.05),\
                                #(1215.67,'Lya',1150.,1290.,'Lya_na',        1,  1e-3,   5e-4,       0.002,      0.01,   0,  0,  0,  0.002),\
                                (1240.14,'Lya',1150.,1290.,'NV1240',        1,  2e-3,   0.001,      0.01,       0.005,  0,  0,  0,  0.002)\
                                ],\
                             formats='float32,a20,float32,float32,a20,float32,float32,float32,float32,\
                             float32,float32,float32,float32,float32',\
                             names='lambda,compname,minwav,maxwav,linename,ngauss,inisig,minsig,maxsig,voff,vindex,windex,findex,fvalue')
        #------header-----------------
        hdr = fits.Header()
        hdr['lambda'] = 'Vacuum Wavelength in Ang'
        hdr['minwav'] = 'Lower complex fitting wavelength range'
        hdr['maxwav'] = 'Upper complex fitting wavelength range'
        hdr['ngauss'] = 'Number of Gaussians for the line'
        hdr['inisig'] = 'Initial guess of linesigma [in lnlambda]'
        hdr['minsig'] = 'Lower range of line sigma [lnlambda]'
        hdr['maxsig'] = 'Upper range of line sigma [lnlambda]'
        hdr['voff  '] = 'Limits on velocity offset from the central wavelength [lnlambda]'
        hdr['vindex'] = 'Entries w/ same NONZERO vindex constrained to have same velocity'
        hdr['windex'] = 'Entries w/ same NONZERO windex constrained to have same width'
        hdr['findex'] = 'Entries w/ same NONZERO findex have constrained flux ratios'
        hdr['fvalue'] = 'Relative scale factor for entries w/ same findex'
        #------save line info-----------
        hdu = fits.BinTableHDU(data=newdata,header=hdr,name='data')
        hdu.writeto('./qsopar.fits',overwrite=True)

    if generate_local_fits_para:
    ## generate fitting parameters for local fit
                                # lambda cpname min/maxwav lname          ngauss inisig  minsig      maxsig      voff  vindex windex findex  fvalue
        newdata = np.rec.array([(3934.78,'CaII',3900.,3960.,'CaII3934_1',     1,  1e-3,   3.333e-4,   0.002,     0.01,   1, 0,  0,  -0.001),\
                                (3934.78,'CaII',3900.,3960.,'CaII3934_2',     1,  1e-3,   3.333e-4,   0.002,     0.01,   1, 0,  0,  -0.001),\

                                (3728.48,'OII',3650.,3800.,'OII3728',       1,  1e-3,   3.333e-4,   0.002,     0.01,   1,  1,  0,  0.001),\

                                (3426.84,'NeV',3380.,3480.,'NeV3426',       1,  1e-3,   3.333e-4,   0.002,     0.01,   1,  0,  0,  0.001),\
                                (3426.84,'NeV',3380.,3480.,'NeV3426_br',    1,  5e-3,   0.0025,     0.01,       0.01,   1,  0,  0,  0.001),\
                                ],\
                             formats='float32,a20,float32,float32,a20,float32,float32,float32,float32,\
                             float32,float32,float32,float32,float32',\
                             names='lambda,compname,minwav,maxwav,linename,ngauss,inisig,minsig,maxsig,voff,vindex,windex,findex,fvalue')
        #------header-----------------
        hdr = fits.Header()
        hdr['lambda'] = 'Vacuum Wavelength in Ang'
        hdr['minwav'] = 'Lower complex fitting wavelength range'
        hdr['maxwav'] = 'Upper complex fitting wavelength range'
        hdr['ngauss'] = 'Number of Gaussians for the line'
        hdr['inisig'] = 'Initial guess of linesigma [in lnlambda]'
        hdr['minsig'] = 'Lower range of line sigma [lnlambda]'
        hdr['maxsig'] = 'Upper range of line sigma [lnlambda]'
        hdr['voff  '] = 'Limits on velocity offset from the central wavelength [lnlambda]'
        hdr['vindex'] = 'Entries w/ same NONZERO vindex constrained to have same velocity'
        hdr['windex'] = 'Entries w/ same NONZERO windex constrained to have same width'
        hdr['findex'] = 'Entries w/ same NONZERO findex have constrained flux ratios'
        hdr['fvalue'] = 'Relative scale factor for entries w/ same findex'
        #------save line info-----------
        hdu = fits.BinTableHDU(data=newdata,header=hdr,name='data')
        hdu.writeto('./qsopar_local.fits',overwrite=True)
