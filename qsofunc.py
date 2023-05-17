#!/usr/bin/python3

import os,sys,time
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
#from astropy.modeling.blackbody import blackbody_lambda
from scipy import interpolate
from scipy import integrate
from astropy import constants as const
warnings.filterwarnings("ignore")
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

fe_uv = np.genfromtxt('./fe_uv.txt')
fe_op = np.genfromtxt('./fe_optical.txt')

def F_poly_conti(xval, pp):
    """Fit the continuum with a polynomial component account for the dust reddening with a*X+b*X^2+c*X^3"""
    xval2 = xval-3000.
    yval = 0.*xval2
    for i in range(len(pp)):
        yval = yval+pp[i]*xval2**(i+1)
    return yval

def Balmer_conti(xval, pp):
    """Fit the Balmer continuum from the model of Dietrich+02"""
    lambda_BE = 3646.  # A
    bbflux = blackbody_lambda(xval, pp[1]).value*3.14  # in units of ergs/cm2/s/A
    tau = pp[2]*(xval/lambda_BE)**3
    result = pp[0]*bbflux*(1.-np.exp(-tau))
    ind = np.where(xval > lambda_BE, True, False)
    if ind.any() == True:
        result[ind] = 0.
    return result

def Onegauss(xval, pp):
    term1 = np.exp(- (xval-pp[1])**2/(2.*pp[2]**2))
    yval = pp[0]*term1/(np.sqrt(2.*np.pi)*pp[2])
    return yval

def Manygauss(xval, pp):
    """The multi-Gaussian model used to fit the emission lines, it will call the onegauss function"""
    ngauss = int(pp.shape[0]/3)
    if ngauss != 0:
        yval = 0.
        for i in range(ngauss):
            yval = yval+Onegauss(xval, pp[i*3:(i+1)*3])
        return yval
    else:
        return np.zeros_like(xval)

def Fe_flux_mgii(xval, pp):
    "Fit the UV Fe compoent on the continuum from 1200 to 3500 A based on the Boroson & Green 1992."
    yval = np.zeros_like(xval)
    wave_Fe_mgii = 10**fe_uv[:, 0]
    flux_Fe_mgii = fe_uv[:, 1]*10**15
    Fe_FWHM = pp[1]
    xval_new = xval*(1.0+pp[2])

    ind = np.where((xval_new > 1200.) & (xval_new < 3500.), True, False)
    if np.sum(ind) > 100:
        if Fe_FWHM < 900.0:
            sig_conv = np.sqrt(910.0**2-900.0**2)/2./np.sqrt(2.*np.log(2.))
        else:
            sig_conv = np.sqrt(Fe_FWHM**2-900.0**2)/2./np.sqrt(2.*np.log(2.))  # in km/s
        # Get sigma in pixel space
        sig_pix = sig_conv/106.3  # 106.3 km/s is the dispersion for the BG92 FeII template
        khalfsz = np.round(4*sig_pix+1, 0)
        xx = np.arange(0, khalfsz*2, 1)-khalfsz
        kernel = np.exp(-xx**2/(2*sig_pix**2))
        kernel = kernel/np.sum(kernel)

        flux_Fe_conv = np.convolve(flux_Fe_mgii, kernel, 'same')
        tck = interpolate.splrep(wave_Fe_mgii, flux_Fe_conv)
        yval[ind] = pp[0]*interpolate.splev(xval_new[ind], tck)
    return yval

def Fe_flux_balmer(xval, pp):
    "Fit the optical FeII on the continuum from 3686 to 7484 A based on Vestergaard & Wilkes 2001"
    yval = np.zeros_like(xval)

    wave_Fe_balmer = 10**fe_op[:, 0]
    flux_Fe_balmer = fe_op[:, 1]*10**15
    ind = np.where((wave_Fe_balmer > 3686.) & (wave_Fe_balmer < 7484.), True, False)
    wave_Fe_balmer = wave_Fe_balmer[ind]
    flux_Fe_balmer = flux_Fe_balmer[ind]
    Fe_FWHM = pp[1]
    xval_new = xval*(1.0+pp[2])
    ind = np.where((xval_new > 3686.) & (xval_new < 7484.), True, False)
    if np.sum(ind) > 100:
        if Fe_FWHM < 900.0:
            sig_conv = np.sqrt(910.0**2-900.0**2)/2./np.sqrt(2.*np.log(2.))
        else:
            sig_conv = np.sqrt(Fe_FWHM**2-900.0**2)/2./np.sqrt(2.*np.log(2.))  # in km/s
        # Get sigma in pixel space
        sig_pix = sig_conv/106.3  # 106.3 km/s is the dispersion for the BG92 FeII template
        khalfsz = np.round(4*sig_pix+1, 0)
        xx = np.arange(0, khalfsz*2, 1)-khalfsz
        kernel = np.exp(-xx**2/(2*sig_pix**2))
        kernel = kernel/np.sum(kernel)
        flux_Fe_conv = np.convolve(flux_Fe_balmer, kernel, 'same')
        tck = interpolate.splrep(wave_Fe_balmer, flux_Fe_conv)
        yval[ind] = pp[0]*interpolate.splev(xval_new[ind], tck)
    return yval

def continuum_all(wave_val, conti_ip_val):
    return conti_ip_val[0]*(wave_val/3000.0)**conti_ip_val[1]\
        + F_poly_conti(wave_val, conti_ip_val[2:5])\
        + Fe_flux_mgii(wave_val, conti_ip_val[5:8])\
        + Fe_flux_balmer(wave_val, conti_ip_val[8:])

def continuum(wave_val, conti_ip_val):
    return conti_ip_val[0]*(wave_val/3000.0)**conti_ip_val[1]\
        + F_poly_conti(wave_val, conti_ip_val[2:5])

def get_line_prop_n(dl, center, line_params, conti_params):
    line_params = line_params.astype(float)
    conti_params = conti_params.astype(float)
    if np.any(line_params[::3]):
        c = 299792.458  # km/s
        n_gauss = int(len(line_params)/3)

        cen = np.zeros(n_gauss)
        sig = np.zeros(n_gauss)
        for i in range(n_gauss):
            cen[i] = line_params[3*i+1]
            sig[i] = line_params[3*i+2]

        left = min(cen-3*sig)
        right = max(cen+3*sig)
        disp = 1.e-4*np.log(10.)
        npix = int((right-left)/disp)

        xx = np.linspace(left, right, npix)
        yy = Manygauss(xx, line_params)
        contiflux = conti_params[0]*(np.exp(xx)/3000.0)**conti_params[1]\
        + F_poly_conti(np.exp(xx), conti_params[2:5])
        + Fe_flux_mgii(np.exp(xx), conti_params[5:8])\
        + Fe_flux_balmer(np.exp(xx), conti_params[8:])

        ypeak = yy.max()
        ypeak_ind = np.argmax(abs(yy))
        peak = np.exp(xx[ypeak_ind])

        # find the FWHM in km/s
        spline = interpolate.UnivariateSpline(xx, abs(yy)-np.max(abs(yy))/2, s=0)
        if len(spline.roots()) > 0:
            fwhm_left, fwhm_right = spline.roots().min(), spline.roots().max()
            #print(fwhm_left, fwhm_right)
            fwhm = abs(np.exp(fwhm_left)-np.exp(fwhm_right))/center*c
            line_flux = Manygauss(xx, line_params)
            line_wave = np.exp(xx)
            lambda0 = integrate.trapz(line_flux, line_wave)
            lambda1 = integrate.trapz(line_flux*line_wave, line_wave)
            lambda2 = integrate.trapz(line_flux*line_wave*line_wave, line_wave)
            ew = integrate.trapz(np.abs(line_flux/contiflux), line_wave)
            area = lambda0
            logL = np.log10(4*np.pi*dl**2*abs(area)*1e-17)
            sigma = np.sqrt(lambda2/lambda0-(lambda1/lambda0)**2)/center*c
        else:
            fwhm, area, logL, ew, peak = 0., 0., 0., 0., 0.
    else:
        fwhm, area, logL, ew, peak = 0., 0., 0., 0., 0.
    return peak, area, logL, fwhm, ew

def MC_1sigma_err(ip_arr):
    return (np.percentile(ip_arr, 84)-np.percentile(ip_arr, 16))/2

def set_mpl_style(fsize=15, tsize=18, tdir='in', major=5.0, minor=3.0, lwidth=1.8, lhandle=2.0):
    """Function to set MPL style"""

    plt.style.use('default')
    plt.rcParams['text.usetex'] = False
    plt.rcParams["axes.axisbelow"] = False
    plt.rcParams['font.size'] = fsize
    plt.rcParams['legend.fontsize'] = tsize
    plt.rcParams['xtick.direction'] = tdir
    plt.rcParams['ytick.direction'] = tdir
    plt.rcParams['ytick.right'] = True
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.major.size'] = major
    plt.rcParams['xtick.minor.size'] = minor
    plt.rcParams['ytick.major.size'] = major
    plt.rcParams['ytick.minor.size'] = minor
    plt.rcParams['xtick.major.width'] = lwidth
    plt.rcParams['xtick.minor.width'] = lwidth
    plt.rcParams['ytick.major.width'] = lwidth
    plt.rcParams['ytick.minor.width'] = lwidth
    plt.rcParams['axes.linewidth'] = lwidth
    plt.rcParams['legend.handlelength'] = lhandle
    return

def get_line_prop(linename, conti_arr, qso_data, qso_data_mc):
    op_linelist = np.array(['HALPHA', 'HALPHA_BR', 'NII6585', 'SII6718', \
                            'HBETA', 'HBETA_BR', 'HEII4687', 'HEII4687_BR', \
                            'OIII5007', 'OIII5007C', 'CAII3934', 'OII3728', 'NEV3426', \
                            'MGII', 'MGII_BR', 'CIII_ALL', 'CIII_BR', \
                            'SIIII1892', 'ALIII1857', 'NIII1750', \
                            'CIV', 'HEII1640', 'HEII1640_BR', \
                            'SIIV_OIV', 'OI1304', 'LYA', 'NV1240'])
    op_linelist_wv = np.array([6564.61, 6564.61, 6585.28, 6718.29, \
                               4862.68, 4862.68, 4687.02, 4687.02, \
                               5008.24, 5008.24, 3934.78, 3728.48, 3426.84, \
                               2798.75, 2798.75, 1908.73, 1908.73, \
                               1892.03, 1857.40 , 1750.26, \
                               1549.06, 1640.42, 1640.42, \
                               (1402.06+1396.76)/2, 1304.35, 1215.67, 1240.42])
    op_sublinelist = np.array([['Halpha_br_1', 'Halpha_br_2', 'Halpha_br_3', 'Halpha_na_1'], ['Halpha_br_1', 'Halpha_br_2', 'Halpha_br_3'], ['NII6585_1'], ['SII6718_1'],\
                               ['Hbeta_br_1', 'Hbeta_br_2', 'Hbeta_br_3', 'Hbeta_na_1'], ['Hbeta_br_1', 'Hbeta_br_2', 'Hbeta_br_3'], ['HeII4687_br_1', 'HeII4687_na_1'], ['HeII4687_br_1'], \
                               ['OIII5007w_1', 'OIII5007c_1'], ['OIII5007c_1'], ['CaII3934_1_1', 'CaII3934_2_1'], ['OII3728_1'], ['NeV3426_1', 'NeV3426_br_1'], \
                               ['MgII_br_1', 'MgII_br_2', 'MgII_na_1'], ['MgII_br_1', 'MgII_br_2'], ['CIII_br1_1', 'CIII_br2_1', 'SiIII1892_1', 'AlIII1857_1'], ['CIII_br1_1', 'CIII_br2_1'], \
                               ['SiIII1892_1'], ['AlIII1857_1'], ['NIII1750_1'], \
                               ['CIV_br_1', 'CIV_br_2', 'CIV_br_3'], ['HeII1640_1', 'HeII1640_br_1'], ['HeII1640_br_1'], \
                               ['SiIV_OIV1_1', 'SiIV_OIV2_1'], ['OI1304_1'], ['Lya_br_1', 'Lya_br_2', 'Lya_br_3'], ['NV1240_1']])
    op_sublinelist_flat = np.array([op_sublinelist[i][j] for i in range(len(op_sublinelist)) for j in range(len(op_sublinelist[i]))])
    op_sublinelist_flat_order = np.array([i for i in range(len(op_sublinelist)) for j in range(len(op_sublinelist[i]))])

    dl = cosmo.luminosity_distance(qso_data.Z_FIT[0]).value*10**6*3.086*10**18
    opl = np.where((op_linelist==linename))[0]
    if len(opl) == 1:
        opl = opl[0]
        subline_lst = op_sublinelist[opl]
        gauss_para_name = np.array([])
        gauss_para_nameMC = np.array([])
        for sbl in range(len(subline_lst)):
            if subline_lst[sbl]+'_scale' in qso_data.names:
                gauss_para_name = np.append(gauss_para_name, \
                                              [subline_lst[sbl]+'_scale', subline_lst[sbl]+'_centerwave', subline_lst[sbl]+'_sigma'])
                gauss_para_nameMC = np.append(gauss_para_nameMC, \
                                              [subline_lst[sbl]+'_scale_MC', subline_lst[sbl]+'_centerwave_MC', subline_lst[sbl]+'_sigma_MC'])
        gauss_para = np.array([qso_data[name][0] for name in gauss_para_name])
        line_pp_arr = get_line_prop_n(dl, op_linelist_wv[opl], gauss_para, conti_arr)
        line_pp_err_arr = np.copy(line_pp_arr)
        MC = len(qso_data_mc)
        pp_mc_arr = np.zeros((MC,5))
        for nmc in range(MC):
            gauss_para_tmp = np.array([qso_data_mc[name][nmc] for name in gauss_para_nameMC])
            pp_mc_arr[nmc,:] = get_line_prop_n(dl, op_linelist_wv[opl], gauss_para_tmp, conti_arr)
        for pp in range(5):
            line_pp_err_arr[pp] = MC_1sigma_err(pp_mc_arr[:,pp])
        print(linename)
        print('Peak wavelength: ', str(np.around(line_pp_arr[0], 2))+'+-'+str(np.around(line_pp_err_arr[0], 2)), )
        print('Flux: ', str(np.around(line_pp_arr[1], 2))+'+-'+str(np.around(line_pp_err_arr[1], 2)), )
        print('LogL: ', str(np.around(line_pp_arr[2], 2))+'+-'+str(np.around(line_pp_err_arr[2], 2)), )
        print('FWHM: ', str(np.around(line_pp_arr[3], 2))+'+-'+str(np.around(line_pp_err_arr[3], 2)), )
        print('EW: ', str(np.around(line_pp_arr[4], 2))+'+-'+str(np.around(line_pp_err_arr[4], 2)))
        print('')
        return line_pp_arr, line_pp_err_arr
    else:
        print('Please type in the right line name!')
