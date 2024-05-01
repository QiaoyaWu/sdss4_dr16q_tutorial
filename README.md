# A Catalog of Quasar Properties from Sloan Digital Sky Survey Data Release 16


We present a catalog of continuum and emission line properties for 750,414 broad-line quasars included in the Sloan Digital Sky Survey Data Release 16 quasar catalog (DR16Q), measured from optical spectroscopy. These quasars cover broad ranges in redshift $(0.1 \lesssim z \lesssim 6)$ and luminosity $(44\lesssim \log (L_{\rm bol}/{\rm erg\,s^{-1}})\lesssim 48)$, and probe lower luminosities than an earlier compilation of SDSS DR7 quasars. Derived physical quantities such as single-epoch virial black hole masses and bolometric luminosities are also included in this catalog. We present improved systemic redshifts and realistic redshift uncertainties for DR16Q quasars using the measured line peaks and correcting for velocity shifts of various lines with respect to the systemic velocity. About 1%, 1.4%, and 11% of the original DR16Q redshifts deviate from the systemic redshifts by $|\Delta V|>1500\,{\rm km\,s^{-1}}$, $|\Delta V|\in [1000,1500]\,{\rm km\,s^{-1}}$, and $|\Delta V|\in [500,1000]\,{\rm km\,s^{-1}}$, respectively; about $1900$ DR16Q redshifts were catastrophically wrong $(|\Delta V|>10,000\,{\rm km\,s^{-1}})$. We demonstrate the utility of this data product in quantifying the spectral diversity and correlations among physical properties of quasars with large statistical samples. 

Arxiv: [2209.03987](https://arxiv.org/abs/2209.03987)

## Tutorial

Here we show the [example demo notebook](https://github.com/QiaoyaWu/sdss4_dr16q_tutorial/blob/main/sdss4_QSOFit_op_tutorial.ipynb) for users to read the [PyQSOFit](https://github.com/legolason/PyQSOFit) output file.

The catalog, PyQSOFit output files, and supplemental materials can be downloaded at: [http://quasar.astro.illinois.edu/paper_data/DR16Q/](http://quasar.astro.illinois.edu/paper_data/DR16Q/)

Due to space limit, we do not provide individual QA plots. You can genererate QA plots using the output fits file and the Jupyter notebook tutorial. If you need bulk QA plots for a specific set of DR16Q quasars, we will be happy to provide them upon requests. 

Files:
- `sdss4_QSOFit_op_tutorial.ipynb`: jupyter notebook for using the PyQSOFit output fits
- `qsofunc.py`: python script to calculate the continuum model, gaussian line profile, line properties, and etc.
- `qsopar_generate.py`: python script to generate the PyQSOFit input parameters
- `fe_optical.txt` and `fe_uv.txt`: text files for the FeII model
- `PyQSOFit_wqy.py`: modified version for PyQSOFit v1.1 used in this work 
- `download`: folder to keep downloaded output fits file
- `op`: folder to keep QA plot

## Catalog format description
Here we provide a detailed description for the main catalog that has two extensions.

### First extension
Main properties of SDSS DR16 quasars.
| Column | Description |
| ---- | ------ |
| `SDSS_NAME` | SDSS DR16 designation (J2000) |
| `PLATE` | Spectroscopic plate number |
| `MJD` | Spectroscopic MJD |
| `FIBERID` | Spectroscopic fiber number |
| `RA` | Right ascension (J2000) |
| `DEC` | Declination (J2000) |
| `OBJID` | PLATE-MJD-FIBERID: PyQSOFit output name |
| `IF_BOSS_SDSS` | Source of the input spectrum: BOSS or SDSS |
| `Z_DR16Q` | Best redshift provided by DR16Q |
| `SOURCE_Z_DR16Q` | Source for DR16Q redshift from Lyke et al. (2020) |
| `Z_FIT` | Input redshift for QSOFit; can differ from Z DR16Q |
| `Z_SYS` | Systemic redshift |
| `Z_SYS_ERR` | Uncertainties of systemic redshift |
| `EBV` | Milky Way extinction E(B − V) from Schlegel et al. (1998) and scaled to match the results in Schlafly & Finkbeiner (2011) |
| `SN_MEDIAN_ALL` | Median S/N per pixel of the raw spectrum |
| `CONTI_FIT` | Best-fit parameters for the continuum model (PL+poly) |
| `CONTI_FIT_ERR` | Uncertainties in the best-fit continuum parameters |
| `CONTI_STAT` | Continuum fitting pixel number, reduced $\chi^2$ |
| `FEII_UV` | Best-fit parameters for the UV FeII model |
| `FEII_UV_ERR` | Uncertainties in the best-fit UV FeII model |
| `FEII_UV_EW` | Rest-frame equivalent width of UV FeII within 2250-2650 Å |
| `FEII_UV_EW_ERR` | Uncertainties in REW_FE_2250_2650 |
| `FEII_OPT` | Best-fit parameters for the optical FeII model |
| `FEII_OPT_ERR` | Uncertainties in the best-fit optical FeII model |
| `FEII_OPT_EW` | Rest-frame equivalent width of optical FeII within 4434-4680Å|
| `FEII_OPT_EW_ERR` | Uncertainties in REW_FE_4434_4684 |
| `LOGL1350` | Continuum luminosity at rest-frame 1350Å |
| `LOGL1350_ERR` | Uncertainty in LOGL1350 |
| `LOGL1700` | ... |
| `LOGL1700_ERR` | ... |
| `LOGL2500` | ... |
| `LOGL2500_ERR` | ... |
| `LOGL3000` | ... |
| `LOGL3000_ERR` | ... |
| `LOGL5100` | ... |
| `LOGL5100_ERR` | ... |
|  | peak wavelength, 50% flux centoid wavelength, flux, logL of lines, FWHM, rest-frame equivalent width |
| `HALPHA` | For the entire Hα profile (narrow and broad lines combined) |
| `HALPHA_BR` | For the broad Hα profile |
| `NII6585` | For the narrow [NII] λ6584 component|
| `SII6718` | ... |
| `HBETA` | For the entire Hβ profile (narrow and broad lines combined) |
| `HBETA_BR` | For the broad Hβ component |
| `HEII4687` | For the entire HeII λ4687 profile (narrow and broad lines combined) |
| `HEII4687_BR` | For the broad HeII λ4687 component | 
| `OIII5007` | For the entire [OIII] λ5007 profile | 
| `OIII5007C` | For the core [OIII] λ5007 profile | 
| `CAII3934` | For the Ca II K absorption line |
| `OII3728` | ... |
| `NEV3426` | For the entire NeV λ3426 profile (narrow and broad lines combined) |
| `MGII` | For the entire MgII profile (narrow and broad lines combined) |
| `MGII_BR` | For the broad MgII component | 
| `CIII_ALL` | For the entire CIII] complex (CIII], SiIII], AlIII) |
| `CIII_BR` | For the broad CIII] component |
| `SIIII1892` | ... |
| `ALIII1857` | ... |
| `NIII1750` | ... |
| `CIV` | ... |
| `HEII1640` | For the entire HeII λ1640 profile (narrow and broad lines combined) |
| `HEII1640_BR` | For the broad HeII λ1640 component | 
| `SIIV_OIV` | For the 1400Å complex|
| `OI1304` | ... |
| `LYA` | ... |
| `NV1240` | ... |
| | Uncertainties in peak wavelength, 50% flux centoid wavelength, flux, logL of lines, FWHM, rest-frame equivalent width |
| `HALPHA_ERR` | ... |
| `HALPHA_BR_ERR` | ... |
| `NII6585_ERR` | ... |
| `SII6718_ERR` | ... |
| `HBETA_ERR` | ... |
| `HBETA_BR_ERR` | ... |
| `HEII4687_ERR` | ... |
| `HEII4687_BR_ERR` | ... |
| `OIII5007_ERR` | ... |
| `OIII5007C_ERR` | ... |
| `CAII3934_ERR` | ... |
| `OII3728_ERR` | ... |
| `NEV3426_ERR` | ... |
| `MGII_ERR` | ... |
| `MGII_BR_ERR` | ... |
| `CIII_ALL_ERR` | ... |
| `CIII_BR_ERR` | ... |
| `SIIII1892_ERR` | ... |
| `ALIII1857_ERR` | ... |
| `NIII1750_ERR` | ... |
| `CIV_ERR` | ... |
| `HEII1640_ERR` | ... |
| `HEII1640_BR_ERR` | ... |
| `SIIV_OIV_ERR` | ... |
| `OI1304_ERR` | ... |
| `LYA_ERR` | ... |
| `NV1240_ERR` | ... |
| | |
| `HA_COMP_STAT` | Complex line pixel number, reduced $\chi^2$|
| `HB_COMP_STAT` | ... |
| `MGII_COMP_STAT` | ... |
| `CIII_COMP_STAT` | ... |
| `CIV_COMP_STAT` | ... |
| `SIIV_COMP_STAT` | ... |
| `LYA_COMP_STAT` | ... |
| `CAII_LOC_STAT` | Local line pixel number, reduced $\chi^2$ |
| `OII_LOC_STAT` | ... |
| `NEV_LOC_STAT` | ... |
| `LOGLBOL` | Bolometric luminosity |
| `LOGLBOL_ERR` | Uncertainties in bolometric luminosity |
| `LOGMBH_HB` | Single-epoch BH mass based on Hβ |
| `LOGMBH_HB_ERR` | Uncertainties in LOGMBH_HB |
| `LOGMBH_MGII` | Single-epoch BH mass based on MgII |
| `LOGMBH_MGII_ERR` | ... |
| `LOGMBH_CIV` | Single-epoch BH mass based on CIV |
| `LOGMBH_CIV_ERR` | ... |
| `LOGMBH` | Fiducial single-epoch BH mass |
| `LOGMBH_ERR` | ... |
| `LOGLEDD_RATIO` | Eddington ratio based on fiducial single-epoch BH mass |
| `LOGLEDD_RATIO_ERR` | Uncertainties in LOGLEDD_RATIO |
| `ZSYS_LINES` | Systematic redshift from individual lines in the order of Hβ_BR, [OIII]5007, CaII3934, [OII]3728, MgII, CIII], CIV, SiIV |
| `ZSYS_LINES_ERR` | Uncertainties in systematic redshift from individual lines |
| `FHOST_5100` | Flux ratio between the host galaxy and total flux at 5100 $\rm \AA$

### Second extension
Quasar properties directly from [DR16Q catalog](https://www.sdss.org/dr16/algorithms/qso_catalog/).
| Column | Description |
| --- | --- | 
| `Z_DLA` | Absorber redshift |
| `NHI_DLA` | Absorber column density |
| `CONF_DLA` | Confidence rating for DLA absorbers |
| `BAL_PROB` | probability an object is a BAL quasar. |
| `BI_CIV` | BALnicity index (BI) for CIV |
| `ERR_BI_CIV` | BI uncertainty for CIV |
| `AI_CIV` | Absorption index for CIV|
| `ERR_AI_CIV` | AI uncertainty for CIV |
| `BI_SIIV` | ... |
| `ERR_BI_SIIV` | ... |
| `AI_SIIV` | ... |
| `ERR_AI_SIIV` | ... |
| `NSPEC_SDSS` | Number of additional spectra for an object from SDSS-I/II  |
| `NSPEC_BOSS` | Number of additional spectra for an object from BOSS/eBOSS|
| `NSPEC` | Total number of additional spectroscopic observations for an object |
| `PSFFLUX` | PSF flux for each of the five SDSS bands: u, g, r, i, and z |
| `PSFFLUX_IVAR` | Inverse variance for PSFFLUX |
| `PSFMAG` | Inverse hyperbolic sine AB magnitudes for each of the five SDSS bands: u, g, r, i, and z |
| `PSFMAGERR` | Inverse variance for PSFMAG |
| `EXTINCTION` | The Galactic extinction values from Schlafly & Finkbeiner (2011a) for the five SDSS bands |
| `M_I` | Absolute i -band magnitude corrected for extinction |
| `SN_MEDIAN_ALL` | Median S/N for all good pixels in the five SDSS bands|
| `GALEX_MATCHED` | Matching flag for objects in the forced-photometry SDSS-DR8/GALEX catalog|
| `FUV` | GALEX FUV flux |
| `FUV_IVAR` | Inerse variance for FUV |
| `NUV` | GALEX NUV flux | 
| `NUV_IVAR` | Inerse variance for NUV |
| `UKIDSS_MATCHED` | Matching flag for objects in the forced-photometry SDSS-DR8/UKIDSS catalog |
| `YFLUX` | Flux density and error in Y band |
| `YFLUX_ERR` | Flux density error in Y band |
| `JFLUX` | ... |
| `JFLUX_ERR` | ... |
| `HFLUX` | ... |
| `HFLUX_ERR` | ... |
| `KFLUX` | ... |
| `KFLUX_ERR` | ... |
| `W1_FLUX` | W1-band (3.4 μm) WISE flux |
| `W1_FLUX_IVAR` | Inverse variance of W1_FLUX |
| `W1_MAG` | W1-band magnitude in Vega magnitude system|
| `W1_MAG_ERR` | Magnitude error in W1_MAG |
| `W1_CHI2` | W1-band profile-weighted $\chi^2$ goodness of fit weighted by the point-spread function in the WISE image|
| `W1_FLUX_SNR` | W1-band flux signal-to-noise ratio |
| `W1_SRC_FRAC` | W1-band profile-weighted number of WISE exposure|
| `W1_EXT_FLUX` | W1-band profile-weighted flux from other sources within the PSF of this source |
| `W1_EXT_FRAC` | Profile-weighted fraction of flux from external sources for this source |
| `W1_NPIX` | W1-band number of pixel |
| `W2_FLUX` | W2-band (4.6 μm) WISE flux |
| `W2_FLUX_IVAR` | ... |
| `W2_MAG` | ... |
| `W2_MAG_ERR` | ... |
| `W2_CHI2` | ... |
| `W2_FLUX_SNR` | ... |
| `W2_SRC_FRAC` | ... |
| `W2_EXT_FLUX` | ... |
| `W2_EXT_FRAC` | ... |
| `W2_NPIX` | ... |
| `FIRST_MATCHED` | Matching flag for objects in the FIRST catalog |
| `FIRST_FLUX` | Peak flux density |
| `FIRST_SNR` | Flux density signal-to-noise ratio |
| `SDSS2FIRST_SEP` | Matching separation in arcsec between the SDSS and FIRST objects|
| `JMAG` | J-band magnitude |
| `JMAG_ERR` | J-band magnitude error |
| `JSNR` | J-band signal-to-noise ratio |
| `JRDFLAG` | J-band photometric read flag |
| `HMAG` | ... |
| `HMAG_ERR` | ... |
| `HSNR` | ... |
| `HRDFLAG` | ... |
| `KMAG` | ... |
| `KMAG_ERR` | ... |
| `KSNR` | ... |
| `KRDFLAG` | ... |
| `SDSS2MASS_SEP` | Matching separation in arcsec between SDSS and 2MASS object |
| `2RXS_ID` | 2RXS ID number designating unique objects | 
| `2RXS_RA` | Right ascension of the 2RXS source in decimal degrees (J2000)|
| `2RXS_DEC` | Declination of the 2RXS source in decimal degrees (J2000)|
| `2RXS_SRC_FLUX` | Source flux in the 0.5–2.0 keV band |
| `2RXS_SRC_FLUX_ERR` | Source flux error in the 0.5–2.0 keV band |
| `SDSS2ROSAT_SEP` | Matching separation in arcsec between SDSS and 2RXS objects|
| `XMM_SRC_ID` | XMM-Newton Source designation ID |
| `XMM_RA` | Right ascension  of the XMM-Newton source in decimal degrees (J2000) |
| `XMM_DEC` | Declination of the XMM-Newton source in decimal degrees (J2000) |
| `XMM_SOFT_FLUX` | The X-ray flux for the full energy range (0.2–2.0 keV) |
| `XMM_SOFT_FLUX_ERR` | The X-ray flux error for the full energy range (0.2–2.0 keV) |
| `XMM_HARD_FLUX` |  The X-ray flux for the full energy range (2.0–12.0 keV) |
| `XMM_HARD_FLUX_ERR` | The X-ray flux error for the full energy range (2.0–12.0 keV) |
| `XMM_TOTAL_FLUX` | The total X-ray flux for the full energy range (0.2–12.0 keV) |
| `XMM_TOTAL_FLUX_ERR` | The total X-ray flux error for the full energy range (0.2–12.0 keV) |
| `XMM_TOTAL_LUM` | Total X-ray luminosity for the full energy range (0.2–12.0 keV)|
| `SDSS2XMM_SEP` | Matching separation in arcsec between SDSS and XMM-Newton objects |
