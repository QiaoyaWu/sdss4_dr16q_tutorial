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
- `download`: folder to keep downloaded output fits file
- `op`: folder to keep QA plot

## Catalog format description
# First extension
Main properties of SDSS DR16 quasars.
| Column | Format | Unit | Description |
| ---- | --- | --- | ------ |
| SDSS_NAME | STRING |  | SDSS DR16 designation (J2000) |
| PLATE | LONG64| | Spectroscopic plate number |
| MJD | LONG64 | | Spectroscopic MJD |
| FIBERID | LONG64 | | Spectroscopic fiber number |
| RA | DOUBLE | degree | Right ascension (J2000) |
| DEC | DOUBLE | degree | Declination (J2000) |
| OBJID | STRING | | PLATE-MJD-FIBERID: PyQSOFit output name |
| IF_BOSS_SDSS | STRING | | Source of the input spectrum: BOSS or SDSS |
| Z_DR16Q | DOUBLE | | Best redshift provided by DR16Q |
| SOURCE_Z_DR16Q | DOUBLE | | Source for DR16Q redshift from Lyke et al. (2020) |
| Z_FIT | DOUBLE | | Input redshift for QSOFit; can differ from Z DR16Q |
| Z_SYS | DOUBLE | | Systemic redshift |
| Z_SYS_ERR | DOUBLE | | Uncertainties of systemic redshift |
| EBV | DOUBLE | | Milky Way extinction E(B − V) from Schlegel et al. (1998) and scaled to match the results in Schlafly & Finkbeiner (2011) |
| SN_MEDIAN_ALL | DOUBLE | | Median S/N per pixel of the raw spectrum |
| CONTI_FIT | DOUBLE[5] | | Best-fit parameters for the continuum model (PL+poly) |
| CONTI_FIT_ERR | DOUBLE[5] | | Uncertainties in the best-fit continuum parameters |
| CONTI_STAT | DOUBLE[2] | | Continuum fitting pixel number, reduced $\chi^2$ |
| FEII_UV | DOUBLE[3]| | Best-fit parameters for the UV FeII model |
| FEII_UV_ERR | DOUBLE[3] | | Uncertainties in the best-fit UV FeII model |
| FEII_UV_EW | DOUBLE | | Rest-frame equivalent width of UV FeII within 2250-2650 Å |
| FEII_UV_EW_ERR | DOUBLE | | Uncertainties in REW_FE_2250_2650 |
| FEII_OPT | DOUBLE[3] | | Best-fit parameters for the optical FeII model |
| FEII_OPT_ERR | DOUBLE[3] | | Uncertainties in the best-fit optical FeII model |
| FEII_OPT_EW | DOUBLE | | Rest-frame equivalent width of optical FeII within 4434-4680Å|
| FEII_OPT_EW_ERR | DOUBLE | | Uncertainties in REW_FE_4434_4684 |
| LOGL1350 | DOUBLE | $\rm [erg\,s^{-1}]$ | Continuum luminosity at rest-frame 1350Å |
| LOGL1350_ERR | DOUBLE |$\rm [erg\,s^{-1}]$ | Uncertainty in LOGL1350 |
| LOGL1700 | DOUBLE | $\rm [erg\,s^{-1}]$ | Continuum luminosity at rest-frame 1700Å |
| LOGL1700_ERR | DOUBLE | $\rm [erg\,s^{-1}]$ | Uncertainty in LOGL1700 |
| LOGL2500 | DOUBLE | $\rm [erg\,s^{-1}]$ | Continuum luminosity at rest-frame 2500Å |
| LOGL2500_ERR | DOUBLE | $\rm [erg\,s^{-1}]$ | Uncertainty in LOGL2500|
| LOGL3000 | DOUBLE | $\rm [erg\,s^{-1}]$ | Continuum luminosity at rest-frame 3000Å |
| LOGL3000_ERR | DOUBLE | $\rm [erg\,s^{-1}]$ | Uncertainty in LOGL3000 |
| LOGL5100 | DOUBLE | $\rm [erg\,s^{-1}]$ | Continuum luminosity at rest-frame 5100Å |
| LOGL5100_ERR | DOUBLE | $\rm [erg\,s^{-1}]$ | Uncertainty in LOGL5100 |
| | | Å,Å,$\rm 10^{-17}erg\,s^{-1}cm^{-2}$,$\rm [erg\,s^{-1}]$,$km\,s^{-1}$,Å | peak wavelength, 50% flux centoid wavelength, flux,$\log L_{\rm line}$, FWHM, rest-frame equivalent width |
| HALPHA | DOUBLE[6] | ... | For the entire H $\alpha$ profile (narrow and broad lines combined) |
| HALPHA_BR | DOUBLE[6] | ... | For the broad H $\alpha$ profile |
| NII6585 | DOUBLE[6] | ... | For the narrow [NII] $\lambda$ 6584 component|
| SII6718 | DOUBLE[6] | ... | For the entire H $\beta$ profile (narrow and broad lines combined) |
| HBETA | DOUBLE[6] | ... | For the broad H $\beta$ component |
| HBETA_BR | DOUBLE[6] | ... | |
| HEII4687 | DOUBLE[6] | ... | |
| HEII4687_BR | DOUBLE[6] | ... | |
| OIII5007 | DOUBLE[6] | ... | |
| OIII5007C | DOUBLE[6] | ... | |
| CAII3934 | DOUBLE[6] | ... | For the Ca II K absorption line |
| OII3728 | DOUBLE[6] | ... | |
| NEV3426 | DOUBLE[6] | ... | |
| MGII | DOUBLE[6] | ... | |
| MGII_BR | DOUBLE[6] | ... | |
| CIII_ALL | DOUBLE[6] | ... | For the entire CIII] complex (CIII], SiIII], AlIII) |
| CIII_BR | DOUBLE[6] | ... | For the broad CIII] component |
| SIIII1892 | DOUBLE[6] | ... | |
| ALIII1857 | DOUBLE[6] | ... | |
| NIII1750 | DOUBLE[6] | ... | |
| CIV | DOUBLE[6] | ... | |
| HEII1640 | DOUBLE[6] | ... | |
| HEII1640_BR | DOUBLE[6] | ... | |
| SIIV_OIV | DOUBLE[6] | ... | For the 1400Å complex|
| OI1304 | DOUBLE[6] | ... | |
| LYA | DOUBLE[6] | ... | |
| NV1240 | DOUBLE[6] | ... | |
| | | Å,Å,$\rm 10^{-17}erg\,s^{-1}cm^{-2}$,$\rm [erg\,s^{-1}]$,$km\,s^{-1}$,Å | Uncertainties in peak wavelength, 50% flux centoid wavelength, flux,$\log L_{\rm line}$, FWHM, rest-frame equivalent width |
| HALPHA_ERR | DOUBLE[6] | ... | |
| HALPHA_BR_ERR | DOUBLE[6] | ... | |
| NII6585_ERR | DOUBLE[6] | ... | |
| SII6718_ERR | DOUBLE[6] | ... | |
| HBETA_ERR | DOUBLE[6] | ... | |
| HBETA_BR_ERR | DOUBLE[6] | ... | |
| HEII4687_ERR | DOUBLE[6] | ... | |
| HEII4687_BR_ERR | DOUBLE[6] | ... | |
| OIII5007_ERR | DOUBLE[6] | ... | |
| OIII5007C_ERR | DOUBLE[6] | ... | |
| CAII3934_ERR | DOUBLE[6] | ... | |
| OII3728_ERR | DOUBLE[6] | ... | |
| NEV3426_ERR | DOUBLE[6] | ... | |
| MGII_ERR | DOUBLE[6] | ... | |
| MGII_BR_ERR | DOUBLE[6] | ... | |
| CIII_ALL_ERR | DOUBLE[6] | ... | |
| CIII_BR_ERR | DOUBLE[6] | ... | |
| SIIII1892_ERR | DOUBLE[6] | ... | |
| ALIII1857_ERR | DOUBLE[6] | ... | |
| NIII1750_ERR | DOUBLE[6] | ... | |
| CIV_ERR | DOUBLE[6] | ... | |
| HEII1640_ERR | DOUBLE[6] | ... | |
| HEII1640_BR_ERR | DOUBLE[6] | ... | |
| SIIV_OIV_ERR | DOUBLE[6] | ... | |
| OI1304_ERR | DOUBLE[6] | ... | |
| LYA_ERR | DOUBLE[6] | ... | |
| NV1240_ERR | DOUBLE[6] | ... | |
| HA_COMP_STAT | DOUBLE[2] |  | |
| HB_COMP_STAT | DOUBLE[2] |  | |
| MGII_COMP_STAT | DOUBLE[2] |  | |
| CIII_COMP_STAT | DOUBLE[2] |  | |
| CIV_COMP_STAT | DOUBLE[2] | | |
| SIIV_COMP_STAT | DOUBLE[2] | | |
| LYA_COMP_STAT | DOUBLE[2] | | |
| CAII_LOC_STAT | DOUBLE[2] | | |
| OII_LOC_STAT | DOUBLE[2] | | |
| NEV_LOC_STAT | DOUBLE[2] | | |
| LOGLBOL | DOUBLE | $\rm [erg s^{-1}]$ | |
| LOGLBOL_ERR | DOUBLE | $\rm [erg s^{-1}]$ | |
| LOGMBH_HB | DOUBLE | | |
| LOGMBH_HB_ERR | DOUBLE | | |
| LOGMBH_MGII | DOUBLE | | |
| LOGMBH_MGII_ERR | DOUBLE | | |
| LOGMBH_CIV | DOUBLE | | |
| LOGMBH_CIV_ERR | DOUBLE | | |
| LOGMBH | DOUBLE | | |
| LOGMBH_ERR | DOUBLE | | |
| LOGLEDD_RATIO | DOUBLE | | |
| LOGLEDD_RATIO_ERR | DOUBLE | | |
| ZSYS_LINES | DOUBLE | | |
| ZSYS_LINES_ERR | DOUBLE | | |


# Second extension
Quasar properties directly from [DR16Q catalog](https://www.sdss.org/dr16/algorithms/qso_catalog/).
| Column | Description |
| --- | --- | 
