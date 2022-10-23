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
| --- | --- | --- | --- |
| SDSS_NAME | STRING |  | SDSS DR16 designation (J2000) |
| PLATE | LONG64| | |
| MJD | LONG64 | | |
| FIBERID | LONG64 | | |
| RA | DOUBLE | | |
| DEC | DOUBLE | | |
| OBJID | STRING | | |
| IF_BOSS_SDSS | STRING | | |
| Z_DR16Q | DOUBLE | | |
| SOURCE_Z_DR16Q | DOUBLE | | |
| Z_FIT | DOUBLE | | |
| Z_SYS | DOUBLE | | |
| Z_SYS_ERR | DOUBLE | | |
| EBV | DOUBLE | | |
| SN_MEDIAN_ALL | DOUBLE | | |
| CONTI_FIT | DOUBLE[5] | | |
| CONTI_FIT_ERR | DOUBLE[5] | | |
| CONTI_STAT | DOUBLE[2] | | |
| FEII_UV | DOUBLE[3]| | |
| FEII_UV_ERR | | | |
| FEII_UV_EW | | | |
| FEII_UV_EW_ERR | | | |
| FEII_OPT | | | |
| FEII_OPT_ERR | | | |
| FEII_OPT_EW | | | |
| FEII_OPT_EW_ERR | | | |
| LOGL1350 | | | |
| LOGL1350_ERR | | | |
| LOGL1700 | | | |
| LOGL1700_ERR | | | |
| LOGL2500 | | | |
| LOGL2500_ERR | | | |
| LOGL3000 | | | |
| LOGL3000_ERR | | | |
| LOGL5100 | | | |
| LOGL5100_ERR | | | |
| HALPHA | | | |
| HALPHA_BR | | | |
| NII6585 | | | |
| SII6718 | | | |
| HBETA | | | |
| HBETA_BR | | | |
| HEII4687 | | | |
| HEII4687_BR | | | |
| OIII5007 | | | |
| OIII5007C | | | |
| CAII3934 | | | |
| OII3728 | | | |
| NEV3426 | | | |
| MGII | | | |
| MGII_BR | | | |
| CIII_ALL | | | |
| CIII_BR | | | |
| SIIII1892 | | | |
| ALIII1857 | | | |
| NIII1750 | | | |
| CIV | | | |
| HEII1640 | | | |
| HEII1640_BR | | | |
| SIIV_OIV | | | |
| OI1304 | | | |
| LYA | | | |
| NV1240 | | | |
| HALPHA_ERR | | | |
| HALPHA_BR_ERR | | | |
| NII6585_ERR | | | |
| SII6718_ERR | | | |
| HBETA_ERR | | | |
| HBETA_BR_ERR | | | |
| HEII4687_ERR | | | |
| HEII4687_BR_ERR | | | |
| OIII5007_ERR | | | |
| OIII5007C_ERR | | | |
| CAII3934_ERR | | | |
| OII3728_ERR | | | |
| NEV3426_ERR | | | |
| MGII_ERR | | | |
| MGII_BR_ERR | | | |
| CIII_ALL_ERR | | | |
| CIII_BR_ERR | | | |
| SIIII1892_ERR | | | |
| ALIII1857_ERR | | | |
| NIII1750_ERR | | | |
| CIV_ERR | | | |
| HEII1640_ERR | | | |
| HEII1640_BR_ERR | | | |
| SIIV_OIV_ERR | | | |
| OI1304_ERR | | | |
| LYA_ERR | | | |
| NV1240_ERR | | | |
| HA_COMP_STAT | | | |
| HB_COMP_STAT | | | |
| MGII_COMP_STAT | | | |
| CIII_COMP_STAT | | | |
| CIV_COMP_STAT | | | |
| SIIV_COMP_STAT | | | |
| LYA_COMP_STAT | | | |
| CAII_LOC_STAT | | | |
| OII_LOC_STAT | | | |
| NEV_LOC_STAT | | | |
| LOGLBOL | | | |
| LOGLBOL_ERR | | | |
| LOGMBH_HB | | | |
| LOGMBH_HB_ERR | | | |
| LOGMBH_MGII | | | |
| LOGMBH_MGII_ERR | | | |
| LOGMBH_CIV | | | |
| LOGMBH_CIV_ERR | | | |
| LOGMBH | | | |
| LOGMBH_ERR | | | |
| LOGLEDD_RATIO | | | |
| LOGLEDD_RATIO_ERR | | | |
| ZSYS_LINES | | | |
| ZSYS_LINES_ERR | | | |


# Second extension
Quasar properties directly from [DR16Q catalog](https://www.sdss.org/dr16/algorithms/qso_catalog/).
| Column | Description |
| --- | --- | 
