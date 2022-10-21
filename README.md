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
| SDSS_NAME | STRING |  | SDSS DR16 designation (J2000) |

# Second extension
Quasar properties directly from [DR16Q catalog](https://www.sdss.org/dr16/algorithms/qso_catalog/).
| Column | Description |
| --- | --- | 
