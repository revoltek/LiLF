# Short explanation of the SkyModels

The models are:
* *calib-hba.skymodel* : Calibrator models for HBA, taken from prefactor.
* *calib-highres.skymodel* : High-resolution calibrator models for LBA.
* *calib-simple.skymodel* : Simple skymodels for calibrators in LBA.
* *demix_all.skymodel* : A-Team models for demixing in LBA.
* *A-Team_lowres.skymodel* : lowres skymodel for A-Team clipping in HBA in prefactor

For convenience, all models are also available in the compiled version.

To re-create compiled version:
% makesourcedb outtype="blob" format="<" in="xxx" out="xxx"

Other models are:
* *models/lotss_dr3_gaus_110325.skymodel* : Gaussian component list of LoTSS dr3 prepared on 11-03-2025
* *3CRR/* : all 3CRR models found by Boxelaar et al. 2025 (in prep.)