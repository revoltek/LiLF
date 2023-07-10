# Short explanation of the SkyModels

The models are:
* *calib-hba.skymodel* : Calibrator models for HBA, taken from prefactor.
* *calib-highres.skymodel* : High-resolution calibrator models for LBA.
* *calib-simple.skymodel* : Simple skymodels for calibrators in LBA.
* *demix_all.skymodel* : A-Team models for demixing in LBA.
* *A-Team_lowres.skymodel* : lowres skymodel for A-Team clipping in HBA in prefactor.

For convenience, all models are also available in the compiled version.

To re-create compiled version:
% makesourcedb outtype="blob" format="<" in="xxx" out="xxx"
