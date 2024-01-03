**Longitudinal QSM analysis in PD**

QSM = quantitative susceptibility mapping 
PD = Parkinson's disease

The code here can be used to reproduce the group level voxel-wise and region of interest (ROI) results from: 10.1002/mds.29702

The voxel-wise data required for this, as well as associated statistical maps, can be found here: https://zenodo.org/records/10436155

Extract the contents of the zip file from zenodo into the master directory of this repository.

All data provided here are anonymised. 

**Repo contents** (after zip extraction)

* ROI
  * code - _R scripts to run ROI analysis from the paper_ 
  * data - _ROI QSM values, motor and cognitive scores, and basic demographic info_
  * results - _output .png and .csv files will be put here_
* Wholebrain
  * code - _shell scripts to run voxel-wise analyses from the paper_ 
  * data - _smoothed voxel-wise QSM maps for subjects, as well as group mean and brain masks_
    * design_matrices - _design and contrast matrices required for running voxel-wise analyses_ 
  * results - _voxel-wise stasistical summarizing group level results (these can also be re-generated using the code provided)_

If you have any questions about this repository, please email george.thomas.14@ucl.ac.uk
    
