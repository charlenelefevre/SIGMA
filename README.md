# SIGMA
Simple Icy Grain Model for Aggregates

SIGMA is a new flexible tool, inherited from DIANA opacity code (Woitke et al. 2016, Min et al. 2016), to compute dust properties of aggregates made of a mixture of distinct materials, including vacuum and ice mantles. Our goal was to deliver a code able to reproduce classic models from the literature on a very fast computation time, but also to deviate from them. Hence, SIGMA is freely available to the community.

We rely on effective medium theory to compute effective refractive index and Mie theory applied to a distribution of hollow spheres to mimic non–spherical dust shapes. This method successfully reproduces exact results of other methods (like discrete dipole approximation) when selecting carefully the parameters related to dust shape. The core of dust grains is built using Bruggeman rule, while ice mantles are added as surrounding material thanks to Maxwell-Garnett rule.

A related paper was submitted to Astronomy & Astrophysics journal.
The reference will be updated in GitHub once accepted.

If your work makes use of SIGMA or any modified version derived from SIGMA, please cite the following publications :
Lefèvre, C., Min, M., Pagani, P., et al. 2019, A&A, SIGMA: Simple Icy Grain Model for Aggregates (reference in latest update in GitHub)
Woitke, P., Min, M., Pinte, C., et al. 2016, A&A, 586, A103, Consistent dust and gas models for protoplanetary disks. I. Disk shape, dust settling, opacities, and PAHs
Min, M., Rab, C., Woitke, P., Dominik, C., & Ménard, F. 2016, A&A, 585, A13, Multiwavelength optical properties of compact dust aggregates in protoplanetary disks


-------------------------------------------------------------

INSTALL SIGMA:

-------------------------------------------------------------
SIGMA is available via git. It can eiter be downloaded from the website directly or with git. 

git clone https://github.com/charlenelefevre/SIGMA.git

Please try to update it on a regular basis with:
git pull

SIGMA can be installed with gfortran or ifort. Version tested so far are the following:
* gfortran GNU Fortran (GCC) 4.8.5 20150623 (Red at 4.8.5-16, Linux)
* gfortran GNU Fortran (GCC) 6.3.0 (Mac)
* ifort (IFORT) 11.1 20091130 (Linux)
* ifort (IFORT) 18.0.0 20170811 (Mac)

If you want to use ifort please make sure to add in your bashrc:
export ifort="true"
or in your .tcshrc:
setenv ifort "true"

Please try to solve software installation with your computer support before reporting.
--------------------------------------------------------------

To install SIGMA:
cd SIGMA
make clean
make

--------------------------------------------------------------
To run SIGMA:
You can test with default value:
SIGMA -nm 6 -na 100 -verbose

--------------------------------------------------------------
Or obtain more information about options with:
SIGMA -help


