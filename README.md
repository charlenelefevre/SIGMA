# SIGMA
Simple Icy Grain Model for Aggreagtes

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


