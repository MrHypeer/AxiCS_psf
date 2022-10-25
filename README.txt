Introduction
====================

This package provides the code employed for using compressive sensing with laser-scanning
microscopy, where data reduction is applied directly in the spatial domain, while including
a model of the imaging system to retain optimal spatial resolution.

This work is described in details in the article:

N. Pavillon and N. I. Smith, "Compressed sensing laser scanning microscopy", 
    Opt. Express 24(26), pp. 30038-300052 (2016), doi: 10.1364/OE.24.030038.
 
This package, as well as demonstration data, can be found at the address of the Biophotonics 
Laboratory, Osaka University.
http://biophotonics.ifrec.osaka-u.ac.jp/en/

Credits
====================

This code makes extensive used of the TV solver TVAL3, originally developed by Chengbo Li and
Yin Zhang from Rice University. Please refer to the readme in TVAL3/Solver/readme.txt and
the webpage of the original package http://www.caam.rice.edu/~optimization/L1/TVAL3/ for more details.

The solver is used mostly without core modifications, with the significant changes being inserted as
new measurement matrices functions. Refer to ChangeLog.txt for details.

Usage
====================

The function TVAL3 executes the solver. Wrapper functions for our applications are given for both
2D microscopy images (ConfocalCS) and hyperspectral 3D Raman data (RamanCS). 

Examples of use are provided for different configurations, as detailed below.

NOTE: Confocal and Raman examples require external data which is distributed separately on the 
Biophotonics Laboratory website.

exampleSimulated
-----------------------
Comparison of performance between standard compressed sensing (based on Fourier ensemble or Hadamard
matrices) and the spatial approach, with and without inclusion of the imaging model.

The level of smoothing in the image as well as noise levels can also be adjusted.

exampleConfocal
-----------------------
This script illustrates the reconstruction from real data measured on a confocal microscope, as well
as the use of a practically measured PSF.

exampleRaman
-----------------------
This script illustrates the reconstruction in the case of hyperspectral data from Raman spectroscopy, 
where the TV regularization also contributes in increasing the signal quality from noisy data.
   
Contact Information
=======================

Biophotonics Laboratory
Immunology Frontier Research Center
Osaka University
Japan

Please address enquiries to either 
Nicolas Pavillon        n-pavillon@ifrec.osaka-u.ac.jp
or
Nicholas I. Smith       nsmith@ap.eng.osaka-u.ac.jp

Copyright Notice
====================

TVAL3 is free software, and you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details at
<http://www.gnu.org/licenses/>. 
