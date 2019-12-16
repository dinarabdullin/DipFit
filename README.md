DipFit
=========
The program DipFit allows the analysis of Pulsed Dipolar EPR Spectroscopy (PDS) data acquired on spin pairs consisting of one isotropic spin-1/2 center and one anisotropic spin-1/2 center.

The program support three modes: simulation, fitting, validation.

In the simulation mode, a dipolar spectrum or a dipolar time trace is simulated for a spin pair with a predefined inter-spin distance distribution and a predefined relative orientation of spin centers.

In the fitting mode, an experimental PDS spectrum of an experimetal PDS time trace is fitted by means of genetic algorithm. In this case, parameters describing the inter-spin distance distribution and the relative orientation of spin centers are used as fitting parameters. In cases when the PDS experiment is done at liquid helium temperatures, the temperature of the experiment can be used an additional fitting parameter.

In the validation mode, the results of the fitting is validated via determining the precision of the optimized fitting parameters.

Further description of the program can be found in our paper (see below).

General Information
=========
The source code of the program is written in Python 3.7 using the Python libraries numpy, scipy, matplotlib and libconf.

The Windows and Linux executables of the program will be uploaded shortly.

The manual of the program will be published shortly.

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite: D. Abdullin, H. Matsuoka, M. Yulikov, N. Fleck, C. Klein, S. Spicher, G. Hagelueken, S. Grimme, A. Lützen, O. Schiemann, Pulsed EPR Dipolar Spectroscopy under the Breakdown of the High-Field Approximation: The High-Spin Iron(III) Case,
Chem. Eur. J. 2019, 25, 8820-8828.
