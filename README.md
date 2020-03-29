DipFit
=========
The program DipFit was developed for the analysis of the Pulsed EPR Dipolar Spectroscopy (PDS) signals, which correspond to spin systems consisting of one isotropic and one anisotropic S = 1/2 centers.

The program has three modes: simulation, fitting, and validation.

In the simulation mode, the PDS time trace / spectrum is simulated using the pre-defined geometric model of a spin system and the spectroscopic parameters of spin centers. 

In the fitting mode, the experimental PDS time trace / spectrum is fitted by means of genetic algorithm. In this case, the geometric model of a spin system is optimized until the simulated PDS time trace / spectrum provides the best fit to the experimental PDS time trace / spectrum.

In the validation mode, the results of the fitting are validated via the calculation of the confidence intervals for the optimized fitting parameters.

Further description of the program can be found in the Manual and the papers below.

General Information
=========
The source code of the program is written in Python 3.7 using the Python libraries numpy, scipy, matplotlib and libconf.

The Windows and Linux executables of the program are available at:
https://github.com/dinarabdullin/DipFit/releases

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite: 
1) D. Abdullin, H. Matsuoka, M. Yulikov, N. Fleck, C. Klein, S. Spicher, G. Hagelueken, S. Grimme, A. Lützen, O. Schiemann, “Pulsed EPR Dipolar Spectroscopy under the Breakdown of the High-Field Approximation: The High-Spin Iron(III) Case”, Chem. Eur. J. 2019, 25, 8820-8828.
2) D. Abdullin, P. Brehm, N. Fleck, S. Spicher, S. Grimme, O. Schiemann, “Pulsed EPR Dipolar Spectroscopy on Spin Pairs with one Highly Anisotropic Spin Center: The Low-Spin Fe(III) Case”, Chem. Eur. J. 2019, 25, 14388-14398.
