DipFit
=========
A program DipFit is developed for the analysis of RIDME spectra acquired on the spin pairs that consist of one isotropic S = 1/2 center and one anisotropic S = 1/2 center.

The program support two modes: a simulation mode and a fitting mode.

In the simulation mode, the RIDME spectrum is simulated for a specified spin pair with a pre-defined inter-spin distance distribution, P(r), and two angular distributions, P(ξ) and P(φ). 

In the fitting mode, the experimental RIDME spectrum of a certain spin pair is fitted by means of genetic algorithm using an inter-spin  distance distribution, P(r), and two angular distributions, P(ξ) and P(φ), as fitting parameters. To reduce the number of fitting parameters, all three distributions are parameterized. Such parameterization is done by assuming that all distributions have the Gaussian shape, i.e., each distribution can be described by two parameters – a mean value and a standard deviation. The correlation between the values of r, ξ, and φ is neglected. In cases, when the RIDME experiment was performed at liquid helium temperatures, the temperature of the experiment can be considered as an additional fitting parameter.

Further description of the program can be found in our paper (see below).

General Information
=========
The source code of the program is written using python and the libraries 'numpy' and 'matplotlib'. 

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite: D. Abdullin, H. Matsuoka, M. Yulikov, N. Fleck, C. Klein, S. Spicher, G. Hagelueken, S. Grimme, A. Lützen, O. Schiemann, Pulsed EPR Dipolar Spectroscopy under the Breakdown of the High-Field Approximation: The High-Spin Iron(III) Case, submitted.
