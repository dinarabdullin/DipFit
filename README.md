DipFit
=========
The program DipFit allows extracting distance and angular distributions out of the RIDME data, which is acquired on the spin pairs consisting of one isotropic spin-1/2 center and one anisotropic spin-1/2 center.

The program support two modes: a simulation mode and a fitting mode.

In the simulation mode, the RIDME spectrum or the RIDME time trace is simulated for a specified spin pair with a pre-defined inter-spin distance distribution, P(r), and two angular distributions, P(ξ) and P(φ). 

In the fitting mode, the experimental RIDME spectrum of the experimetal RIDME time trace is fitted by means of genetic algorithm using an inter-spin  distance distribution, P(r), and two angular distributions, P(ξ) and P(φ), as fitting parameters. To reduce the number of fitting parameters, all three distributions are assumed to have the Gaussian shape, i.e., each distribution can be described by two parameters – a mean value and a standard deviation. The correlation between the values of r, ξ, and φ is neglected. In cases when the RIDME experiment is performed at liquid helium temperatures, the temperature of the experiment has to be included as an additional fitting parameter.

Further description of the program can be found in our paper (see below).

General Information
=========
The source code of the program is written in Python (2.7 and 3.7) using the libraries numpy, scipy, matplotlib and libconf. 

Copyright
=========
This program can be distributed under GNU General Public License.

If you use this code please cite: D. Abdullin, H. Matsuoka, M. Yulikov, N. Fleck, C. Klein, S. Spicher, G. Hagelueken, S. Grimme, A. Lützen, O. Schiemann, Pulsed EPR Dipolar Spectroscopy under the Breakdown of the High-Field Approximation: The High-Spin Iron(III) Case,
Chem. Eur. J. 2019, doi: 10.1002/chem.201900977
