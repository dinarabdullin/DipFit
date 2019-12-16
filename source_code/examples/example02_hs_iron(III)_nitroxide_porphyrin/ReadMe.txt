NAME
--------------------------------------
Example 2


DESCRIPTION
--------------------------------------
In this example, DipFit performs the fitting of an experimental dipolar spectra of a model compound (see Pic_1) which contains two spin centers, a high-spin Fe(III) and a nitroxide.

The ZFS tensor of the high-spin Fe(III) ion is almost axial, D = 6.465 cm-1 and E = 0.020 cm-1. 
The actual principal g-values of the high-spin Fe(III) ion are gxx = gyy = 2.00 and gzz = 1.95.

The principal g-values of the nitroxide spin center are set to gxx = 2.009, gyy = 2.006, and gzz = 2.002.

The parameters of the fitting include the distance between the spin centers (mean value and std), the angle xi (mean value and std), and the temperature of the experiment.


HOW TO RUN  
--------------------------------------
1) Open Terminal or Command Prompt.

2) Navigate to the root directory of the DipFit program:
   cd …/DipFit/source_code
   
3) Run the program by the following command:
   python main.py examples/example02_hs_iron(III)_nitroxide_porphyrin/config_ex02.cfg

4) The results of the DipFit analysis will be saved in the folder:
   .../DipFit/examples/example02_hs_iron(III)_nitroxide_porphyrin/