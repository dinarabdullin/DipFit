NAME
--------------------------------------
Example 1


DESCRIPTION
--------------------------------------
In this example, the dipolar spectra is simulated for a spin system consisting of a high-spin Fe(III) ion and a nitroxide. 

The ZFS tensor of the high-spin Fe(III) ion is set to be axial, D = 10 cm-1 and E = 0 cm-1. 
The actual principal g-values of the high-spin Fe(III) ion are set to the free electron g-factor.

For simplicity, the principal g-values of the nitroxide spin center are set to the free electron g-factor too.

The distance between the spin centers is set to 2.5 nm. 
The orientation of the distance vector relative to the g-frame of the high-spin Fe(III) ion is described by two angles, xi - the polar angle and phi - the azimuthal angle.
Due to the axial ZFS tensor of the high-spin Fe(III), the phi angle has no effect on the dipolar spectrum and thus was set to 0. 


HOW TO RUN  
--------------------------------------
1) Open Terminal or Command Prompt.

2) Navigate to the root directory of the DipFit program:
   cd …/DipFit/
   
3) Run the program by the following command:
   python main.py /examples/example01_hs_iron(III)_nitroxide/config_ex01_1.cfg (simulates the dipolar spectrum vs angle xi)
   python main.py /examples/example01_hs_iron(III)_nitroxide/config_ex01_2.cfg (simulates the dipolar spectrum for xi = 90°)
   python main.py /examples/example01_hs_iron(III)_nitroxide/config_ex01_3.cfg (simulates the dipolar spectrum for xi = 0°)

4) The results of the DipFit analysis will be saved in the folder:
   .../DipFit/examples/example01_hs_iron(III)_nitroxide/
