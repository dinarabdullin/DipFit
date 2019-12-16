NAME
--------------------------------------
Example 4


DESCRIPTION
--------------------------------------
In this example, the dipolar spectra is simulated for a spin system consisting of a low-spin Fe(III) ion and a nitroxide. 

The  principal g-values of the low-spin Fe(III) ion are set to gxx = 1.56, gyy = 2.28, and gzz = 2.91.

For simplicity, the principal g-values of the nitroxide spin center are set to the free electron g-factor.

The distance between the spin centers is set to 2.5 nm. 
The orientation of the distance vector relative to the g-frame of the low-spin Fe(III) ion is described by two angles, xi - the polar angle and phi - the azimuthal angle.


HOW TO RUN  
--------------------------------------
1) Open Terminal or Command Prompt.

2) Navigate to the root directory of the DipFit program:
   cd …/DipFit/source_code
   
3) Run the program by the following command:
   python main.py examples/example04_ls_iron(III)_nitroxide/config_ex04_1.cfg (simulates the dipolar spectrum vs angle xi, phi = 0°)
   python main.py examples/example04_ls_iron(III)_nitroxide/config_ex04_2.cfg (simulates the dipolar spectrum vs angle xi, phi = 30°)
   python main.py examples/example04_ls_iron(III)_nitroxide/config_ex04_3.cfg (simulates the dipolar spectrum vs angle xi, phi = 60°)
   python main.py examples/example04_ls_iron(III)_nitroxide/config_ex04_3.cfg (simulates the dipolar spectrum vs angle xi, phi = 90°)

4) The results of the DipFit analysis will be saved in the folder:
   .../DipFit/examples/example04_ls_iron(III)_nitroxide/
