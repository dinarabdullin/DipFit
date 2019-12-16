NAME
--------------------------------------
Example 5


DESCRIPTION
--------------------------------------
In this example, DipFit performs the fitting of an experimental dipolar spectra of a model compound (see Pic_3) which contains two spin centers, a low-spin Fe(III) and a trityl.

The  principal g-values of the low-spin Fe(III) ion are set to gxx = 1.56, gyy = 2.28, and gzz = 2.91.

The principal g-values of the trityl center are set to 2.0032.

The parameters of the fitting include the distance between the spin centers (mean value and std) and the angles xi (mean value and std) and phi (mean value and std).


HOW TO RUN  
--------------------------------------
1) Open Terminal or Command Prompt.

2) Navigate to the root directory of the DipFit program:
   cd …/DipFit/source_code
   
3) Run the program by the following command:
   python main.py examples/example05_ls_iron(III)_trityl_porphyrin/config_ex05.cfg

4) The results of the DipFit analysis will be saved in the folder:
   .../DipFit/examples/example05_ls_iron(III)_trityl_porphyrin/
