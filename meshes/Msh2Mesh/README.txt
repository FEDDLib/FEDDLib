# May 2025 ################################################

With respect to FAIR principles:

###########################################################
Files from GMSH have to be saved in ASCII version 2 format!
###########################################################

This folder contains two mesh converters that transform Gmsh .msh files into INRIA .mesh files:

MATLAB converter (converter_MATLAB): -------------------------------------------------------------------------------------------------------------------------
The converter tested extensively and used for mesh generations in recent years (until May 2025)

Python converter (converter_python): ------------------------------------------------------------------------------------------------------------------------
Implementation that translates the MATLAB logic into Python. Goal is to enable users without MATLAB the conversion of *.msh Files into *.mesh Files. 
Users with MATLAB are encouraged to test it with their meshes and report any inconsistencies or issues with the MATLAB converter. The python script should be the new Default.

-> With the file test_compare_converter.py one can compare the generated meshes (ignoring the white spaces) between MATLAB and python 
