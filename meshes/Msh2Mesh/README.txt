# May 2025

This folder contains two mesh converters that transform Gmsh .msh files into INRIA .mesh files:

MATLAB converter (converter_MATLAB):
The DEFAULT converter tested extensively and used for mesh generation.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Python converter (converter_python):
Implementation that translates the MATLAB logic into Python. Goal is to enable users without MATLAB the conversion of msh Files. 
It is currently UNDER DEVELOPMENT and should be used for testing and comparison purposes only. Users are encouraged to test it with their meshes and report any inconsistencies or issues.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


With test_compare_converter.py one can compare the generated meshes (ignoring the white spaces) included
