import os
import difflib

"""
Mesh File Comparison Script
===========================

Date: 2025-05-26
Generated with the help of GPT-4o-mini (OpenAI)

Description:
------------
This script compares all `.mesh` files with identical names located in two separate 
directories: `converter_MATLAB/tests` and `converter_python/tests`.

For each pair of matching files, it checks if the files are identical, ignoring 
differences caused by trailing whitespace and blank lines at the end of the files.
"""


def compare_mesh_files(file1, file2):
    # Compare ignoring whitespace and trailing newlines
    def remove_whitespaces(file_path):
        with open(file_path, 'r') as f:
            lines = f.readlines()
        # Strip trailing whitespace and empty lines at end
        lines = [line.rstrip() + '\n' for line in lines]
        while lines and lines[-1].strip() == '':
            lines.pop()
        return lines
    
    
    modified_file1 = remove_whitespaces(file1)
    modified_file2 = remove_whitespaces(file2)
    #modified_file1=file1
    #modified_file2=file2

    if modified_file1 == modified_file2:
        print(f"[PASS] Files are identical (ignoring trailing whitespace): {file1} and {file2}")
        return True
    else:
        print(f"[FAIL] Files differ: {file1} and {file2}")
        diff = difflib.unified_diff(modified_file1, modified_file2, fromfile=file1, tofile=file2)
        print(''.join(diff))
        return False

if __name__ == "__main__":
    matlab_folder = 'converter_MATLAB/tests'
    python_folder = 'converter_python/tests'

    # Get set of .mesh filenames in both folders
    matlab_files = {f for f in os.listdir(matlab_folder) if f.endswith('.mesh')}
    python_files = {f for f in os.listdir(python_folder) if f.endswith('.mesh')}

    # Find common files by name
    common_files = matlab_files.intersection(python_files)

    if not common_files:
        print("No matching .mesh files found in both folders.")
    else:
        for filename in sorted(common_files):
            matlab_path = os.path.join(matlab_folder, filename)
            python_path = os.path.join(python_folder, filename)
            compare_mesh_files(matlab_path, python_path)

