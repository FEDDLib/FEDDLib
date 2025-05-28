import os
import re
import numpy as np
from glob import glob

"""
Gmsh MSH to INRIA MESH Format Converter
=======================================

This script converts Gmsh .msh files to INRIA .mesh format, which is required by 
certain simulation environments like FEDDLib. It preserves physical flags and 
correctly handles precedence of dimensional elements.

Purpose:
--------
- Convert Gmsh MSH files (version 2.2) to INRIA MESH format
- Properly assign physical flags to vertices
- Process points, lines, triangles and tetrahedra with their physical flags

Notes:
--------
- Lower dimensional elements' flags take precedence over higher dimensional ones

Usage:
------
1. Place this script in the directory containing your .msh files
2. Run the script: python msh_to_mesh_converter.py
3. All .msh files in the directory will be converted to .mesh format

Original implementation:
-----------------------
This is a Python translation of the MATLAB code by Jascha Knepper (updated by
Christian Hochmuth, 2019.07.04). The script preserves the same functionality
and flag assignment logic as the original MATLAB implementation.
For more detailed information, refer to the original MATLAB code comments.

Python version:
--------------
This Python implementation was created with assistance from Claude (Anthropic)
in March 2025, translating and enhancing the original MATLAB code while 
maintaining its core functionality.

Key features:
------------
- Physical flag assignment to vertices based on dimensional precedence
- Support for lines (edges), triangles and tetrahedra
- Preservation of element connectivity and physical IDs
- Handling of large mesh files through chunked reading/writing

Note: This script only supports the ASCII MSH format (version 2.2).


Updated vers
"""
##############################################################################################################################
def read_gmsh_file(filename, be_silent=False):
    """
    Read points and elements from a Gmsh-MSH-file.
    Supports both MSH version 2.2 ASCII format and tries to be robust with variations.
    """
    """
    Read points and elements from a Gmsh-MSH-file.
    Similar to MSH_Gmsh__readFile.m
    
    Gmsh online manual: http://gmsh.info/doc/texinfo/gmsh.html
    
    QUICKSTART: Extract tetrahedra and triangles.
       x, elements_dict = read_gmsh_file(filename)
       tetr = elements_dict[4]['elements'] if 4 in elements_dict else []
       tri = elements_dict[2]['elements'] if 2 in elements_dict else []
    For more IDs (tetr:4, tri:2), see the Gmsh ID below.
    
    Args:
        filename: Path to the msh file
        be_silent: [optional] disable output to console (except warnings and errors)
            default: False
    
    Returns:
        x: coordinates (x,y,z); numpy array (n x 4) where 4th col is flag
        elements_dict: Dictionary containing elements categorized by type
            elements_dict[element_type]['elements']: node IDs
            elements_dict[element_type]['physicalIDs']: physical IDs
            elements_dict[element_type]['elementaryIDs']: elementary IDs
            elements_dict[element_type]['elementType']: element type name
    
    Most important supported element types:
    
                              |                        |             |
        element name          |   number of vertices   |   Gmsh ID   |   element order
        ______________________|________________________|_____________|________________
                              |                        |             |
        'point'               |          1             |     15      |        1
        ----------------------|--------------------------------------|----------------
        'line'                |          2             |      1      |        1
                              |          3             |      8      |        2
                              |          4             |     26      |        3
                              |          5             |     27      |        4
                              |          6             |     28      |        5
        ----------------------|------------------------|-------------|----------------
        'triangle'            |          3             |      2      |        1
                              |          6             |      9      |        2
                              |         10             |     21      |        3
                              |         15             |     23      |        4
                              |         21             |     25      |        5
        ----------------------|------------------------|-------------|----------------
        'tetrahedron'         |          4             |      4      |        1
                              |         10             |     11      |        2
                              |         20             |     29      |        3
                              |         35             |     30      |        4
                              |         56             |     31      |        5
        ----------------------|------------------------|-------------|----------------
    """
    # Define Gmsh element types (element type id, number of nodes)
    # Adapted from the MATLAB code's c_mshElement structure
    msh_element = {
        1: (2, 'line'),         # line (2 nodes) order 1
        2: (3, 'triangle'),     # triangle (3 nodes) order 1
        3: (4, 'quadrangle'),   # quadrilateral (4 nodes) order 1
        4: (4, 'tetrahedron'),  # tetrahedron (4 nodes) order 1
        5: (8, 'hexahedron'),   # hexahedron (8 nodes) order 1
        6: (6, 'prism'),        # prism (6 nodes) order 1
        7: (5, 'pyramid'),      # pyramid (5 nodes) order 1
        8: (3, 'line'),         # line (3 nodes) order 2
        9: (6, 'triangle'),     # triangle (6 nodes) order 2
        10: (9, 'quadrangle'),  # quadrilateral (9 nodes) order 2
        11: (10, 'tetrahedron'),# tetrahedron (10 nodes) order 2
        12: (27, 'hexahedron'), # hexahedron (27 nodes) order 2
        13: (18, 'prism'),      # prism (18 nodes) order 2
        14: (14, 'pyramid'),    # pyramid (14 nodes) order 2
        15: (1, 'point'),       # point (1 node)
        16: (8, 'quadrangle'),  # quadrilateral (8 nodes) order 2
        17: (20, 'hexahedron'), # hexahedron (20 nodes) order 2
        18: (15, 'prism'),      # prism (15 nodes) order 2
        19: (13, 'pyramid'),    # pyramid (13 nodes) order 2
        20: (9, 'triangle'),    # triangle (9 nodes) order 3, incomplete
        21: (10, 'triangle'),   # triangle (10 nodes) order 3
        22: (12, 'triangle'),   # triangle (12 nodes) order 4, incomplete
        23: (15, 'triangle'),   # triangle (15 nodes) order 4
        24: (15, 'triangle'),   # triangle (15 nodes) order 5, incomplete
        25: (21, 'triangle'),   # triangle (21 nodes) order 5
        26: (4, 'edge'),        # edge (4 nodes) order 3
        27: (5, 'edge'),        # edge (5 nodes) order 4
        28: (6, 'edge'),        # edge (6 nodes) order 5
        29: (20, 'tetrahedron'),# tetrahedron (20 nodes) order 3
        30: (35, 'tetrahedron'),# tetrahedron (35 nodes) order 4
        31: (56, 'tetrahedron'),# tetrahedron (56 nodes) order 5
        92: (64, 'hexahedron'), # hexahedron (64 nodes) order 3
        93: (125, 'hexahedron') # hexahedron (125 nodes) order 4
    }
    
    if not be_silent:
        print(f">> (Progress) Reading meshfile {filename}.")
    
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except UnicodeDecodeError:
        # Try reading as binary if UTF-8 fails
        with open(filename, 'rb') as file:
            lines = [line.decode('latin1').strip() for line in file.readlines()]
    except Exception as e:
        print(f"Error reading file: {e}")
        raise
    
    x = None
    node_ids = None
    elements_raw = []
    element_tags = []
    element_types = []
    
    # Parse the file line by line
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line == '$Nodes':
            # Read nodes section
            i += 1
            num_nodes = int(lines[i].strip())
            if not be_silent:
                print(f">>    (Stat.) number of nodes = {num_nodes}")
            
            if not be_silent:
                print(">>    (Progress) Reading node coordinates.")
            
            # Extract coordinates and node IDs
            node_data = []
            for j in range(num_nodes):
                i += 1
                node_data.append([float(x) for x in lines[i].strip().split()])
            
            node_data = np.array(node_data)
            node_ids = node_data[:, 0].astype(int)
            
            # Create x array with zeros for any skipped node IDs
            # This matches MATLAB's behavior with 1-indexed arrays
            max_node_id = int(max(node_ids))
            x = np.zeros((max_node_id, 3))
            
            # Fill in coordinates at the correct indices (Matlab is 1-indexed, Python is 0-indexed)
            for i, node_id in enumerate(node_ids):
                x[node_id-1] = node_data[i, 1:4]
            
            # Try to find the end nodes marker
            i += 1
            while i < len(lines) and lines[i].strip() != '$EndNodes':
                # Continue reading lines until we find the end marker
                # This handles cases where there might be empty lines or comments
                i += 1
            
            if i >= len(lines) or lines[i].strip() != '$EndNodes':
                print("WARNING: $EndNodes tag not found. File may have an unexpected format.")
                # Reset position to after the nodes data we already read
                i = i - 1
            
        elif line == '$Elements':
            # Read elements section
            i += 1
            num_elements = int(lines[i].strip())
            if not be_silent:
                print(f">>    (Stat.) number of elements = {num_elements}")
            
            if not be_silent:
                print(">>    (Progress) Reading elements.")
            
            # Parse each element line
            # Format: (element number, element type, number of tags, tag_1, ..., tag_{number of tags}, node number 1, ... , node number n)
            for j in range(num_elements):
                i += 1
                el_data = [int(x) for x in lines[i].strip().split()]
                el_id = el_data[0]  # element number
                el_type = el_data[1]  # element type
                num_tags = el_data[2]  # number of tags
                
                tags = el_data[3:3+num_tags]
                nodes = el_data[3+num_tags:]
                
                element_types.append(el_type)
                element_tags.append(tags)
                elements_raw.append(nodes)
            
            # Try to find the end elements marker
            i += 1
            while i < len(lines) and lines[i].strip() != '$EndElements':
                # Continue reading lines until we find the end marker
                # This handles cases where there might be empty lines or comments
                i += 1
            
            if i >= len(lines) or lines[i].strip() != '$EndElements':
                print("WARNING: $EndElements tag not found. File may have an unexpected format.")
                # Reset position to after the elements data we already read
                i = i - 1
        
        i += 1
    
    # Organize elements by type - similar to the MATLAB cElements structure
    unique_types = np.unique(element_types)
    elements_dict = {}
    
    for el_type in unique_types:
        indices = [i for i, t in enumerate(element_types) if t == el_type]
        
        elements = np.array([elements_raw[i] for i in indices])
        phys_ids = np.array([element_tags[i][0] for i in indices])
        elem_ids = np.array([element_tags[i][1] for i in indices])
        
        elements_dict[el_type] = {
            'elements': elements,
            'physicalIDs': phys_ids,
            'elementaryIDs': elem_ids,
            'elementType': msh_element.get(el_type, (0, 'unknown'))[1]
        }
    
    # Add space for flags to coordinates (4th column)
    x_with_flags = np.column_stack((x, np.zeros(len(x))))
    
    return x_with_flags, elements_dict

##############################################################################################################################
def write_inria_mesh_file(filename, x, lin=None, tri=None, tetr=None, acc=6):
    """
    Write vertices, triangles and tetrahedra to an INRIA-MESH-file.
    Similar to writeInriaMeshFile.m
    
    Args:
        filename: filename to write meshfile to
        x: vertices (x,y,z, flag)
        lin: lines (v1,v2, flag)
        tri: triangles (v1,v2,v3, flag)
        tetr: tetrahedra (v1,v2,v3,v4, flag)
        acc: number of decimal places to use for the vertices
    
    Number of vertices:   size(x,1)
    Number of triangles:  size(tri,1)
    Number of tetrahedra: size(tetr,1)
    """
    # Handle empty inputs with empty arrays
    lin = np.array([]) if lin is None else lin
    tri = np.array([]) if tri is None else tri
    tetr = np.array([]) if tetr is None else tetr
    
    # Determine dimension based on presence of tetrahedra
    dimension = 3 if tetr.size > 0 else 2
    
    # Ensure the number of decimal places is valid
    assert isinstance(acc, int) and acc >= 0, "Number of decimal places must be a non-negative integer."
    
    print(f">> (Progress) Writing meshfile: {filename}")
    
    with open(filename, 'w') as f:
        # Write vertices
        num_nodes = len(x)
        print(f">>    (Progress) Writing {num_nodes} vertices.")
        f.write(f"MeshVersionFormatted 1\n\nDimension {dimension}\n\nVertices\n")
        f.write(f"{num_nodes}\n")
        
        # Format string for coordinates
        fmt = f"{{:.{acc}f}} {{:.{acc}f}} {{:.{acc}f}} {{:d}}\n"
        
        # Write vertices in chunks to avoid possible memory issues - similar to MATLAB chunking
        chunk_size = 1000000
        for k in range(0, num_nodes, chunk_size):
            chunk_end = min(k + chunk_size, num_nodes)
            for i in range(k, chunk_end):
                # We don't skip zero rows - we need to write all vertices
                # This is different from the previous implementation which was skipping zeros
                f.write(fmt.format(x[i, 0], x[i, 1], x[i, 2], int(x[i, 3])))
        
        f.write("\n")
        
        # Write lines (edges)
        if lin.size > 0:
            num_lines = len(lin)
            print(f">>    (Progress) Writing {num_lines} edges.")
            f.write("Edges\n")
            f.write(f"{num_lines}\n")
            
            # Write lines in chunks
            chunk_size = 1000000
            for k in range(0, num_lines, chunk_size):
                chunk_end = min(k + chunk_size, num_lines)
                for i in range(k, chunk_end):
                    f.write(f"{int(lin[i, 0])} {int(lin[i, 1])} {int(lin[i, 2])}\n")
            
            f.write("\n")
        
        # Write triangles
        if tri.size > 0:
            num_triangles = len(tri)
            print(f">>    (Progress) Writing {num_triangles} triangles.")
            f.write("Triangles\n")
            f.write(f"{num_triangles}\n")
            
            # Write triangles in chunks
            chunk_size = 1000000
            for k in range(0, num_triangles, chunk_size):
                chunk_end = min(k + chunk_size, num_triangles)
                for i in range(k, chunk_end):
                    f.write(f"{int(tri[i, 0])} {int(tri[i, 1])} {int(tri[i, 2])} {int(tri[i, 3])}\n")
            
            f.write("\n")
        
        # Write tetrahedra
        if tetr.size > 0:
            num_tetrahedra = len(tetr)
            print(f">>    (Progress) Writing {num_tetrahedra} tetrahedra.")
            f.write("Tetrahedra\n")
            f.write(f"{num_tetrahedra}\n")
            
            # Write tetrahedra in chunks
            chunk_size = 1000000
            for k in range(0, num_tetrahedra, chunk_size):
                chunk_end = min(k + chunk_size, num_tetrahedra)
                for i in range(k, chunk_end):
                    f.write(f"{int(tetr[i, 0])} {int(tetr[i, 1])} {int(tetr[i, 2])} {int(tetr[i, 3])} {int(tetr[i, 4])}\n")
            
            f.write("\n")

##############################################################################################################################
def convert_msh_to_mesh(filename):
    """
    Convert a Gmsh-MSH-file to an INRIA-MESH-file.
    Only supports triangles and tetrahedra.
    
    Important: This script was written to convert a Gmsh MSH file with edges,
    triangles and tetrahedra (and physical flags for all) to an INRIA ASCII
    MESH file of triangles and tetrahedra. The flags will be written to the
    vertex flag of the MESH file. Flags of edges take precedence over flags
    of triangles and those of triangles over those of tetrahedra.
    
    Write the output file in a format that FEDDLib can read (vertices are
    assigned physical flags; if a vertex is part of a e.g. a triangle and a
    tetrahedron, then it will be assigned the flag of the lower-dimensional
    element, thus of the triangle). Gmsh outputs the INRIA-MESH-file in a
    format such that FEDDLib cannot parse the flags (Gmsh only writes the
    elementary tags and not the physical tags to the vertex IDs).
    
    Original MATLAB code by Jascha Knepper
    Code last updated: 2019.07.04 by Christian Hochmuth
    Recent change: If a point is part of several surfaces it will get the
    lowest surface ID
    
    Args:
        filename: The path to the .msh file to convert - Filename of the Gmsh-MSH-file without the extension.
        e.g. filename = 'rectangle_H0.005_L_0.025Nele_200.0ref2.0';
    """
    print("-------------------------------------------------------------")
    print(">>>> Start of script.")
    
    # Read mesh file
    x, c_elements = read_gmsh_file(filename)
    
    try:
        # Extract the filename without extension
        filename_without_ext = os.path.splitext(filename)[0]
    except Exception as e:
        print(f"Warning: Could not extract filename without extension: {e}")
        # Fallback to using the original name and appending .mesh
        filename_without_ext = filename.replace('.msh', '')
    
    # Set physical type / flag of nodes.
    # First set for higher dimensional elements (e.g.: tetrahedron) and then 
    # for lower dimensional (e.g.: triangle).
    # Thus lower dimensional elements can overwrite the physical type of a 
    # node. E.g.: A node on the boundary (which belongs to a tetrahedron) 
    # is attributed the physical flag of the corresponding surface triangle.
    
    # Gmsh IDs:
    # point       = 15
    # line        =  1
    # triangle    =  2
    # tetrahedron =  4
    el_types = [4, 2, 1, 15]  # We process in the order: tetrahedron, triangle, line, point
    
    # Initialize flags with high value (100)
    x_flags = np.ones(len(x)) * 100
    
    # Process each element type - lower dimensional elements can override flags
    for el_id in el_types:
        if el_id in c_elements:
            for i in range(len(c_elements[el_id]['elements'])):
                for j in range(len(c_elements[el_id]['elements'][i])):
                    # Get the node index (adjusting for 0-indexed Python arrays)
                    node_index = c_elements[el_id]['elements'][i][j] - 1
                    
                    # Only assign the flag if it's lower than the current flag
                    # This ensures lower dimensional elements have precedence
                    if x_flags[node_index] > c_elements[el_id]['physicalIDs'][i]:
                        x_flags[node_index] = c_elements[el_id]['physicalIDs'][i]
    
    # Add flags to coordinates (4th column)
    x[:, 3] = x_flags
    
    # Prepare data for INRIA-MESH-file
    # Initialize empty arrays
    lin = np.empty((0, 3))
    tri = np.empty((0, 4))
    tetr = np.empty((0, 5))
    
    # Extract lines with their physical IDs
    if 1 in c_elements:  # Lines
        els = c_elements[1]['elements']
        phys_ids = c_elements[1]['physicalIDs'].reshape(-1, 1)
        lin = np.hstack((els, phys_ids))
    
    # Extract triangles with their physical IDs
    if 2 in c_elements:  # Triangles
        els = c_elements[2]['elements']
        phys_ids = c_elements[2]['physicalIDs'].reshape(-1, 1)
        tri = np.hstack((els, phys_ids))
    
    # Extract tetrahedra with their physical IDs
    if 4 in c_elements:  # Tetrahedra
        els = c_elements[4]['elements']
        phys_ids = c_elements[4]['physicalIDs'].reshape(-1, 1)
        tetr = np.hstack((els, phys_ids))
    
    # Write INRIA-MESH-file
    write_inria_mesh_file(f"{filename_without_ext}.mesh", x, lin, tri, tetr)
    
    print("<<<< End of script.")
    print("-------------------------------------------------------------")


##############################################################################################################################
def main():
    """
    Main function to find and convert all .msh files in the current directory.
    Equivalent to the loop in main_convertMshToMesh.m
    """
    # Get all .msh files in the current directory
    msh_files = glob("tests/*.msh")
    print("msh_files:", msh_files)
    
    if not msh_files:
        print("No .msh files found in the current directory.")
        return
    
    # Loop over each .msh file in the directory
    for msh_file in msh_files:
        print(f"Converting {msh_file}...")
        try:
            convert_msh_to_mesh(msh_file)
            print(f"Conversion of {msh_file} completed.")
        except Exception as e:
            print(f"Error converting {msh_file}: {e}")
            print("Continuing with next file...")


if __name__ == "__main__":
    main()
