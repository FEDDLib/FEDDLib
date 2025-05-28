% Convert a Gmsh-MSH-file to an INRIA-MESH-file.
% Only supports triangles and tetrahedra.
%
% Important: This script was written to convert a Gmsh MSH file with edges,
% triangles and tetrahedra (and physical flags for all) to an INRIA ASCII
% MESH file of triangles and tetrahedra. The flags will be written to the
% vertex flag of the MESH file. Flags of edges take precedence over flags
% of triangles and those of triangles over those of tetrahedra.
%
% Write the output file in a format that FEDDLib can read (vertices are
% assigned physical flags; if a vertex is part of a e.g. a triangle and a
% tetrahedron, then it will be assigned the flag of the lower-dimensional
% element, thus of the triangle). Gmsh outputs the INRIA-MESH-file in a
% format such that FEDDLib cannot parse the flags (Gmsh only writes the
% elementary tags and not the physical tags to the vertex IDs).
%
% Code by Jascha Knepper
% Code last updated: 2019.07.04
% by Christian Hochmuth
% recenct change: If a point is part of several surfaces it will get the
% lowest surface ID
clear;clc
fprintf('-------------------------------------------------------------\n')
fprintf('>>>> Start of script.\n')

%% User defined settings
% Filename of the Gmsh-MSH-file without the extension.
filename = 'tests/square';

%% read meshfile
[x,~,cElements] = MSH_Gmsh__readFile([filename '.msh']);

%% Convert Gmsh-MSH format to INRIA-MESH format.

% Set physical type / flag of nodes.
% First set for higher dimensional elements (e.g.: tetrahedron) and then 
% for lower dimensional (e.g.: triangle).
% Thus lower dimensional elements can overwrite the physical type of a 
% node. E.g.: A node on the boundary (which belongs to a tetrahedron) 
% is attributed the physical flag of the corresponding surface triangle.
% Gmsh IDs:
% point       = 15
% line        =  1
% triangle    =  2
% tetrahedron =  4
elTypes = [4,2,1,15];
x_flags = 100*ones(size(x,1),1);  % initialize flags

for k = 1:length(elTypes)
    elID = elTypes(k);
    if not(isempty(cElements{elID})) %&& elID~=1 %we do an etxra sweep for the line elements below
        for i = 1:size(cElements{elID}.elements,1)
            for j=1:size(cElements{elID}.elements,2)                              
                if x_flags(cElements{elID}.elements(i,j)) > cElements{elID}.physicalIDs(i) 
                    x_flags(cElements{elID}.elements(i,j)) = cElements{elID}.physicalIDs(i);                                                   
                end     
            end
        end
    end
end
%We do an extra sweep for the line elements, which we want to prefer.
%Therefore we first set an arbitrary high value again
% elID = 1;
% if not(isempty(cElements{elID}))
%     for i = 1:size(cElements{elID}.elements,1)
%         for j=1:size(cElements{elID}.elements,2)                              
%             x_flags(cElements{elID}.elements(i,j)) = 100;           
%         end
%     end
% end
% if not(isempty(cElements{elID}))
%     for i = 1:size(cElements{elID}.elements,1)
%         for j=1:size(cElements{elID}.elements,2)                                          
%             if x_flags(cElements{elID}.elements(i,j)) > cElements{elID}.physicalIDs(i) 
%                 x_flags(cElements{elID}.elements(i,j)) = cElements{elID}.physicalIDs(i);                                                   
%             end     
%         end
%     end
% end

x(:,4) = x_flags;
if not(isempty(cElements{1}))
    lin  = [cElements{1}.elements,cElements{1}.physicalIDs];
else
    lin = [];
end
if not(isempty(cElements{2}))
    tri  = [cElements{2}.elements,cElements{2}.physicalIDs];
else
    tri = [];
end
%cElements{2} = [];
if not(isempty(cElements{4}))
    tetr = [cElements{4}.elements,cElements{4}.physicalIDs];
else 
    tetr =[];
end

%clear('cElements')

%% Write INRIA-MESH-file.
writeInriaMeshFile([filename '.mesh'],x,lin,tri,tetr);

%% Done.
fprintf('<<<< End of script.\n')
fprintf('-------------------------------------------------------------\n')
