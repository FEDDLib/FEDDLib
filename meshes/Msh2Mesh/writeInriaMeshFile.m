function writeInriaMeshFile(filename,x,lin,tri,tetr,acc)
% Write vertices, triangles and tetrahedra to an INRIA-MESH-file.
% Code last updated: 2019.07.04
% 
% INPUT:
%   filename: filename to write meshfile to.
%   x:    vertices   (x,y,z,       flag)
%   lin:  lines      (v1,v2,       flag)
%   tri:  triangles  (v1,v2,v3,    flag)
%   tetr: tetrahedra (v1,v2,v3,v4, flag)
%   acc:  number of decimal places to use for the vertices
% Number of vertices:   size(x,1)
% Number of triangles:  size(tri,1)
% Number of tetrahedra: size(tetr,1)
    %
    if nargin < 6
        acc = 6;
    else
        assert(round(acc)-acc == 0,'Number of decimal places must be an integer.')
        assert(acc >= 0,'Number of decimal places must not be negative.')
    end
    dimension = 3;
    if size(tetr,1)==0
        dimension = 2;
    end
    
    fprintf('>> (Progress) Writing meshfile: %s\n',filename)
    
    % open file for writing
    fid = fopen(filename,'w');
    
    % write vertices
    numberOfNodes = size(x,1);  
    fprintf('>>    (Progress) Writing %i vertices.\n',numberOfNodes)
    fprintf(fid,'MeshVersionFormatted 1\n\nDimension %i\n\nVertices\n',dimension);
    fprintf(fid,'%d\n',numberOfNodes);
    str = ['%.',num2str(acc),'f '];
    %
    % Write in chunkgs to avoid possible 32-bit limits on the length of
    % strings.
    m = ceil(numberOfNodes/1000000);
    for k = 1:m
        m1 = (k-1)*1000000+1;
        m2 = min(k*1000000,numberOfNodes);
        oneLine = sprintf([str,str,str,'%d\n'],x(m1:m2,:)');        
        fprintf(fid,'%s',oneLine);
    end
    fprintf(fid,'\n');
    
     % write lines
    numberOfLines= size(lin,1);
    if numberOfLines>0
        fprintf('>>    (Progress) Writing %i edges.\n',numberOfLines)
        fprintf(fid,'Edges\n');
        fprintf(fid,'%d\n',numberOfLines);
        %
        % Write in chunkgs to avoid possible 32-bit limits on the length of
        % strings.
        m = ceil(numberOfLines/1000000);
        for k = 1:m            
            m1 = (k-1)*1000000+1;
            m2 = min(k*1000000,numberOfLines);
            oneLine = sprintf('%d %d %d \n',lin(m1:m2,:)');    
            fprintf(fid,'%s',oneLine);
        end
        fprintf(fid,'\n');
    end
    % write triangles
    numberOfTriangles = size(tri,1);
    if numberOfTriangles>0
        fprintf('>>    (Progress) Writing %i triangles.\n',numberOfTriangles)
        fprintf(fid,'Triangles\n');
        fprintf(fid,'%d\n',numberOfTriangles);
        %
        % Write in chunkgs to avoid possible 32-bit limits on the length of
        % strings.
        m = ceil(numberOfTriangles/1000000);
        for k = 1:m        
            m1 = (k-1)*1000000+1;
            m2 = min(k*1000000,numberOfTriangles);
            oneLine = sprintf('%d %d %d %d\n',tri(m1:m2,:)');
            fprintf(fid,'%s',oneLine);
        end
        fprintf(fid,'\n');
    end
    
    % write tetrahedra
    numberOfTetrahedra = size(tetr,1);   
    if numberOfTetrahedra>0
        fprintf('>>    (Progress) Writing %i tetrahedra.\n',numberOfTetrahedra)
        fprintf(fid,'Tetrahedra\n');
        fprintf(fid,'%d\n',numberOfTetrahedra);
        %
        % Write in chunkgs to avoid possible 32-bit limits on the length of
        % strings.
        m = ceil(numberOfTetrahedra/1000000);
        for k = 1:m        
            m1 = (k-1)*1000000+1;
            m2 = min(k*1000000,numberOfTetrahedra);
            oneLine = sprintf('%d %d %d %d %d\n',tetr(m1:m2,:)');            
            fprintf(fid,'%s',oneLine);
        end
        fprintf(fid,'\n');
    end
    % close file
    fclose(fid);
end
