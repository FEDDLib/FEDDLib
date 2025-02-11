function [x,elementTypes,cElements,cElementsRaw,cElementTags,elementTypesRaw,nodeIDs,c_mshElement] = MSH_Gmsh__readFile( ...
        filename,disableClustering,beSilent)
% Read points and elements from a Gmsh-MSH-file.
% Code last updated: 2018.10.18
% Gmsh online manual: http://gmsh.info/doc/texinfo/gmsh.html
%
% QUICKSTART: Extract tetrahedra and triangles.
%    [x,cElements,~,elementTypes] = MSH_Gmsh__readFile(filename);
%    tetr = cell2mat(cElements(elementTypes==4));
%    tri = cell2mat(cElements(elementTypes==2));
% For more IDs (tetr:4, tri:2), see the Gmsh ID below.
%
% INPUT:
%   filename: msh file
%   beSilent: [optional] disable output to console (except for warnings 
%      and errors).
%      default: false
% 
% OUTPUT:
%   x: coordinates (x,y,z); matrix (n x 3)
%   cElements: (v1,v2,...,vn) node IDs
%       elements can be points(v1), lines(v1,v2), triangles(v1,v2,v3), 
%       tetrahedra(v1,v2,v3,v3) etc., see list below.
%       cell array
%   cElementTags: (t1,t2,...); t1:physical id; t2:elementary id;
%       number of tags is variable;
%       cell array
%   elementTypes: vector with Gmsh-IDs for {point,line,triangle,...}
%       see list below
%   nodeIDs: Generally, 'x' will have the IDs 1,2,3,... . However, it is
%       permitted in the Gmsh format that IDs are skipped.
%       In that case 'x' will contain rows of zeros and 'nodeIDs' will
%       contain the row numbers with data.
%
% Number of nodes:    size(x,1)
% Number of elements: length(cElements) = length(cElementTags)
%                                       = length(elementTypes)
%
%
% Most important supported element types:
% 
%                          |                        |             |
%    element name          |   number of vertices   |   Gmsh ID   |   element order
%    ______________________|________________________|_____________|________________
%                          |                        |             |
%    'point'               |          1             |     15      |        1
%    ----------------------|--------------------------------------|----------------
%    'line'                |          2             |      1      |        1
%                          |          3             |      8      |        2
%                          |          4             |     26      |        3
%                          |          5             |     27      |        4
%                          |          6             |     28      |        5
%    ----------------------|------------------------|-------------|----------------
%    'triangle'            |          3             |      2      |        1
%                          |          6             |      9      |        2
%                          |         10             |     21      |        3
%                          |         15             |     23      |        4
%                          |         21             |     25      |        5
%    ----------------------|------------------------|-------------|----------------
%    'incomplete triangle' |          9             |     20      |        3
%                          |         12             |     22      |        4
%                          |         15             |     24      |        5
%    ----------------------|------------------------|-------------|----------------
%    'rectangle'           |          4             |      3      |        1
%    ----------------------|------------------------|-------------|----------------
%    'quadrangle'          |          4             |      3      |        1
%                          |          9             |     10      |        2
%                          |          8             |     16      |        2
%    ----------------------|------------------------|-------------|----------------
%    'tetrahedron'         |          4             |      4      |        1
%                          |         10             |     11      |        2
%                          |         20             |     29      |        3
%                          |         35             |     30      |        4
%                          |         56             |     31      |        5
%    ----------------------|------------------------|-------------|----------------
%    'cuboid'              |          8             |      5      |        1
%    ----------------------|------------------------|-------------|----------------
%    'hexahedron'          |          8             |      5      |        1
%                          |         27             |     12      |        2
%                          |         20             |     17      |        2
%                          |         64             |     92      |        3
%                          |        125             |     93      |        4
%    ----------------------|------------------------|-------------|----------------
%    'prism'               |          6             |      6      |        1
%                          |         18             |     13      |        2
%                          |         15             |     18      |        2
%    ----------------------|------------------------|-------------|----------------
%    'pyramid'             |          5             |      7      |        1
%                          |         14             |     14      |        2
%                          |         13             |     19      |        2
%    ______________________|________________________|_____________|________________
%
%
% TODO: Binary format
%   Even though the Gmsh binary format allows for efficient reading and 
%   writing of binary files, the tool Gmsh itself (v3) does not make use 
%   of this format properly which makes reading the files inefficient: 
%   The binary file format permits to store a header for each element. 
%   The header stores the element type and the number of tags the element has. 
%   As, almost always, all elements will have a fixed number of tags, the 
%   file format also permits to store a header for a block of elements of 
%   the same type. Therefore, most of the time, you can store and read, e.g., 
%   all triangles in one block. Unfortunately, Gmsh itself does not make use 
%   of this feature but stores an element header for each element which 
%   makes the files larger and less efficient to read. 
%   Furthermore, it is much more difficult to read these files somewhat 
%   efficiently. For this, you need to read in chunks from the file and then 
%   analyze them. \n
%   Currently, this somewhat efficient reading method is not implemented, 
%   so reading binary MSH-files, exported by Gmsh, is not recommended 
%   (and very slow).
%
    %
    
    if (exist('beSilent','var') ~= 1) || isempty(beSilent)
        beSilent = false;
    end
    
    if (exist('disableClustering','var') ~= 1) || isempty(disableClustering)
        disableClustering = false;
    end
    
    
    % Gmsh element (element type id, number of nodes)
    % http://geuz.org/gmsh/doc/texinfo/gmsh.html#MSH-ASCII-file-format
    c_mshElement = {
         1,  2, 'line';        % line          (  2 nodes) order 1
         2,  3, 'triangle';    % triangle      (  3 nodes) order 1
         3,  4, 'quadrangle';  % quadrilateral (  4 nodes) order 1
         4,  4, 'tetrahedron'; % tetrahedron   (  4 nodes) order 1
         5,  8, 'hexahedron';  % hexahedron    (  8 nodes) order 1
         6,  6, 'prism';       % prism         (  6 nodes) order 1
         7,  5, 'pyramid';     % pyramid       (  5 nodes) order 1
         8,  3, 'line';        % line          (  3 nodes) order 2
         9,  6, 'triangle';    % triangle      (  6 nodes) order 2
        10,  9, 'quadrangle';  % quadrilateral (  9 nodes) order 2
        11, 10, 'tetrahedron'; % tetrahedron   ( 10 nodes) order 2
        12, 27, 'hexahedron';  % hexahedron    ( 27 nodes) order 2
        13, 18, 'prism';       % prism         ( 18 nodes) order 2
        14, 14, 'pyramid';     % pyramid       ( 14 nodes) order 2
        15,  1, 'point';       % point         (  1 node )
        16,  8, 'quadrangle';  % quadrilateral (  8 nodes) order 2
        17, 20, 'hexahedron';  % hexahedron    ( 20 nodes) order 2
        18, 15, 'prism';       % prism         ( 15 nodes) order 2
        19, 13, 'pyramid';     % pyramid       ( 13 nodes) order 2
        20,  9, 'triangle';    % triangle      (  9 nodes) order 3, incomplete
        21, 10, 'triangle';    % triangle      ( 10 nodes) order 3
        22, 12, 'triangle';    % triangle      ( 12 nodes) order 4, incomplete
        23, 15, 'triangle';    % triangle      ( 15 nodes) order 4
        24, 15, 'triangle';    % triangle      ( 15 nodes) order 5, incomplete
        25, 21, 'triangle';    % triangle      ( 21 nodes) order 5
        26,  4, 'edge';        % edge          (  4 nodes) order 3
        27,  5, 'edge';        % edge          (  5 nodes) order 4
        28,  6, 'edge';        % edge          (  6 nodes) order 5
        29, 20, 'tetrahedron'; % tetrahedron   ( 20 nodes) order 3
        30, 35, 'tetrahedron'; % tetrahedron   ( 35 nodes) order 4
        31, 56, 'tetrahedron'; % tetrahedron   ( 56 nodes) order 5
        92, 64, 'hexahedron';  % hexahedron    ( 64 nodes) order 3
        93,125, 'hexahedron'}; % hexahedron    (125 nodes) order 4
    ids = cell2mat(c_mshElement(:,1));
    mshElement(ids,1:2) = cell2mat(c_mshElement(:,1:2));
    c_mshElement(ids,:) = c_mshElement;
    
    if not(beSilent)
        fprintf('>> (Progress) Reading meshfile %s.\n',filename)
    end
    
    [sFormat,~,fileID] = checkHeader(filename,beSilent);
    
    if strcmpi(sFormat,'ascii')
        [x,cElementsRaw,cElementTags,elementTypesRaw,nodeIDs] = readASCII(fileID,beSilent,mshElement);
    elseif strcmpi(sFormat,'binary')
        [x,cElementsRaw,cElementTags,elementTypesRaw,nodeIDs] = readBinary(fileID,beSilent,mshElement);
    end
    
    if not(disableClustering)
        elementTypes = unique(elementTypesRaw);
        cElements = cell(max(mshElement(:,1)),1);
        for i = 1:length(elementTypes)
            id = elementTypes(i);
            ind = (elementTypesRaw == id);
            tags = cell2mat(cElementTags(ind));
            cElements{id}.elements = cell2mat(cElementsRaw(ind));
            cElements{id}.physicalIDs = tags(:,1);
            cElements{id}.elementaryIDs = tags(:,2);
            cElements{id}.elementType = c_mshElement{id,3};
        end
    else
        elementTypes = [];
        cElements = [];
    end
end

function [x,cElements,cElementTags,elementTypes,nodeIDs] = readBinary(fileID,beSilent,mshElement)
    bFoundNodes = false;
    bFoundElements = false;
    label_Nodes = '$Nodes';
    label_Elements = '$Elements';
    while not(feof(fileID))
        line = fgetl(fileID);
        if isempty(line)
            % ignore
        elseif strcmpi(line,label_Nodes)
            [x,nodeIDs] = readNodes_binary(fileID,beSilent);
            bFoundNodes = true;
        elseif strcmpi(line,label_Elements)
            [cElements,cElementTags,elementTypes] = readElements_binary(fileID,beSilent,mshElement);
            bFoundElements = true;
        else
            if not(beSilent)
                fprintf('>> (Warning) Ignoring line %s.\n',line)
            end
        end
    end
    
    assert(bFoundNodes,'Nodes were not found.')
    assert(bFoundElements,'Elements were not found.')
end

function [cElements,cElementTags,elementTypes] = readElements_binary(fileID,beSilent,mshElement)
    
    % number of elements
    numberOfElements = fscanf(fileID,'%d\n',1);
    if not(beSilent)
        fprintf('>>    (Stat.) number of elements = %i\n',numberOfElements)
    end
    
    % read in elements
    % (element number, element type, number of tags, tag_1, ..., tag_{number of tags}, 
    %   node number 1, ... , node number n)
    if not(beSilent)
        fprintf('>>    (Progress) Reading elements.\n');
    end
    
    cElements = cell(numberOfElements,1);
    cElementTags = cell(numberOfElements,1);
    elementTypes = zeros(numberOfElements,1);
    
    cntEl = 0;
    while cntEl < numberOfElements
        % Read element header.
        header_el = fread(fileID,3,'int32');
        elementType = header_el(1);
        nel = header_el(2);     % number of elements of that type
        nTags = header_el(3);   % number of tags for each element
        nNodesPerElement = mshElement(elementType,2);
        sz = 1 + ...            % element ID
             nTags + ...        % number of tags
             nNodesPerElement;  % number of nodes per element
        elements = fread(fileID,[sz,nel],'int32');
        elements = elements';
        
        if nel == 1
            cElements{cntEl+1} = elements(end-nNodesPerElement+1:end);
            cElementTags{cntEl+1} = elements(2:1+nTags);
        else
            cElements(cntEl+1:cntEl+nel) = mat2cell( ...
                elements(:,end-nNodesPerElement+1:end), ...
                ones(1,nel), ...
                nNodesPerElement);
            cElementTags(cntEl+1:cntEl+nel) = mat2cell( ...
                elements(:,2:1+nTags), ...
                ones(1,nel), ...
                nTags);
        end
        elementTypes(cntEl+1:cntEl+nel) = elementType;
        
        cntEl = cntEl + nel;
    end
    fgetl(fileID);  % skip new line character
    
    if not(all(mshElement(elementTypes,1) == elementTypes))
        warning('Cannot read msh file. Unknown element type(s).')
    end
    
    lineEndElements = fgetl(fileID);
    label_Elements_End = '$EndElements';
    assert(strcmp(lineEndElements,label_Elements_End), ...
        'Unknown file format: Expected %s. Found %s.', ...
        label_Elements_End,lineEndElements)
end

function [x,nodeIDs] = readNodes_binary(fileID,beSilent)
    % number of nodes
    numberOfNodes = fscanf(fileID,'%d\n',1);
    if not(beSilent)
        fprintf('>>    (Stat.) number of nodes = %i\n',numberOfNodes)
    end
    
    % read in node coordinates (x,y,z) and node IDs
    if not(beSilent)
        fprintf('>>    (Progress) Reading node coordinates.\n')
    end
    %
    % The format is (nodeID,x     ,y     ,z     ) with the data types 
    %               (int32 ,double,double,double).
    % Since int32 is 4 bytes and double is 8 bytes, each node requires 
    % (1+3*2)*4 bytes of storages.
    % We cannot read vectors of mixed data types with Matlab.
    % Therefore we read (1+3*2) int32 for each node, extract the node ID
    % and then convert the coordinates to the double format.
    tmp = int32(zeros(numberOfNodes,1+3*2));
    tmp(:,:) = transpose(fread(fileID,[1+3*2,numberOfNodes],'int32'));
    nodeIDs = double(tmp(:,1));
    tmp = transpose(tmp(:,2:end));
    points = typecast(tmp(:),'double');
    points = transpose(reshape(points,[3,numberOfNodes]));
    x = zeros(max(nodeIDs),size(points,2));
    x(nodeIDs,:) = points;
    %
    fgetl(fileID);  % skip new line character
    
    lineEndNodes = fgetl(fileID);
    label_Nodes_End = '$EndNodes';
    assert(strcmp(lineEndNodes,label_Nodes_End), ...
        'Unknown file format: Expected %s. Found %s.', ...
        label_Nodes_End,lineEndNodes)
end

function [x,cElements,cElementTags,elementTypes,nodeIDs] = readASCII(fileID,beSilent,mshElement)
    
    % Read file as text.
    tmp = textscan(fileID,'%s','Delimiter','\n');
    lines = tmp{1};
    
    % Close file.
    fclose(fileID);
    
    i_line = 1;
    
    bFoundNodes = false;
    bFoundElements = false;
    label_Nodes = '$Nodes';
    label_Elements = '$Elements';
    while i_line < length(lines)
        if isempty(lines{i_line})
            % ignore
        elseif strcmpi(lines{i_line},label_Nodes)
            [i_line,x,nodeIDs] = readNodes_ASCII(lines,i_line+1,beSilent);
            bFoundNodes = true;
        elseif strcmpi(lines{i_line},label_Elements)
            [i_line,cElements,cElementTags,elementTypes] = readElements_ASCII(lines,i_line+1,beSilent,mshElement);
            bFoundElements = true;
        else
            if not(beSilent)
                fprintf('>> (Warning) Ignoring line %s.\n',lines{i_line})
            end
            i_line = i_line + 1;
        end
    end
    
    assert(bFoundNodes,'Nodes were not found.')
    assert(bFoundElements,'Elements were not found.')
end

function [sFormat,i_line,fileID] = checkHeader(filename,beSilent)
    % Open file.
    fileID = fopen(filename,'r');
    assert(fileID ~= -1,'Could not open file %s.',filename)
    
    i_line = 1;
    line1 = fgetl(fileID);
    meshFormat_Label = '$MeshFormat';
    assert(strcmp(line1,meshFormat_Label), ...
        'Unknown file format: line %i. Expected %s. Found %s.', ...
        i_line,meshFormat_Label,line1)
    i_line = i_line + 1;
    
    line2 = fgetl(fileID);
    meshFormat__ASCII  = '2.2 0 8';
    meshFormat__BINARY = '2.2 1 8';
    if strcmp(line2,meshFormat__ASCII)
        sFormat = 'ascii';
        if not(beSilent)
            fprintf('>>    (Info) File is in ASCII format.\n')
        end
    elseif strcmp(line2,meshFormat__BINARY)
        sFormat = 'binary';
        if not(beSilent)
            fprintf('>>    (Info) File is in BINARY format.\n')
        end
    else
        error('Unknown mesh format: line %i. Expected %s (ASCII) or %s (binary). Found %s.', ...
                i_line,meshFormat__ASCII,meshFormat__BINARY,line2)
    end
    i_line = i_line + 1;
    
    if strcmpi(sFormat,'binary')
        % little/big endian test
        one = fread(fileID,1,'int32=>double');
        assert(one == 1,'Expected to read 1, but read %g instead. Endian test failed.',one)
        fgetl(fileID);  % read new line character
        i_line = i_line + 1;
    end
    
    lineLast = fgetl(fileID);
    meshFormat_Label_End = '$EndMeshFormat';
    assert(strcmp(lineLast,meshFormat_Label_End), ...
        'Unknown file format: line i_line. Expected %s. Found %.', ...
        i_line,meshFormat_Label_End,lineLast)
end

function [i_line,x,nodeIDs] = readNodes_ASCII(lines,i_line,beSilent)
    % number of nodes
    numberOfNodes = str2double(lines{i_line});
    i_line = i_line + 1;
    if not(beSilent)
        fprintf('>>    (Stat.) number of nodes = %i\n',numberOfNodes)
    end
    
    % read in node coordinates (x,y,z)
    if not(beSilent)
        fprintf('>>    (Progress) Reading node coordinates.\n')
    end
    oneLine = sprintf('%s ',lines{i_line:(i_line+numberOfNodes-1)});
    x = sscanf(oneLine,'%f ',4*numberOfNodes);
    x = reshape(x,4,numberOfNodes)';
    nodeIDs = x(:,1);
    x = x(:,2:end);
    
    i_line = i_line + numberOfNodes;
    
    lineEndNodes = lines{i_line};
    label_Nodes_End = '$EndNodes';
    assert(strcmp(lineEndNodes,label_Nodes_End), ...
        'Unknown file format: line %i. Expected %s. Found %s.', ...
        i_line,label_Nodes_End,lineEndNodes)
    i_line = i_line + 1;
end

function [i_line,cElements,cElementTags,elementTypes] = readElements_ASCII(lines,i_line,beSilent,mshElement)
    
    % number of elements
    numberOfElements = str2double(lines{i_line});
    if not(beSilent)
        fprintf('>>    (Stat.) number of elements = %i\n',numberOfElements)
    end
    
    % read in elements
    % (element number, element type, number of tags, tag_1, ..., tag_{number of tags}, 
    %   node number 1, ... , node number n)
    if not(beSilent)
        fprintf('>>    (Progress) Reading elements.\n');
    end
    
    i_line = i_line + 1;
    
    cElements = cell(numberOfElements,1);
    cElementTags = cell(numberOfElements,1);
    elementTypes = zeros(numberOfElements,1);
    
    % Convert elements in chunks.
    % On 64-bit systems / Matlab it should not be necessary.
    % But on 32-bit systems / Matlab, the string may grow beyond the 32-bit
    % limit.
    % The chunking should not affect the speed but is only a little harder
    % to read/maintain.
    it = ceil(numberOfElements/1000000);
    for k = 1:it
        m1 = (k-1)*1000000+1;
        m2 = min(k*1000000,numberOfElements);
        oneLine = sprintf('%s ',lines{(i_line+m1-1):(i_line+m2-1)});
        data = sscanf(oneLine,'%d ');
        
        cnt = 1;
        for i = m1:m2
            cnt = cnt + 1;  % skip element ID
            
            elementType = data(cnt);
            elementTypes(i) = elementType;
            cnt = cnt + 1;

            numberOfTags = data(cnt);
            cnt = cnt + 1;

            cElementTags{i} = data(cnt:(cnt+numberOfTags-1))';
            cnt = cnt + numberOfTags;

            numberOfElementNodes = mshElement(elementType,2);
            cElements{i} = data(cnt:(cnt+numberOfElementNodes-1))';
            cnt = cnt + numberOfElementNodes;
        end
    end
    if not(all(mshElement(elementTypes,1) == elementTypes))
        warning('Cannot read msh file. Unknown element type(s).')
    end
    
    i_line = i_line + numberOfElements;

    lineEndElements = lines{i_line};
    label_Elements_End = '$EndElements';
    assert(strcmp(lineEndElements,label_Elements_End), ...
        'Unknown file format: line %i. Expected %s. Found %s.', ...
        i_line,label_Elements_End,lineEndElements)
    i_line = i_line + 1;
end