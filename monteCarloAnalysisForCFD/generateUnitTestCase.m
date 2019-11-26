function testCase = generateUnitTestCase()
    %% Parameters
    Lx          = 2;
    Ly          = 1;
    R           = 1/6;

    nx          = 80;
    ny          = ceil(Ly/Lx * nx);

    %% Allocate
    nodes       = zeros(nx*ny,2);
    elements    = zeros((ny-1)*2*(nx-1),3);

    %% Generate node coordinates
    nodes(:,1)  = repmat(linspace(0,Lx,nx),1,ny)';
    nodes(:,2)  = repelem(linspace(0,Ly,ny),nx)';

    %% Generate elements
    for ky=1:ny-1
        for kx=1:2*(nx-1)

            elemIndex   = (ky-1)*2*(nx-1)+kx;
            BLIndex     = (ky-1)*nx+ceil(kx/2);

            if mod(kx+1,2)
                elements(elemIndex,:)   = [BLIndex+1,BLIndex+nx,BLIndex+nx+1];
            else
                elements(elemIndex,:)   = [BLIndex,BLIndex+1,BLIndex+nx];
            end

        end
    end
    
    %% Delete elements and nodes inside the circle
    % Find nodes that are outside the circle
    outside     = ones(length(nodes),1,'logical');
    for k=1:length(nodes)
        outside(k) = ~isInside(nodes(k,:));
    end
    % Delete nodes inside the circle
    nodes       = nodes(outside,:);
    % Create index map
    indexMap    = zeros(size(outside));
    indexMap(outside) = 1:length(indexMap(outside));
    % Delete elements inside the circle
    outsideElem = ones(length(elements),1,'logical');
    for k=1:length(outsideElem)
        outsideElem(k) = all( outside(elements(k,:)) );
    end
    elements = elements(outsideElem,:);
    % Remap elements
    for k1=1:size(elements,1)
        for k2=1:size(elements,2)
            elements(k1,k2) = indexMap(elements(k1,k2));
        end
    end
    
    %% Homogenous DBC on body
    bNodes          = zeros(length(nodes),1,'logical');
    dx              = 1.5*Lx/(nx-1);
    dy              = 1.5*Ly/(ny-1);
    for k=2:length(nodes)-1
        % Check x
        if nodes(k+1,1)-nodes(k,1) > dx
            bNodes(k)=true;
        elseif nodes(k,1)-nodes(k-1,1) > dx
            bNodes(k)=true;
        end
        % Check y
        column  = abs(nodes(:,1)-nodes(k,1))<1e-14;
        top     = nodes(:,2)>nodes(k,2);
        bottom  = nodes(:,2)<nodes(k,2);
        if min( nodes(column & top,2)-nodes(k,2) ) > dy
            bNodes(k) = true;
        elseif min( nodes(k,2) - nodes(column & bottom,2) ) > dy
            bNodes(k) = true;
        end
    end
    
    temp    = 1:length(nodes);
    bNodes  = temp(bNodes);
    
    %% Homogenous DBC on boundaries
    bottom  = indexMap( 1:nx );
    top     = indexMap( (ny-1)*nx + (1:nx) );
    
    %% Inlet
    inlet   = indexMap( nx+1:nx:(ny-1)*nx );
    
    %% Create structure
    testCase.fldMsh.nodes       = nodes;
    testCase.fldMsh.elements    = elements;
    testCase.homDBC             = 3*[bNodes,bottom',top']';
    testCase.inhomDBC           = 3*inlet;
    testCase.valuesInhomDBC     = ones(size(inlet));

    %% FUNCTION DEFINITIONS
    function out = isInside(point)
        center  = 0.5*[Lx,Ly];
        out = false;
        if center(1)-R<point(1) && point(1)<center(1)+R
            if center(2)-R<point(2) && point(2)<center(2)+R
                out = true;
            end
        end
    end

end