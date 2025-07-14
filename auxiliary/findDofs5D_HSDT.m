function homDOFs = findDofs5D_HSDT(homDOFs, xiSup, etaSup, dirSupp, CP)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    [Your Name] (based on Andreas Apostolatos)
%
%% Function documentation
%
% Returns the global numbering of the DOFs where homogeneous Dirichlet 
% boundary conditions are applied for HSDT formulation with 5 DOF per 
% control point: [u, v, w, θx, θy]
%
%                Input :
%              homDOFs : The existing array with the global numbering of 
%                        the DOFs where homogeneous Dirichlet boundary 
%                        conditions are applied
%                xiSup : Extension of the support in xi-direction
%               etaSup : Extension of the support in eta-direction  
%              dirSupp : The direction of the support:
%                        1 : u-displacement
%                        2 : v-displacement  
%                        3 : w-displacement
%                        4 : θx-rotation
%                        5 : θy-rotation
%                   CP : The Control Point coordinates and weights
%
%               Output :
%              homDOFs : The updated array with the global numbering of the
%                        DOFs where homogeneous Dirichlet boundary 
%                        conditions are applied
%
% Function layout :
%
% 1. Read input
%
% 2. Find the Control Points to be supported
%
% 3. Find the corresponding global numbering (5 DOF per control point)
%
%% Function main body

%% 1. Read input
nxi = length(CP(:, 1, 1));
neta = length(CP(1, :, 1));

%% 2. Find the Control Points to be supported

% Initialize counters
counter = 1;

% Find the parametric locations in xi-direction
if xiSup(1) == xiSup(2)
    % Support along a xi = constant line
    if xiSup(1) == 0
        xiLoc = 1;
    elseif xiSup(1) == 1
        xiLoc = nxi;
    else
        error('Support location xi = %f not supported. Use xi = 0 or xi = 1', xiSup(1));
    end
    xiIndices = xiLoc;
else
    % Support along an xi interval
    error('Interval supports in xi-direction not yet implemented for HSDT');
end

% Find the parametric locations in eta-direction  
if etaSup(1) == etaSup(2)
    % Support along an eta = constant line
    if etaSup(1) == 0
        etaLoc = 1;
    elseif etaSup(1) == 1
        etaLoc = neta;
    else
        error('Support location eta = %f not supported. Use eta = 0 or eta = 1', etaSup(1));
    end
    etaIndices = etaLoc;
else
    % Support along an eta interval
    if etaSup(1) == 0 && etaSup(2) == 1
        etaIndices = 1:neta;
    else
        error('Only full eta interval [0,1] supported for HSDT');
    end
end

%% 3. Find the corresponding global numbering (5 DOF per control point)

% Loop over the identified control points
for etaIndex = etaIndices
    for xiIndex = xiIndices
        % Compute the global control point number
        globalCPIndex = (etaIndex - 1) * nxi + xiIndex;
        
        % Compute the global DOF number for the specified direction
        % HSDT DOF ordering: [u, v, w, θx, θy] per control point
        globalDOFIndex = 5 * (globalCPIndex - 1) + dirSupp;
        
        % Add to the homogeneous Dirichlet boundary conditions array
        homDOFs(counter) = globalDOFIndex;
        counter = counter + 1;
    end
end

% Remove duplicates and sort
homDOFs = unique(homDOFs);

end