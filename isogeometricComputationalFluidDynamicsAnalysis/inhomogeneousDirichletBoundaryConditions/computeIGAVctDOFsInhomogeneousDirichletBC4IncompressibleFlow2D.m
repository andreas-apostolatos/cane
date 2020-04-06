%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function irb = computeIGAVctDOFsInhomogeneousDirichletBC4IncompressibleFlow2D(irb,xi,eta,direction,CP)
%% Function main body
%
% Returns the global numbering of the DOFs on which inhomogeneous Dirichlet
% boundary conditions are prescribed.
%
%     Input :
%       irb : The existing vector with the global numbering of the DOFs where
%             inhomogeneous Dirichlet boundary conditions are applied
%    xi,eta : The surface parameters on the NURBS patch
% direction : The direction of the DOF in the global coordinate system
%        CP : The set of Control Point coordinates and weights
%
%    Output :
%       irb : The updated vector with the global numbering of the DOFs 
%             where inhomogeneous Dirichlet boundary conditions are applied
%
% Function layout :
%
% 0. Read input
%
% 1. Iterate and add new supports preserving the old ones
%
% 2. Sort the output array and delete double entries
%
%% Function main body

%% 0. Read input

% counters
prbCounter = length(irb) + 1;
    
% number of control points in u,v-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

%% 1. Iterate and add new supports preserving the old ones
for j = eta(1)*(neta-1)+1:eta(2)*(neta-1)+1
    for i = xi(1)*(nxi-1)+1:xi(2)*(nxi-1)+1
        % Get the DoF to be prescribed
        irb(prbCounter) = 3*((j-1)*nxi + i-1) + direction;
        
        % Round to nearest integer
        irb(prbCounter) = round(irb(prbCounter));
        
        % Update counter
        prbCounter = prbCounter + 1;
    end
end

%% 2. Sort the output array and delete double entries

% Sort the values
irb = sort(irb);

% Delete double entries
irb = unique(irb);

end

