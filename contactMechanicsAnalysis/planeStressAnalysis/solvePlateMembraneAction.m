%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. K.-U. Bletzinger                 %
%   _____________________________________________________                 %
%                                                                         %
%   Authors                                                               %
%   _______                                                               %
%                                                                         %
%   Dr.-Ing. Roland Wüchner                                               %
%   Dipl.-Math. Andreas Apostolatos (andreas.apostolatos@tum.de)          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displacement = solvePlateMembraneAction(mesh,rb,F,materialProperties,analysis)
%% Function documentation
%
% Returns the displacement field corresponding to a plain stress/strain
% analysis for the given mesh of the geometry together with its Dirichlet
% and Neumann boundary conditions.
%
%              Input :
%               mesh : Elements and nodes of the mesh
%                 rb : Vector of the prescribed DoFs (by global numbering)
%                  F : Global load vector
% materialProperties : The material properties of the structure
%           analysis : Analysis type (plane stress or plane strain)
%      
%             Output :
%       displacement : The resulting displacement field
%
% Function layout :
%
% 1. Compute the master stiffness matrix of the structure
%
% 2. Reduce the system according to the given constraints
%
% 3. Compute reduced displacement vector
%
% 4. Assemble to the complete displacement vector
%
% 5. compute the complete load vector and verify the results
%
%% Function main body
if strcmp(analysis.physics,'plain_stress')
    fprintf('Plain stress analysis has been initiated \n');
elseif strcmp(analysis.physics,'plain_strain')
    fprintf('Plain strain analysis has been initiated \n');
end
fprintf('\n');

%% 1. Compute the master stiffness matrix of the structure

fprintf('\t Computing master stiffness matrix... \n');

% Master stiffness matrix
K = computeStiffnessMatrixPlateMembraneActionLinear(mesh,materialProperties,analysis);

%% 2. Reduce the system according to the given constraints

fprintf('\t Reducing the system according to the constraints... \n');

Kred = K;
Fred = F;

for i = length(rb):-1:1
    Kred(:,rb(i)) = [];
    Kred(rb(i),:) = [];
    Fred(rb(i)) = [];
end

%% 3. Compute reduced displacement vector

fprintf('\t Solving the linear system of %d equations... \n',length(F));
dred = Kred\Fred;

%% 4. Assemble to the complete displacement vector

fprintf('\t Re-arranging the global displacement vector... \n');

% Dimension of the complete displacement vector
displacement = zeros(length(F),1);

% Initialize counters
i=1;
k=1;
for l = 1:length(F)
    if (i<=length(rb))
        % if we are in a Dirichlet boundary condition location add 0
        if (l==rb(i))
            displacement(l,1) = 0;
            % update counter
            i=i+1;
        % if not add the Control Point displacement 
        else
            displacement(l,1) = dred(k);
            % update counter
            k=k+1;
        end
    else
        displacement(l,1) = dred(k);
        % update counter
        k=k+1;
    end
end

%% 5. compute the complete load vector and verify the results

fprintf('\t Verifying the results... \n');
fprintf('\n');

F_verification = K*displacement;

%  A tolerance value
tolerance = 10;

% Compute a residual
residual = norm(F_verification-F);

if residual>=tolerance
   fprintf('\t The computed load vector has been detected to be different \n');
   fprintf('\t than the global load vector in the range of %d \n',tolerance);
end

end