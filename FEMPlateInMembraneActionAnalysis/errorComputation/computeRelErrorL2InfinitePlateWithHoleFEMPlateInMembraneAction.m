function errorL2 = computeRelErrorL2InfinitePlateWithHoleFEMPlateInMembraneAction...
    (strMsh,dHat,parameters,radiusHole,forceAmplitude,propError,int,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the relative error in the L2-norm for the plate in membrane
% action problem of an infinite plate with a hole subject to constant
% pressure at -infty.
%
%               Input :
%              Output :
%              strMsh : On the structural mesh   
%                           .nodes : The nodes in the FE mesh
%                        .elements : The elements in the FE mesh
%                dHat : The solution vector
%          parameters : Technical parameters for the structure
%                               .E : Young's modulus
%                             .nue : Poisson's ratio
%          radiusHole : The radius of the hole
%      forceAmplitude : The amplitude of the applied force at -infty
%           propError : Error computation properties :
%                          .resultant : The resultant on which to compute
%                                       the relative error
%                          .component : The component of the resultant on 
%                                       which to compute the relative error
%                 int : On the integration
%                           .type : 'default' or 'user'
%                           .noGP : Number of Gauss points in case 'user'
%                                   defined quadrature is selected
%              outMsg : Enables outputting information on the command when 
%                       chosen as "outputEnabled"
%
%              Output :
%             errorL2 : The relative error in the L2-norm
%
% Function layout :
%
% 0. Read input
%
% 1. Get integration rule for the computation of the relative error
%
% 2. Get the coordinates of all nodes at once
%
% 3. Create the element freedom tables for all elements at once
%
% 4. Loop over all the quadrature points
% ->
%    4i. Transform the Gauss Point location from the parameter to the physical space
%
%   4ii. Compute the basis functions and their derivatives at the Gauss point page-wise
%
%  4iii. Get the displacement vector at each element page-wise
%
%   4iv. Compute the B-operator matrix for the Voigt strain vector for each element page-wise
%
%    4v. Compute the Voigt strain vector for each element page-wise
%
%   4vi. Compute the Voigt stress vector for each element page-wise
%
%  4vii. Compute the displacement field at the Gauss point for each element pagewise
%
% 4viii. Compute the curvilinear coordinates of the Gauss point for each element pagewise
%
%   4ix. Compute the analytical displacement field
%
%    4x. Compute the analytical stress resultants in the curvilinear basis for each element pagewise
%
%   4xi. Transform the analytical stress resultants into the Cartesian coordinate system for each element pagewise
%
%  4xii. Compute the difference of the numerical and the exact fields on the Gauss point for each element pagewise
%
% 4xiii. Compute the L2 norm of the exact field at the Gauss point
% <-
%
% 5. Compute the relative error in the L2-norm
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________________\n');
    fprintf('##################################################################\n');
    fprintf('Computation of the relative error in the L2-norm for the benchmark\n');
    fprintf('problem of a curved beam-like plate in membrane action discretized\n');
    fprintf('with Finite Elements subject to tip shear pressure load has been \n');
    fprintf('initiated \n\n');
    fprintf('Relative error in the L2-norm for the ');
    if strcmp(propError.resultant,'displacement')
        if strcmp(propError.component,'x')
            fprintf('displacement component d_x\n');
        elseif strcmp(propError.component,'y')
            fprintf('displacement component d_y\n');
        elseif strcmp(propError.component,'2norm')
            fprintf('displacement vector d\n');
        else
            error('Select component of the displacement')
        end
    elseif strcmp(propError.resultant,'strain')
        error('The computation of the relative error in the L2-norm for the strain has not yet been implemented');
        if strcmp(propError.component,'x')
            fprintf('strain component epsilon_xx\n');
        elseif strcmp(propError.component,'y')
            fprintf('strain component epsilon_yy\n');
        elseif strcmp(propError.component,'xy')
            fprintf('strain component epsilon_xy\n');
        elseif strcmp(propError.component,'2norm')
            fprintf('2-norm of the strain tensor\n');
        end
    elseif strcmp(propError.resultant,'stress')
        if strcmp(propError.component,'x')
            fprintf('stress component sigma_xx\n');
        elseif strcmp(propError.component,'y')
            fprintf('stress component sigma_yy\n');
        elseif strcmp(propError.component,'xy')
            fprintf('stress component sigma_xy\n');
        elseif strcmp(propError.component,'2norm')
            fprintf('2-norm of the stress tensor\n');
        end
    end
    fprintf('Number of Gauss Points : %d\n',int.noGP);
    fprintf('__________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Constant values for the computation of the analytical resultants
mue = parameters.E/(2*(1 + parameters.nue));
kappa = (3 - parameters.nue)/(1 + parameters.nue);

% Number of nodes at the element level
noNodesEl = 3;

% Number of DOFs per node
noDOFsPerNode = 2;

% Number of degrees of freedom per element
noDOFsEl = noDOFsPerNode*noNodesEl;

% Total number of elements in the mesh
noElmnts = length(strMsh.elements(:,1));

% Initialize variables
errorInL2 = 0;
exactInL2 = 0;

% Compute the material matrix pagewise
matMtx = zeros(noElmnts,3,3);
for counterElmnts = 1:noElmnts
    matMtx(counterElmnts,:,:) = parameters.E/(1-parameters.nue^2)*...
        [1              parameters.nue 0
         parameters.nue 1              0 
         0              0              (1-parameters.nue)/2];
end

%% 1. Get integration rule for the computation of the relative error
if strcmp(int.type,'default')
    noGP = 1;
elseif strcmp(int.type,'user')
    noGP = int.noGP;
end
[GP,GW] = getGaussRuleOnCanonicalTriangle(noGP);

%% 2. Get the coordinates of all nodes at once
nodes1 = strMsh.nodes(strMsh.elements(:,1),:);
nodes2 = strMsh.nodes(strMsh.elements(:,2),:);
nodes3 = strMsh.nodes(strMsh.elements(:,3),:);

%% 3. Create the element freedom tables for all elements at once
EFT = zeros(noDOFsEl,noElmnts);
for counterEFT = 1:noNodesEl
    for counterDOFsPerNode = 1:noDOFsPerNode-1
        EFT(noDOFsPerNode*counterEFT, :) = noDOFsPerNode*strMsh.elements(:,counterEFT)';
        EFT(noDOFsPerNode*counterEFT-(noDOFsPerNode-counterDOFsPerNode), :) = ...
            EFT(noDOFsPerNode*counterEFT, :) - (noDOFsPerNode-counterDOFsPerNode);
    end
end

%% 4. Loop over all the quadrature points
for iGP = 1:noGP
    %% 4i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(iGP,1)*nodes1 + GP(iGP,2)*nodes2 + (1 - GP(iGP,1) - GP(iGP,2))*nodes3;
    
    %% 4ii. Compute the basis functions and their derivatives at the Gauss point page-wise
    [dN,area] = computeCST2DBasisFunctionsAndFirstDerivatives...
            (nodes1,nodes2,nodes3,xGP(:,1,:),xGP(:,2,:));
    detJxxi = 2.*area;
        
	%% Form the basis functions matrix at the Gauss Point page-wise
    if strcmp(propError.resultant,'displacement')
        NMtx = zeros(noElmnts,2,noDOFsEl);
        for i = 1:noNodesEl
            NMtx(:,1,2*i-1) = dN(:,i,1);
            NMtx(:,2,2*i) = dN(:,i,1);
        end
    end

    %% 4iii. Get the displacement vector at each element page-wise
    dHatEl = dHat(EFT');
    
    %% 4iv. Compute the B-operator matrix for the Voigt strain vector for each element page-wise
    if strcmp(propError.resultant,'strain') || strcmp(propError.resultant,'stress')
        BMtx = zeros(noElmnts,3,noDOFsEl);
        for i = 1:noNodesEl
            BMtx(:,1,2*i-1) = dN(:,i,2);
            BMtx(:,2,2*i) = dN(:,i,3);
            BMtx(:,3,2*i-1) = dN(:,i,3);
            BMtx(:,3,2*i) = dN(:,i,2);
        end
    end
    
    %% 4v. Compute the Voigt strain vector for each element page-wise
    if strcmp(propError.resultant,'strain') || strcmp(propError.resultant,'stress')
        strainVoigt = pmtimes(BMtx,dHatEl);
        if strcmp(propError.component,'x')
            resultant = strainVoigt(:,1);
        elseif strcmp(propError.component,'y')
            resultant = strainVoigt(:,2);
        elseif strcmp(propError.component,'xy')
            resultant = strainVoigt(:,3);
        elseif strcmp(propError.component,'2norm')
            resultant  = [strainVoigt(:,1) strainVoigt(:,2) 2*strainVoigt(:,3)];
        end
    end
    
    %% 4vi. Compute the Voigt stress vector for each element page-wise
    if strcmp(propError.resultant,'stress')
        stressVoigt = pmtimes(matMtx,strainVoigt);
        if strcmp(propError.component,'x')
            resultant = stressVoigt(:,1);
        elseif strcmp(propError.component,'y')
            resultant = stressVoigt(:,2);
        elseif strcmp(propError.component,'xy')
            resultant = stressVoigt(:,3);
        elseif strcmp(propError.component,'2norm')
            resultant  = [stressVoigt(:,1) stressVoigt(:,2) 2*stressVoigt(:,3)];
        end
    end
    
    %% 4vii. Compute the displacement field at the Gauss point for each element pagewise
    if strcmp(propError.resultant,'displacement')
        disp = pmtimes(NMtx,dHatEl);
        if strcmp(propError.component,'x')
            resultant = disp(:,1);
        elseif strcmp(propError.component,'y')
            resultant = disp(:,2);
        elseif strcmp(propError.component,'2norm')
            resultant = disp;
        end
    end
    
    %% 4viii. Compute the curvilinear coordinates of the Gauss point for each element pagewise
    theta1 = sqrt(xGP(:,1).^2 + xGP(:,2).^2);
    theta2 = 2*pi() - atan(xGP(:,2)./(-xGP(:,1)));
    
    %% 4ix. Compute the analytical displacement field
    if strcmp(propError.resultant,'displacement')
        dx = -forceAmplitude*radiusHole*(theta1*(kappa + 1).*cos(theta2)/radiusHole + 2*radiusHole*((1 + kappa)*cos(theta2) + cos(3*theta2))./theta1 - 2*radiusHole^3*cos(3*theta2)./theta1.^3)/mue/8;
        dy = -forceAmplitude*radiusHole*(theta1*(kappa - 3).*sin(theta2)/radiusHole + 2*radiusHole*((1 - kappa)*sin(theta2) + sin(3*theta2))./theta1 - 2*radiusHole^3*sin(3*theta2)./theta1.^3)/mue/8;
        if strcmp(propError.component,'x')
            resultantEx = dx;
        elseif strcmp(propError.component,'y')
            resultantEx = dy;
        elseif strcmp(propError.component,'2norm')
            resultantEx = [dx dy];
        end
    end
    
    %% 4x. Compute the analytical stress resultants in the curvilinear basis for each element pagewise
    if strcmp(propError.resultant,'stress')
        pr = forceAmplitude/2*(1 - (radiusHole./theta1).^2) + forceAmplitude/2*(1 - 4*(radiusHole./theta1).^2 + 3*(radiusHole./theta1).^4).*cos(2*theta2);
        pt = forceAmplitude/2*(1 + (radiusHole./theta1).^2) - forceAmplitude/2*(1 + 3*(radiusHole./theta1).^4).*cos(2*theta2);
        prt = -forceAmplitude/2*(1 + 2*(radiusHole./theta1).^2 - 3*(radiusHole./theta1).^4).*sin(2*theta2);
    end
    
    %% 4xi. Transform the analytical stress resultants into the Cartesian coordinate system for each element pagewise
    if strcmp(propError.resultant,'stress')
        px = pr.*cos(theta2).^2 + pt.*sin(theta2).^2 - 2*prt.*sin(theta2).*cos(theta2);
        py = pr.*sin(theta2).^2 + pt.*cos(theta2).^2 + 2*prt.*sin(theta2).*cos(theta2);
        pxy = pr.*sin(theta2).*cos(theta2) - pt.*sin(theta2).*cos(theta2) - prt.*(sin(theta2).^2 - cos(theta2).^2);
        if strcmp(propError.component,'x')
            resultantEx = px;
        elseif strcmp(propError.component,'y')
            resultantEx = py;
        elseif strcmp(propError.component,'xy')
            resultantEx = pxy;
        elseif strcmp(propError.component,'2norm')
            resultantEx = [px py 2*pxy];
        end
    end
    
    %% 4xii. Compute the difference of the numerical and the exact fields on the Gauss point for each element pagewise
    errorInL2 = errorInL2 + sum((resultant - resultantEx).^2,2).*detJxxi*GW(iGP,1);
    
    %% 4xiii. Compute the L2 norm of the exact field at the Gauss point
    exactInL2 = exactInL2 + sum(resultantEx.^2,2).*detJxxi*GW(iGP,1);
end

%% 5. Compute the relative error in the L2-norm
errorInL2 = sqrt(sum(errorInL2));
exactInL2 = sqrt(sum(exactInL2));
errorL2 = errorInL2/exactInL2;

end
