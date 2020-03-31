function errorL2 = computeRelErrorL2CurvedBeamTipShearFEMPlateInMembraneAction...
    (strMsh, dHat, parameters, internalRadius, externalRadius, forceAmplitude, ...
    propError, propInt, outMsg)
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
% action problem of a curved beam-like plate subject to end tip shear force
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
%      internalRadius : The internal radius of the beam-like curved plate
%      forceAmplitude : The amplitude of the applied force at the tip
%           propError : Error computation properties :
%                          .resultant : The resultant on which to compute
%                                       the relative error
%                          .component : The component of the resultant on 
%                                       which to compute the relative error
%             propInt : Structure containing properties on the numerical
%                       integration,
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
%  4vii. Compute the curvilinear coordinates of the Gauss point for each element pagewise
%
% 4viii. Transform the analytical stress resultants into the Cartesian coordinate system for each element pagewise
%
%   4ix. Compute the difference of the numerical and the exact fields on the Gauss point for each element pagewise
%
%    4x. Compute the L2 norm of the exact field at the Gauss point
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
        error('The computation of the relative error in the L2-norm for the displacement has not yet been implemented');
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
        end
    elseif strcmp(propError.resultant,'stress')
        if strcmp(propError.component,'x')
            error('The computation of the relative error in the L2-norm for the tensor component sigma_xx has not yet been implemented');
            fprintf('stress component sigma_xx\n');
        elseif strcmp(propError.component,'y')
            error('The computation of the relative error in the L2-norm for the tensor component sigma_yy has not yet been implemented');
            fprintf('stress component sigma_yy\n');
        elseif strcmp(propError.component,'xy')
            error('The computation of the relative error in the L2-norm for the tensor component sigma_xy has not yet been implemented');
            fprintf('stress component sigma_xy\n');
        elseif strcmp(propError.component,'tensor')
            fprintf('stress tensor sigma\n');
        end
    end
    fprintf('\nNumber of Gauss Points : %d\n',propInt.noGP);
    fprintf('__________________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Constant value for the computation of the analytical resultants
N = internalRadius^2 - externalRadius^2 + (internalRadius^2 + externalRadius^2)*...
    log(externalRadius/internalRadius);

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
if strcmp(propInt.type,'default')
    noGP = 1;
elseif strcmp(propInt.type,'user')
    noGP = propInt.noGP;
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
for counterGP = 1:noGP
    %% 4i. Transform the Gauss Point location from the parameter to the physical space
    xGP = GP(counterGP,1)*nodes1 + GP(counterGP,2)*nodes2 + (1-GP(counterGP,1)-GP(counterGP,2))*nodes3;
    
    %% 4ii. Compute the basis functions and their derivatives at the Gauss point page-wise
    [dN,~] = computeCST2DBasisFunctionsAndFirstDerivatives...
            (nodes1,nodes2,nodes3,xGP(:,1,:),xGP(:,2,:));
        
% 	%% Form the basis functions matrix at the Gauss Point page-wise
%     NMtx = zeros(noElmnts,2,noDOFsEl);
%     for i = 1:noNodesEl
%         NMtx(:,1,noDOFsPerNode*i-noDOFsPerNode+1) = dN(:,i,1);
%         NMtx(:,2,noDOFsPerNode*i-noDOFsPerNode+2) = dN(:,i,1);
%     end

    %% 4iii. Get the displacement vector at each element page-wise
    dHatEl = dHat(EFT');
    
    %% 4iv. Compute the B-operator matrix for the Voigt strain vector for each element page-wise
    BMtx = zeros(noElmnts,3,noDOFsEl);
    for i = 1:noNodesEl
        BMtx(:,1,2*i-1) = dN(:,i,2);
        BMtx(:,2,2*i) = dN(:,i,3);
        BMtx(:,3,2*i-1) = dN(:,i,3);
        BMtx(:,3,2*i) = dN(:,i,2);
    end
    
    %% 4v. Compute the Voigt strain vector for each element page-wise
    strainVoigt = pmtimes(BMtx,dHatEl);
    
    %% 4vi. Compute the Voigt stress vector for each element page-wise
    stressVoigt = pmtimes(matMtx,strainVoigt);
    
%     %% Compute the displacement field at the Gauss point for each element pagewise
%     dElOnGP = pmtimes(NMtx,dHatEl);
    
    %% 4vii. Compute the curvilinear coordinates of the Gauss point for each element pagewise
    r = sqrt(xGP(:,1).^2 + xGP(:,2).^2);
    theta = atan(xGP(:,2)./xGP(:,1));
    
    %% 4vii. Compute the analytical stress resultants in the curvilinear basis for each element pagewise
    pr = forceAmplitude*(r + internalRadius^2*externalRadius^2./r.^3 - (internalRadius^2+externalRadius^2)./r).*sin(theta)/N;
    pt = forceAmplitude*(3*r - internalRadius^2*externalRadius^2./r.^3 - (internalRadius^2 + externalRadius^2)./r).*sin(theta)/N;
    prt = - forceAmplitude*(r + internalRadius^2*externalRadius^2./r.^3 - (internalRadius^2 + externalRadius^2)./r).*cos(theta)/N;
    
    %% 4viii. Transform the analytical stress resultants into the Cartesian coordinate system for each element pagewise
    px = pr.*cos(theta).^2 + pt.*sin(theta).^2 - 2*prt.*sin(theta).*cos(theta);
    py = pr.*sin(theta).^2 + pt.*cos(theta).^2 + 2*prt.*sin(theta).*cos(theta);
    pxy = pr.*sin(theta).*cos(theta) - pt.*sin(theta).*cos(theta) - prt.*(sin(theta).^2 - cos(theta).^2);
    stressVoigtEx = zeros(noElmnts,3);
    stressVoigtEx(:,1) = px;
    stressVoigtEx(:,2) = py;
    stressVoigtEx(:,3) = pxy;
    
    %% 4ix. Compute the difference of the numerical and the exact fields on the Gauss point for each element pagewise
    errorInL2 = errorInL2 + ...
        ((stressVoigt(:,1) - stressVoigtEx(:,1)).^2 + ...
        (stressVoigt(:,2) - stressVoigtEx(:,2)).^2 + ...
        (stressVoigt(:,3) - stressVoigtEx(:,3)).^2)*GW(counterGP,1);
    
    %% 4x. Compute the L2 norm of the exact field at the Gauss point
    exactInL2 = exactInL2 + ...
        (stressVoigtEx(:,1).^2 + ...
        stressVoigtEx(:,2).^2 + ...
        stressVoigtEx(:,3).^2)*GW(counterGP,1);
end

%% 5. Compute the relative error in the L2-norm
errorInL2 = sqrt(sum(errorInL2));
exactInL2 = sqrt(sum(exactInL2));
errorL2 = errorInL2/exactInL2;

end
