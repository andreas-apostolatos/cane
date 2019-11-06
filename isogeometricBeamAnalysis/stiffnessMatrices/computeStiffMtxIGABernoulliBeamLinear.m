function [KGlobal,F,minElSize] = computeStiffMtxIGABernoulliBeamLinear...
    (KConstant,dHat,dHatSaved,dHatDot,dHatDotSaved,BSplinePatch,connections,...
    propCoupling,propTransientAnalysis,t,tab,loadFactor,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the master stiffness matrix corresponding to the isogeometric
% Bernoulli beam formulation.
%
%                 Input :
%             KConstant : Dummy variable for this function
%                  dHat : The displacement field from the previous iteration 
%                         step (dummy variable for this function)
%             dHatSaved : The displacement field from the previous time 
%                         step (dummy variable for this function)
%               dHatDot : The velocity field from the previous iteration 
%                         step (dummy variable for this function)
%          dHatDotSaved : The velocity field from the previous time step 
%                         (dummy variable for this function)
%          BSplinePatch : Structure containing all geometrical, technical 
%                         and load information for the B-Spline patch
%           connections : Structure on the way that several multipatches
%                         are connected to each other (dummy variable for 
%                         this function)
%          propCoupling : Dummy variable for this function
% propTransientAnalysis : Dummy variable for this function
%                     t : Time instance
%                   tab : Tabulation for the outputting on the command
%                         window
%            loadFactor : Dummy variable for this function
%                outMsg : Dummy variable for this function 
%
%
%                Output :
%               KGlobal : The global master stiffness matrix corresponding 
%                         to the Timoshenko beam formulation
%                     F : The right-hand side load vector
%             minElSize : The minimum element size in the isogeometric mesh
% 
% Function layout :
%
% 0. Read input
%
% 1. Choose the integration scheme
%
% 2. Loop over all the elements
% ->
%    2i. Initialize the element area size
%
%   2ii. Compute the determinant of the Jacobian to the transformation from the parameter to the integration space
%
%  2iii. Create an element freedom table
%
%   2iv. Loop over all the Gauss Points
%   ->
%        2iv.1. Linear transformation from the quadrature domain to the knot span
%
%        2iv.2. Find the correct knot span where xi lives
%
%        2iv.3. Compute the NURBS basis functions and their first and second derivatives on GP
%
%        2iv.4. Compute the base vectors and their derivatives on GP
%
%        2iv.5. Compute the determinant of the Jacobian to the transformation from the NURBS parameter to the physical space
%
%        2iv.6. Compute the element size on the Gauss Point and sum up the contribution to it
%
%        2iv.7. Compute the Timoshenko element stiffness matrix on GP
%
%        2iv.8. Add the contribution from the Gauss point and assemble to the global stiffness matrix
%   <-
%    2v. Check if the current minimum element area size is smaller than for the previous element
% <-
%
% 3. Compute the exernally applied load vector
%
%% Function main body

%% 0. Read input

% Check the NURBS geometry input
if iscell(BSplinePatch)
    if length(BSplinePatch) > 1
        error('Multipatch NURBS surface is given as input to the computation of the stiffness matrix for a single patch NURBS surface');
    else
        BSplinePatch = BSplinePatch{1};
    end
end

% Get the B-Spline patch parameters
p = BSplinePatch.p;
Xi = BSplinePatch.Xi;
CP = BSplinePatch.CP;
isNURBS = BSplinePatch.isNURBS;
parameters = BSplinePatch.parameters;
int = BSplinePatch.int;
NBC = BSplinePatch.NBC;

% Number of Control Points
nxi = length(CP(:,1));

% Initialization of the Global Stiffness matrix
noDOFs = 2*nxi;

% Number of Control Points in the element level
noCPsEl = p + 1;

% Number of Degrees of Freedom in the element level
noDOFsEl = 2*noCPsEl;

% Initialize minimum element area size
minElementSize = norm(CP(1,1:3) - CP(nxi,1:3));
minElSize = minElementSize;

% Initialize output array
KGlobal = zeros(noDOFs,noDOFs);

%% 1. Choose the integration scheme

% choose the number of the Gauss points accordingly
if strcmp(int.type,'default')
    noGP = ceil(0.5*(2*p - 1)) ;
elseif strcmp(int.type,'user')
    noGP = int.noGP;
end

% Get Gauss Points and weights
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(noGP);

%% 2. Loop over all the elements
for i = p+1:length(Xi)-p-1
    % Check if we are in a non zero span
    if Xi(i+1)-Xi(i)~=0
        %% 2i. Initialize the element area size
            
        % Initialize the new one
        minElementSize = 0;
        
        %% 2ii. Compute the determinant of the Jacobian to the transformation from the parameter to the integration space
        detJxiu = (Xi(i+1)-Xi(i))/2.0;
        
        %% 2iii. Create an element freedom table
        
        % Initialize element freedom table
        EFT = zeros(1,noDOFsEl);

        % Initialize counter
        k = 1;
        
        % Relation global-local DoFs
        for cpi = i-BSplinePatch.p:i
            EFT(k)   = BSplinePatch.DOFNumbering(cpi,1);
            EFT(k+1) = BSplinePatch.DOFNumbering(cpi,2);

            % update counter
            k = k + 2;
        end
        
        %% 2iv. Loop over all the Gauss Points
        for j = 1:length(GW)
            %% 2iv.1. Linear transformation from the quadrature domain to the knot span
            xi = ((1-GP(j))*Xi(i)+(1+GP(j))*Xi(i+1))/2;
            
            %% 2iv.2. Find the correct knot span where xi lives
            kontSpan = findKnotSpan(xi,Xi,nxi);
            
            %% 2iv.3. Compute the NURBS basis functions and their first and second derivatives on GP
            dR = computeIGABasisFunctionsAndDerivativesForCurve(kontSpan,p,xi,Xi,CP,isNURBS,2);
            
            %% 2iv.4. Compute the base vectors and their derivatives on GP
            [G,dG] = computeBaseVectorNormalToNURBSCurveAndDeivativesForCurve2D(kontSpan,p,CP,dR);
            
            %% 2iv.5. Compute the determinant of the Jacobian to the transformation from the NURBS parameter to the physical space
            detJxxi = norm(G(:,1));
            
            %% 2iv.6. Compute the element size on the Gauss Point and sum up the contribution to it
            minElementSizeOnGP = abs(detJxxi)*abs(detJxiu)*GW(j);
            minElementSize = minElementSize + minElementSizeOnGP;
            
            %% 2iv.7. Compute the Timoshenko element stiffness matrix on GP
            KelementOnGP = computeElStiffMtxIGABernoulliBeamLinear(p,G(:,1),dG(:,1),dR,parameters);

            %% 2iv.8. Add the contribution from the Gauss point and assemble to the global stiffness matrix
            KGlobal(EFT,EFT) = KGlobal(EFT,EFT) + KelementOnGP*minElementSizeOnGP;
        end
    end
    
    %% 2v. Check if the current minimum element area size is smaller than for the previous element
    if minElementSize <= minElSize
        minElSize = minElementSize;
    end
end

%% 3. Compute the exernally applied load vector
F = zeros(noDOFs,1);
for counterNBC = 1:NBC.noCnd
    funcHandle = str2func(NBC.computeLoadVct{counterNBC});
    F = funcHandle(F,NBC.xiLoadExtension{counterNBC},p,Xi,...
        CP,isNURBS,NBC.loadAmplitude(counterNBC,1),...
        NBC.loadDirection(counterNBC,1),t,int,'');
end
 
end
