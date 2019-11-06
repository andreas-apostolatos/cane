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
function [X,Y,Z,sigma0] = createBSplineSurfaceOnCartesianSpace...
    (p,q,Xi,Eta,CP,isNURBS,prestress,compPrestress,xiGrid,etaGrid)
%% Function documentation
%
% Returns three arrays, containing the Cartesian coordinates of the points 
% on the NURBS surface in a grid of xiGrid Times etaGridv lines
%
%            Input :
%              p,q : Polynomial degrees
%           Xi,Eta : Knot vectors in xi,eta-direction
%               CP : Control point coordinates and weights
%          isNURBS : Flag on whether the geometrical basis is NURBS or B-Spline
%        prestress : The Voigt tensor of the prestress [px; py; pxy] or 
%                    'undefined' if not defined
%    compPrestress : '1' '2' or '12' for the component of the prestress to
%                    be computed
%           xiGrid : Number of lines in xi-direction
%          etaGrid : Number of lines in eta-direction
%
%           Output :
%                X : Array containing the x-coordinates of the points on 
%                    the surface
%                Y : Array containing the y-coordinates of the points on 
%                    the surface
%                Z : Array containing the z-coordinates of the points on 
%                    the surface
%           sigma0 : The prestress tensor computed as a tensor product at 
%                    the evaluation points
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the parametric coordinate locations
% ->
%    1i. Compute start coordinate on the xi parameter lines and initialize counter
%
%   1ii. Compute the knot span index in eta-direction
%
%  1iii. Loop over all the coordinates in xi-direction
%  ->
%        1iii.1. Compute the knot span index in xi-direction
%
%        1iii.2. Compute the IGA basis functions at (xi,eta)
%
%        1iii.3. Compute the Cartesian coordinates of (xi,eta)
%
%        1iii.4. Compute the local to the convective space of the surface coordinates
%
%        1iii.5. Compute the base vectors of the patch at (xi,eta) in case the prestress is defined over a user defined coordinate system
%
%        1iii.6. Compute the transformation matrix from the user-defined coordinate system of the prestress to the local Cartesian system
%
%        1iii.7. Compute the prestress components according to the options
%
%        1iii.8. Update counter for the lines in xi-direction
%
%        1iii.9. Update the xi-parametric coordinate
%  <-
%
%   1iv. Update counter for the lines in eta-direction
%
%    1v. Update the eta-parametric coordinate
% <-
%
% 2. Write the coordinates into the individual arrays
%
%% Function main body

%% 0. Read input

% Number of knots in xi,eta-direction
mxi = length(Xi);
meta = length(Eta);

% Number of control points in xi,eta-direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Check the compatibility of the NURBS parameters
checkInputForBSplineSurface(p,mxi,nxi,q,meta,neta);

% Assign a tolerance value
eps = 10e-10;

% Initialize counter for the lines in eta-direction
counterEta = 1;  

% Incremental step for the lines in eta-direction
deta = (Eta(meta)-Eta(1))/etaGrid;    

% Incremental step for the lines in xi-direction
dxi = (Xi(mxi)-Xi(1))/xiGrid;  

% Start coordinate on the eta parameter lines
eta = Eta(1);

% Read input
if ~ischar(prestress)
    if ~strcmp(compPrestress,'1') && ~strcmp(compPrestress,'2') && ~strcmp(compPrestress,'12')
        error('The value of optionsPrestress.componentPrestress can be "1", "2" or "12"');
    end
end

% Check if prestress is defined over a user defined coordinates system
if ~ischar(prestress)
    isPrestressOverDefinedSystem = false;
    if isfield(prestress,'computeBaseVectors')
        if ~isfield(prestress,'computeParametricCoordinates')
            error('When defining the prestress over a user-defined coordinate system the function handle prestress.computeParametricCoordinates needs to also be defined');
        end
        isPrestressOverDefinedSystem = true;
    end
end

% Initialize output arrays
S = zeros(xiGrid,etaGrid,3);
if ~ischar(prestress)
    sigma0 = zeros(xiGrid,etaGrid);
else
    sigma0 = 'undefined';
end

%% 1. Loop over all the parametric coordinate locations
while eta <= Eta(meta)+eps
    %% 1i. Compute start coordinate on the xi parameter lines and initialize counter
    xi = Xi(1);
    
    % Initialize counter for the lines in u-direction
    counterXi = 1;
    
    %% 1ii. Compute the knot span index in eta-direction
    etaSpan = findKnotSpan(eta,Eta,neta);
    
    %% 1iii. Loop over all the coordinates in xi-direction
    while xi <= Xi(mxi) + eps
        %% 1iii.1. Compute the knot span index in xi-direction
    	xiSpan = findKnotSpan(xi,Xi,nxi);
    	
        %% 1iii.2. Compute the IGA basis functions at (xi,eta)
        if ~ischar(prestress)
            if isPrestressOverDefinedSystem
                nDrv = 1;
            else
                nDrv = 0;
            end
        else
            nDrv = 0;
        end
        dR = computeIGABasisFunctionsAndDerivativesForSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,isNURBS,nDrv);
        
        %% 1iii.3. Compute the Cartesian coordinates of (xi,eta)
        S(counterXi,counterEta,1:3) = computeCartesianCoordinatesOfAPointOnBSplineSurface(xiSpan,p,xi,Xi,etaSpan,q,eta,Eta,CP,dR(:,1));

        %% 1iii.4. Compute the local to the convective space of the surface coordinates
        if ~ischar(prestress)
            if isPrestressOverDefinedSystem || isa(prestress.voigtVector,'function_handle')
                theta = real(prestress.computeParametricCoordinates(squeeze(S(counterXi,counterEta,1:3))));
            end
        end
        
        %% 1iii.5. Compute the base vectors of the patch at (xi,eta) in case the prestress is defined over a user defined coordinate system
        if ~ischar(prestress)
            if isPrestressOverDefinedSystem
                mixedDerivOrder = 0;
                [dGXi,dGEta] = computeBaseVectorsAndDerivativesForBSplineSurface...
                    (xiSpan,p,etaSpan,q,CP,mixedDerivOrder,dR);
            end
        end
        
        %% 1iii.6. Compute the transformation matrix from the user-defined coordinate system of the prestress to the local Cartesian system
        if ~ischar(prestress)
            if isPrestressOverDefinedSystem
                GabCovariant = [dGXi(:,1) dGEta(:,1)]'*[dGXi(:,1) dGEta(:,1)];
                GContravariant = GabCovariant\[dGXi(:,1) dGEta(:,1)]';
                GContravariant = GContravariant';
                eLC = computeLocalCartesianBasis4BSplineSurface([dGXi(:,1) dGEta(:,1)],GContravariant);
                prestressBaseVct = prestress.computeBaseVectors(theta(1,1),theta(2,1));
                T2LC = computeT2LocalCartesianBasis(prestressBaseVct,eLC);
            end
        end
        
        %% 1iii.7. Compute the prestress components according to the options
        if ~ischar(prestress)
            if isa(prestress.voigtVector,'function_handle')
                prestressVoigtVector = prestress.voigtVector(theta);
            else
                prestressVoigtVector = prestress.voigtVector;
            end
            if isPrestressOverDefinedSystem
                prestressVoigtVector = T2LC*prestressVoigtVector;
            end
            if strcmp(compPrestress,'1')
                sigma0(counterXi,counterEta) = prestressVoigtVector(1,1);
            elseif strcmp(compPrestress,'2')
                sigma0(counterXi,counterEta) = prestressVoigtVector(2,1);
            elseif strcmp(compPrestress,'12')
                sigma0(counterXi,counterEta) = prestressVoigtVector(3,1);
            end
        end
    	
        %% 1iii.8. Update counter for the lines in xi-direction
    	counterXi = counterXi + 1;
        
        %% 1iii.9. Update the xi-parametric coordinate
        xi = xi + dxi;
    end
  %% 1iv. Update counter for the lines in eta-direction
  counterEta = counterEta + 1;
  
  %% 1v. Update the eta-parametric coordinate
  eta = eta + deta;
end

%% 2. Write the coordinates into the individual arrays
X = S(:,:,1);
Y = S(:,:,2);
Z = S(:,:,3);

end