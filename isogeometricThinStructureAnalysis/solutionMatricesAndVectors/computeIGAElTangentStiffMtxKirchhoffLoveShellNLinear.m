function [KTel, FIEl] = ...
    computeIGAElTangentStiffMtxKirchhoffLoveShellNLinear ...
    (p, q, dR, ACovariant, dACovariant, A3Tilde, aCovariant, ...
    daCovariant, a3Tilde, Dm, Db)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix and the internal residual load
% vector at the element level.
% 
%        Input : 
%          p,q : The polynomial degrees of the B-Spline patch
%           dR : The B-Spline basis functions, their first and second
%                derivatives
%   ACovariant : The covariant base vectors of the reference configuration
%  dACovariant : The derivatives of the covariant base vectors
%      A3Tilde : The not normalized surface normal vector of the reference
%                configuration
%   aCovariant : The covariant base vectors of the current configuration
%  daCovariant : The derivatives of the covariant base vectors of the 
%                current configuration
%      a3Tilde : The not normalized surface normal vector of the current
%                configuration
%           Dm : The membrane material stiffness matrix
%           Db : The bending material stiffness matrix
%
%       Output :
%         KTel : Element tangent stiffness matrix
%         FIEl : Element internal residual load vector
%
% Function layout :
%
% THIS FUNCTION NEEDS IMPROVEMENT WITH RESPECT TO THE COMPUTATIONS (EFFICIENCY)
%
%% Function main body

%% 0. Read input

% Number of DoFs at the element level
ndof = (p+1)*(q+1)*3;

% Initialize element tangent stiffness matrix
Kem = zeros(ndof,ndof);
Keb = zeros(ndof,ndof);

% Initialize element internal residual load vector
FInternalElementMembrane = zeros(ndof,1);
FInternalElementBending = zeros(ndof,1);

% Initialize auxuiliary arrays
dg = zeros(3,2,ndof);
dg3Tilde = zeros(3,ndof);
g3Tildedg3Tilde = zeros(ndof);
g3dg3lg3_3 = zeros(ndof);
dn = zeros(3,ndof);

% Initialize arrays on the first variation of the strain and the curvature
% w.r.t. the DoFs
dECartesian = zeros(3,ndof);
dKCartesian = zeros(3,ndof);

%% 1. Compute metrics for the reference and the current configuration

% For the reference configuration :
% _________________________________

% Compute the covariant metric coefficients
GabCovariant = zeros(2,2);
GabCovariant(1,1) = ACovariant(:,1)'*ACovariant(:,1);
GabCovariant(1,2) = ACovariant(:,1)'*ACovariant(:,2);
GabCovariant(2,1) = GabCovariant(1,2);
GabCovariant(2,2) = ACovariant(:,2)'*ACovariant(:,2);

% Compute the contravariant basis
GContravariant = GabCovariant\ACovariant';
GContravariant = GContravariant';

% Compute the surface unit normal vector
G3 = A3Tilde/norm(A3Tilde);

% Compute the local Cartesian basis
eLC = computeLocalCartesianBasis4BSplineSurface(ACovariant,GContravariant);

% Compute the transformation matrix from the contravariant basis to the
% local Cartesian one
TFromContraToLC = computeTFromContra2LocalCartesian4VoigtStrainIGAKLShell(eLC,GContravariant);

% Compute the curvature in Voigt notation
BV = dACovariant'*G3;

% For the current configuration :
% _______________________________

% Compute the covariant metric coefficients
gabCovariant = zeros(2,2);
gabCovariant(1,1) = aCovariant(:,1)'*aCovariant(:,1);
gabCovariant(1,2) = aCovariant(:,1)'*aCovariant(:,2);
gabCovariant(2,1) = gabCovariant(1,2);
gabCovariant(2,2) = aCovariant(:,2)'*aCovariant(:,2);

% Normalize it to unit length
g3 = a3Tilde/norm(a3Tilde);

% Compute the curvature in Voigt notation
bV = daCovariant'*g3;

%% 2. Compute non-linear strain and curvature in the Cartesian coordinate system

% Get the non-linear strain vector [E11,E22,E12] referred to curvilinear 
% coordinate system
ECurvilinear(1,1) = 1/2*(gabCovariant(1,1)-GabCovariant(1,1));
ECurvilinear(2,1) = 1/2*(gabCovariant(2,2)-GabCovariant(2,2));
ECurvilinear(3,1) = 1/2*(gabCovariant(1,2)-GabCovariant(1,2));

% Transformation of the strain vector [E11,E22,2*E12] from the curvilinear
% to the Cartesian coordinate system
ELocalCartesian = TFromContraToLC*ECurvilinear;

% Get the non-linear curvature vector [K11,K22,K12] referred to curvilinear 
% coordinate system
KCurvilinear = - (bV - BV);

% Transformation of the curvature vector [K11,K22,K12] from the curvilinear
% to the Cartesian coordinate system
KLocalCartesian = TFromContraToLC*KCurvilinear;

%% 3. Compute the first variation of the axial strain and change in the curvature with respect to the DoFs

% Compute first variation of strain and curvature w.r.t. the DoFs
for r = 1:ndof
    %% 3i. Get local node number kr and dof direction dirr
    kr = ceil(r/3);
    dirr = r-3*(kr-1);
    dg(dirr,1,r) = dR(kr,2);
    dg(dirr,2,r) = dR(kr,4);
  
    %% 3ii. Compute strain variation
    
    % Compute strain variation on the curvilinear system
    dEContravariant(1,1) = dR(kr,2)*aCovariant(dirr,1);
    dEContravariant(2,1) = dR(kr,4)*aCovariant(dirr,2);
    dEContravariant(3,1) = 1/2*(dR(kr,2)*aCovariant(dirr,2) + aCovariant(dirr,1)*dR(kr,4));
    
    % Compute strain variation on the Cartesian system
    dECartesian(:,r) = TFromContraToLC*dEContravariant;
  
    %% 3iii. Compute curvature variation
    
    % Compute the derivative of the normal to the surface vector
    dg3Tilde(1,r) = dg(2,1,r)*aCovariant(3,2)-dg(3,1,r)*aCovariant(2,2) + aCovariant(2,1)*dg(3,2,r)-aCovariant(3,1)*dg(2,2,r);
    dg3Tilde(2,r) = dg(3,1,r)*aCovariant(1,2)-dg(1,1,r)*aCovariant(3,2) + aCovariant(3,1)*dg(1,2,r)-aCovariant(1,1)*dg(3,2,r);
    dg3Tilde(3,r) = dg(1,1,r)*aCovariant(2,2)-dg(2,1,r)*aCovariant(1,2) + aCovariant(1,1)*dg(2,2,r)-aCovariant(2,1)*dg(1,2,r);
    
    % Compute product g3Tilde*dg3Tilde
    g3Tildedg3Tilde(r) = a3Tilde(1)*dg3Tilde(1,r)+a3Tilde(2)*dg3Tilde(2,r)+a3Tilde(3)*dg3Tilde(3,r);
    
    % Compute vector g3*dg3/lg3_3
    g3dg3lg3_3(r) = g3Tildedg3Tilde(r)/norm(a3Tilde)^3;
    
    % Compute the derivative of the unit normal to the surface vector
    dn(1,r) = dg3Tilde(1,r)/norm(a3Tilde) - a3Tilde(1)*g3dg3lg3_3(r);
    dn(2,r) = dg3Tilde(2,r)/norm(a3Tilde) - a3Tilde(2)*g3dg3lg3_3(r);
    dn(3,r) = dg3Tilde(3,r)/norm(a3Tilde) - a3Tilde(3)*g3dg3lg3_3(r);
    
    % Compute curvature variation on the curvilinear system
    dKCurvilinear(1,1) = -(dR(kr,3)*g3(dirr) + daCovariant(1,1)*dn(1,r)+daCovariant(2,1)*dn(2,r)+daCovariant(3,1)*dn(3,r));
    dKCurvilinear(2,1) = -(dR(kr,6)*g3(dirr) + daCovariant(1,2)*dn(1,r)+daCovariant(2,2)*dn(2,r)+daCovariant(3,2)*dn(3,r));
    dKCurvilinear(3,1) = -(dR(kr,5)*g3(dirr) + daCovariant(1,3)*dn(1,r)+daCovariant(2,3)*dn(2,r)+daCovariant(3,3)*dn(3,r));
    
    % Compute curvature variation on the Cartesian system
    dKCartesian(:,r) = TFromContraToLC*dKCurvilinear;
end

%% 4. Compute the second variation of the axial strain and change in the curvature with respect to the DoFs

% Initialize arrays on the second variation of the strain and the curvature
% w.r.t. the DoFs
ddECartesian = zeros(3,ndof,ndof);
ddKCartesian = zeros(3,ndof,ndof);

% Compute second variation of strain and curvature w.r.t. dofs ur and us

% Loop over r DoFs
for r = 1:ndof
    %% 4i. Get local node number kr and dof direction dirr for DoF ur
    kr = ceil(r/3);
    dirr = r-3*(kr-1);
    
    % Loop over s DoFs for each r DoF
    for s = 1:r
        %% 4i.1. Get local node number ks and dof direction dirs for DoF us
        ks = ceil(s/3);
        dirs = s-3*(ks-1);
    
        %% 41.2. Compute strain variation
        ddECurvilinear = zeros(3,1);
        
        % Compute strain variation on the curvilinear system
        if (dirr==dirs)
            ddECurvilinear(1) = dR(kr,2)*dR(ks,2);
            ddECurvilinear(2) = dR(kr,4)*dR(ks,4);
            ddECurvilinear(3) = 1/2*(dR(kr,2)*dR(ks,4)+dR(kr,4)*dR(ks,2));
        end

        % Compute strain variation on the Cartesian system
        ddECartesian(:,r,s) = TFromContraToLC*ddECurvilinear;
    
        %% 4i.3. Compute curvature variation
        
        % Initialize the second derivative of the normal to the surface vector
        ddg3 = zeros(3,1);
        
        % Compute the second derivative of the normal to the surface vector
        dirt = 6-dirr-dirs;
        ddir = dirr-dirs;
        if ddir==-1
            ddg3(dirt) =  dR(kr,2)*dR(ks,4)-dR(ks,2)*dR(kr,4);
        elseif ddir==2
            ddg3(dirt) =  dR(kr,2)*dR(ks,4)-dR(ks,2)*dR(kr,4);
        elseif ddir==1
            ddg3(dirt) = -dR(kr,2)*dR(ks,4)+dR(ks,2)*dR(kr,4);
        elseif (ddir==-2)
            ddg3(dirt) = -dR(kr,2)*dR(ks,4)+dR(ks,2)*dR(kr,4);
        end
        
        % Compute auxuliary variables C and D
        C = - (ddg3'*a3Tilde + dg3Tilde(:,r)'*dg3Tilde(:,s))/norm(a3Tilde)^3;
        D = 3*g3Tildedg3Tilde(r)*g3Tildedg3Tilde(s)/norm(a3Tilde)^5;
        
        % Compute the second derivative of the unit surface normal vector
        ddn = ddg3/norm(a3Tilde) - g3dg3lg3_3(s)*dg3Tilde(:,r) - g3dg3lg3_3(r)*dg3Tilde(:,s) + C*a3Tilde + D*a3Tilde;
        
        % On the second variation of the curvature in the Curvilinear
        % coordinate system
        ddKCurvilinear(1,1) = -(dR(kr,3)*dn(dirr,s) + dR(ks,3)*dn(dirs,r)...
                      + daCovariant(1,1)*ddn(1) + daCovariant(2,1)*ddn(2) + daCovariant(3,1)*ddn(3));
        ddKCurvilinear(2,1) = -(dR(kr,6)*dn(dirr,s) + dR(ks,6)*dn(dirs,r)...
                      + daCovariant(1,2)*ddn(1) + daCovariant(2,2)*ddn(2) + daCovariant(3,2)*ddn(3));
        ddKCurvilinear(3,1) = -(dR(kr,5)*dn(dirr,s) + dR(ks,5)*dn(dirs,r)...
                      + daCovariant(1,3)*ddn(1) + daCovariant(2,3)*ddn(2) + daCovariant(3,3)*ddn(3));
    
        
        % On the second variation of the curvature in the Cartesian
        % coordinate system
        ddKCartesian(:,r,s) = TFromContraToLC*ddKCurvilinear;
    end
end

%% 5. Compute the element membrane and the bending tangent stiffness matrices

%% 5i. Compute the first and the second variation of the force and the moment w.r.t. the DoFs
for r = 1:ndof
    % Compute the first variation of the axial force w.r.t. the DoFs ur 
    % on the Cartesian coordinate system: N_ca = Dm*E_ca
    NCartesian = Dm*ELocalCartesian;
    
    % Compute the first variation of the normal force w.r.t. the DoFs ur 
    % on the Cartesian coordinate system: M_ca = Db*K_ca
    MCartesian = Db*KLocalCartesian;
    
    % Compute the second variation of the axial force w.r.t. the DoFs ur
    % and us on the Cartesian coordinate system: dN_ca = Dm*dE_ca(:,r) 
    dNCartesian = Dm*dECartesian(:,r);
    
    % Compute the second variation of the moment w.r.t. the DoFs ur
    % on the Cartesian coordinate system: dM_ca = Db*dK_ca(:,r)
    dMCartesian = Db*dKCartesian(:,r);
    
    %% 5ii. Compute the membrane and the bending tangent stiffness matrices
    for s = 1:r
        % Compute membrane stiffness
        Kem(r,s) = dNCartesian'*dECartesian(:,s) + NCartesian'*ddECartesian(:,r,s);
        
        % Compute bending stiffness
        Keb(r,s) = dMCartesian'*dKCartesian(:,s) + MCartesian'*ddKCartesian(:,r,s);
    end
    
    % Compute the residual forces from the membrane and the bending parts
    FInternalElementMembrane(r,1) = - NCartesian'*dECartesian(:,r);
    FInternalElementBending(r,1) = - MCartesian'*dKCartesian(:,r);
end

%% 6. Compute the element tangent stiffness matrix and element internal residual load vector

% Compute element tangent stiffness matrix
KTel = Kem + Keb;

% Compute element internal residual load vector
FIEl = FInternalElementMembrane + FInternalElementBending;

end
