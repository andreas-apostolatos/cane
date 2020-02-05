function [tanMtxLinEl,tanMtxNLinEl,massMtxEl,FBodyEl] = ...
    computeFEMVMSStabElTangentStiffMtxMassMtxLoadVctNLinear4NSE2D...
    (x,y,z,t,upEl,uMeshALEEL,dN,computeBodyForces,parameters,h,...
    propFldDynamics,isAnalysis3D)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the tangent stiffness matrix, the mass matrix and the body force
% vector at the element level, corresponding to the variational muiltiscale 
% stabilized (VMS) finite element approximation of the transient
% Navier-Stokes problem in 2D using the linear triangle (CST) element.
%
% Reference :
%
% Ramon Codina, "A stabilized finite element method for generalized 
% stationary incompressible flows", Computer Methods in Applied Mechanics 
% and Engineering, Vol. 190 (2001), 2681-2706
%
% Implementation :
%
% KRATOS opensourse project, Riccardo Rossi
%
%           Input :
%        analysis : Information on the analysis
%                       .type : analysis type
%         x,y,z,t : Chorochronical location of the Gauss Point in the physical
%                   space-time interval
%            upEl : The nodal solution vector at the element level from the
%                   previous nonlinear iteration step
%      uMeshALEEL : The mesh velocity at the element level
%              dN : matrix containing the basis functions and their 
%                   derivatives in the physical space page-wise for all 
%                   elements as follows,
%                               | N1 dN1/dx dN1/dy | 
%                          dN = | N2 dN2/dx dN2/dy |
%                               | N3 dN3/dx dN3/dy |
%      parameters : Technical parameters of the flow problem :
%                            .rho : flow density
%                            .nue : flow dynamic viscosity
%               h : The minimum element edge size
% propFldDynamics : Transient analysis parameters :
%                           .scheme : The time integration method
%                        .alphaBeta : (parameter for the Bossak scheme)
%                            .gamma : (parameter for the Bossak scheme)
%                               .T0 : Start time of the simulation
%                             .TEnd : End time of the simulation
%                               .nT : Number of time steps
%                               .dt : Time step
%    isAnalysis3D : Flag on whether the analysis is 3D
%
%          Output :
%     tanMtxLinEl : The element linearized stiffness matrix
%    tanMtxNLinEl : The element nonlinear stiffness matrix
%       massMtxEl : The element mass matrix
%         FBodyEl : The element body force vector
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the velocity basis functions matrix
%
% 2. Compute the convection velocity
%
% 3. Compute the stabilization parameters
%
% 4. Compute the derivatives of the velocity field
%
% 5. Compute the divergence of the velocity field
%
% 6. Compute the pressure basis functions matrix
%
% 7. Compute the B-operator matrix for the convective term
%
% 8. Compute terms for the tangent matrix of the convection term
%
% 9. Compute the tangent matrix for the convection term (SOMETHING DOES NOT FIT HERE FOR THE 3D IMPLMENTATION)
%
% 10. Compute the convection-tangent matrix
%
% 11. Compute the B-operator matrix for the pressure grandient term
%
% 12. Compute the B-operator matrix for the diffusion term
%
%% Function main body

%% 0. Read input

isBodyForceNonZero = false;
% b = computeBodyForces(x,y,z,t);
% if norm(b) == 0
%     isBodyForceNonZero = false;
% end

% Compute the number of nodes and DOFs in the element level
if isAnalysis3D
    noNodesEl = 4;
    noDOFsNode = 4;
else
    noNodesEl = 3;
    noDOFsNode = 3;
end

% Number of DOFs in the element level
noDOFsEl = noDOFsNode*noNodesEl;

% Number of elements in total
noElmnts = size(dN,1);

% Initialize the velocities basis function matrix
RVelocities = zeros(noElmnts,3,noDOFsEl);

% Initialize a zero matrix
zero = zeros(noElmnts,1,1);

% Initialize the B-operator matrix related to each row of the tangent
% convection matrix
xiBoperatorTangentConvection = zeros(noElmnts,noDOFsEl,noDOFsEl);

% Compute the relative (ALE) velocity
relVelocity = (upEl - uMeshALEEL)';

% Initialize the velocities basis function matrix only if the source term
% (body forces) are non-zero as well as the matrix due to the stabilization 
% of the body forces-convection term which goes on the left-hand side
if isBodyForceNonZero
    Bc = zeros(noDOFsEl,noDOFsEl);
end

%% 1. Compute the velocity basis functions matrix
% for i=1:noNodesEl
%     RVelocities(1,3*i-2) = dN(i,1);
%     RVelocities(2,3*i-1) = dN(i,1);   
% end
for i = 1:noNodesEl
    RVelocities(:,1,noDOFsNode*i-noDOFsNode+1) = dN(:,i,1);
    RVelocities(:,2,noDOFsNode*i-noDOFsNode+2) = dN(:,i,1);
    if isAnalysis3D
        RVelocities(:,3,noDOFsNode*i-noDOFsNode+3) = dN(:,i,1);
    end
end

%% 2. Compute the convection velocity
% convectionVelocity = RVelocities*relVelocity;
convectionVelocity = pmtimes(RVelocities,relVelocity);
euclideanNorm = @(vector) sqrt(vector(:,1,1).^2 + vector(:,2,1).^2 + vector(:,3,1).^2)';

%% 3. Compute the stabilization parameters

% Stabilization constants
CI = 4;
Ct = 4;

% Compute the stabilization parameter for the momentum equation
if strcmp(propFldDynamics.timeDependence,'STEADY_STATE')
    tauM = (2*euclideanNorm(convectionVelocity)'./h +...
    ((CI*parameters.nue)./(h.^2))).^(-1);
elseif strcmp(propFldDynamics.timeDependence,'TRANSIENT')
    tauM = ( Ct/propFldDynamics.dt + 2*euclideanNorm(convectionVelocity)'./h +...
    ((CI*parameters.nue)./(h.^2)) ).^(-1);
else
    error('Neither steady-state nor transient analysis is selected, see input file');
end

% Compute the stabilization parameter for the continuity equation
tauC = (parameters.nue + .5*h.*euclideanNorm(convectionVelocity)');

%% 4. Compute the derivatives of the velocity field
% dudx(1,3*i-2) = dN(i,2);
% dudy(1,3*i-2) = dN(i,3);
% dvdx(1,3*i-1) = dN(i,2);
% dvdy(1,3*i-1) = dN(i,3);

% Nonsingleton matrix containing the derivatives of the field wrt
% x-coordinate
dudx = reshape(ptranspose(phorzcat(dN(:,:,2),zeros(noElmnts,noDOFsNode,noDOFsNode-1))),noElmnts,1,noDOFsEl);
dvdx = phorzcat(zero,dudx(:,1,1:end-1));
if isAnalysis3D
    dwdx = phorzcat(zero,zero,dudx(:,1,1:end-2));
end

% Nonsingleton matrix containing the derivatives of the field wrt
% y-coordinate
dudy = reshape(ptranspose(phorzcat(dN(:,:,3),zeros(noElmnts,noDOFsNode,noDOFsNode-1))),noElmnts,1,noDOFsEl);
dvdy = phorzcat(zero,dudy(:,1,1:end-1));
if isAnalysis3D
    dwdy = phorzcat(zero,zero,dudy(:,1,1:end-2));
end

% Nonsingleton matrix containing the derivatives of the field wrt
% z-coordinate
if isAnalysis3D
    dudz = reshape(ptranspose(phorzcat(dN(:,:,4),zeros(noElmnts,noDOFsNode,noDOFsNode-1))),noElmnts,1,noDOFsEl);
    dvdz = phorzcat( zero, dudz(:,1,1:end-1));
    dwdz = phorzcat( zero, zero, dudz(:,1,1:end-2));
end

%% 5. Compute the divergence of the velocity field
% dVelocities(1,3*i-2) = dN(i,2);
% dVelocities(1,3*i-1) = dN(i,3); 
dVelocities = dudx + dvdy;
if isAnalysis3D
    dVelocities = dVelocities + dwdz;
end

%% 6. Compute the pressure basis functions matrix
% RPressure(1,3*i) = dN(i,1);
RPressure = reshape(ptranspose(phorzcat(zeros(noElmnts,noDOFsNode,noDOFsNode-1),dN(:,:,1))),noElmnts,1,noDOFsEl);

%% 7. Compute the B-operator matrix for the convective term
% BoperatorConvection(1,3*i-2) = convectionVelocity(1,1)*dN(i,2) + convectionVelocity(2,1)*dN(i,3);
% BoperatorConvection(2,3*i-1) = BoperatorConvection(1,3*i-2);
BoperatorConvection = zeros(noElmnts,3,noDOFsEl);
if isAnalysis3D
    BoperatorConvection(:, 1, 1:noDOFsNode:end) = ptranspose(pstimes(dN(:,:,2),convectionVelocity(:,1,1)) + ...
        pstimes(dN(:,:,3), convectionVelocity(:,2,1)) + pstimes(dN(:,:,4), convectionVelocity(:,3,1)));
else
    BoperatorConvection(:, 1, 1:noDOFsNode:end) = ptranspose(pstimes(dN(:,:,2),convectionVelocity(:,1,1)) + ...
        pstimes(dN(:,:,3), convectionVelocity(:,2,1)));
end
BoperatorConvection(:,2,2:noDOFsNode:end) = BoperatorConvection(:,1,1:noDOFsNode:end);
if isAnalysis3D
    BoperatorConvection(:,3,3:noDOFsNode:end) = BoperatorConvection(:,1,1:noDOFsNode:end);
end

%% 8. Compute terms for the tangent matrix of the convection term
% xiBoperatorTangentConvection(3*j-2,3*i-2) = dN(j,1)*dN(i,2);
% xiBoperatorTangentConvection(3*j-1,3*i-2) = dN(j,1)*dN(i,3);
% etaBoperatorTangentConvection(3*j-2,3*i-1) = dN(j,1)*dN(i,2);
% etaBoperatorTangentConvection(3*j-1,3*i-1) = dN(j,1)*dN(i,3);

% Derivative wrt xi-direction
xiBoperatorTangentConvection(:,1:noDOFsNode:end,1:noDOFsNode:end) = ...
    pmtimes(dN(:,:,1),ptranspose(dN(:,:,2)));
xiBoperatorTangentConvection(:,2:noDOFsNode:end,1:noDOFsNode:end) = ...
    pmtimes(dN(:,:,1),ptranspose(dN(:,:,3)));
if isAnalysis3D
    xiBoperatorTangentConvection(:,3:noDOFsNode:end,1:noDOFsNode:end) = ...
        pmtimes(dN(:,:,1),ptranspose(dN(:,:,4)));
end

% Derivative wrt eta-direction
etaBoperatorTangentConvection = phorzcat(zeros(noElmnts,noDOFsEl,1),...
    xiBoperatorTangentConvection(:,:,1:end-1));

% Derivative wrt zeta-direction
if isAnalysis3D
    zetaBoperatorTangentConvection = phorzcat(zeros(noElmnts,noDOFsEl,2),...
        xiBoperatorTangentConvection(:,:,1:end-2));
end

% The contribution from the body forces is missing,
% Bc(3*j-2,3*i-2) = - bodyForcesRearranged(1)*dN(i,1)*dN(j,2);
% Bc(3*j-2,3*i-1) = - bodyForcesRearranged(1)*dN(i,1)*dN(j,3);
% Bc(3*j-1,3*i-2) = - bodyForcesRearranged(2)*dN(i,1)*dN(j,2);
% Bc(3*j-1,3*i-1) = - bodyForcesRearranged(2)*dN(i,1)*dN(j,3);

%% 9. Compute the tangent matrix for the convection term (SOMETHING DOES NOT FIT HERE FOR THE 3D IMPLMENTATION)
% CTangentBOperator = [relVelocity'*xiBoperatorTangentConvection' 
%                      relVelocity'*etaBoperatorTangentConvection' 
%                      zeros(1,noDOFsEl)];
if isAnalysis3D
    CTangentBOperator = pvertcat(pmtimes(ptranspose(relVelocity),ptranspose(xiBoperatorTangentConvection)),...
                                 pmtimes(ptranspose(relVelocity),ptranspose(etaBoperatorTangentConvection)),...
                                 pmtimes(ptranspose(relVelocity),ptranspose(zetaBoperatorTangentConvection)));
else
    CTangentBOperator = pvertcat(pmtimes(ptranspose(relVelocity),ptranspose(xiBoperatorTangentConvection)),...
                                 pmtimes(ptranspose(relVelocity),ptranspose(etaBoperatorTangentConvection)),...
                                 zeros(noElmnts,1,noDOFsEl));
end

%% 10. Compute the convection-tangent matrix
% KNLinearEl = RVelocities'*CTangentBOperator;
tanMtxNLinEl = pmtimes(ptranspose(RVelocities),CTangentBOperator);
                         
%% 11. Compute the B-operator matrix for the pressure grandient term
% gradientPressure(1,3*i) = dN(i,2);
% gradientPressure(2,3*i) = dN(i,3);
if isAnalysis3D
    gradientPressure = phorzcat(zeros(noElmnts,3,3),pvertcat(dudx(:,1,1:end-3),...
        dudy(:,1,1:end-3),dudz(:,1,1:end-3)));
else
    gradientPressure = phorzcat(zeros(noElmnts,3,2),[dudx(:,1,1:end-2),...
        dudy(:,1,1:end-2),zeros(noElmnts,1,7)]);
end

%% 12. Compute the B-operator matrix for the diffusion term
% De = parameters.nue*(dudx'*dudx + dudy'*dudy + dvdx'*dvdx + dvdy'*dvdy);
if isAnalysis3D
    BVelocityDivergence = pvertcat(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz);
else
    BVelocityDivergence = pvertcat(dudx,dudy,dvdx,dvdy);
end
De = parameters.nue*pmtimes(ptranspose(BVelocityDivergence),BVelocityDivergence);

%% 13. Compute the B-operator matrix for the velocity-pressure coupling
% Pe = - dVelocities'*RPressure;
Pe = - pmtimes(ptranspose(dVelocities),RPressure);

%% 14. Compute the B-operator matrix for convection term
% Ce = RVelocities'*BoperatorConvection;
Ce = pmtimes(ptranspose(RVelocities),BoperatorConvection);

%% 15. Compute the stabilization matrices for the stiffness
%% 15i. Matrices related to convection v.Gradient(v)

% Convection-convection stabilization:
% Cc = BoperatorConvection'*BoperatorConvection;
Cc = pmtimes(ptranspose(BoperatorConvection),BoperatorConvection);

% Convection-pressure stabilization:
% Cp = BoperatorConvection'*gradientPressure;
Cp = pmtimes(ptranspose(BoperatorConvection),gradientPressure);

%% 15ii. Matrices related to pressure gradient Nabla(p)

% Pressure-convection stabilization: (Continuity equation)
% Pc = Cp';
Pc = ptranspose(Cp);

% Pressure-pressure stabilization: (Continuity equation)
% Pp = gradientPressure'*gradientPressure;
Pp = pmtimes(ptranspose(gradientPressure),gradientPressure);

%% 15iii. Matrix related to the divergence-free condition (actual elimination of the saddle-point problem)
% Div = dVelocities'*dVelocities;
Div = pmtimes(ptranspose(dVelocities),dVelocities); % (Continuity equation)

%% 16. Compute the stabilization matrices on the mass
%% 16i. Mass-convection stabilization
% Mc = BoperatorConvection'*RVelocities;
Mc =  pmtimes(ptranspose(BoperatorConvection),RVelocities);

%% 16ii. Mass-pressure stabilization
% Mp = gradientPressure'*RVelocities;
Mp = pmtimes(ptranspose(gradientPressure), RVelocities);

%% 17. Compute the stabilization vectors of the right-hand side (stabilizing the right-hand side linear functional)
if isBodyForceNonZero
    % Stabilization due to pressure gradient Nabla(p): (Continuity equation)
    Fp = gradientPressure';
end
       
%% 18. Compute the element stiffness matrix
if isBodyForceNonZero
    % Linear stiffness matrix needed for the residual computation
    tanMtxLinEl = De + Ce + Pe - ptranspose(Pe) + pstimes((Cc + Cp + Pc + Pp + Bc),tauM) + pstimes(Div,tauC);
else
    % Linear stiffness matrix needed for the residual computation
    tanMtxLinEl = De + Ce + Pe - ptranspose(Pe) + pstimes((Cc + Cp + Pc + Pp),tauM) + pstimes(Div,tauC);
end
% clear Ce Pe Cc Pc Pp Div tauC

%% 19. Compute the element stabilized mass matrix
if strcmp(propFldDynamics.timeDependence,'TRANSIENT')
    massMtxEl = pmtimes(ptranspose(RVelocities),RVelocities) + pstimes(Mc + Mp, tauM);
else
    massMtxEl = 'undefined';
end

%% 20. Compute the stabilized force vector due to source at the Gauss point
if isBodyForceNonZero
    FBodyEl = (RVelocities' + tauM*Fp)*b(1:2,1);
else
	FBodyEl = zeros(noDOFsEl,1);
end

end
