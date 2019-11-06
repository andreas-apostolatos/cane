function V_nc = computeNonConservativeLoadContributionToStiffnessKLShellNLinear(ub,vb,p,q,U,V,CP,mload,dir,CPd,int)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the contribution to the stiffness matrix from a non-consrvative
% externally applied load for the non-linear analysis of the Kirchhoff-Love
% shell
%
%   Input :
%   ub,vb : Load extension
%     p,q : B-Spline polynomial degrees in u-,v-direction
%     U,V : B-Spline knot vectors in u-,v-direction
%      CP : Set of Control Point coordinates and weights
%   mload : The load magnitude (or the load function)
%     dir : The direction of the load
%     CPd : The Control Point Coordinates of the deformed state
%     int : The quadrature rule
%
%  Output :
%    V_nc : The contribution to the non-Linear stiffness matrix for the
%           Kirchhoff-Love shell
%
% Function Layout :
%
% 0. Read input
%
% 1. On the integration
%
% 2. On the load application area
%
% 3. On the DoF numbering
%
% 4. Loop over all elements on which the load is applied
%
%    4i. Initialize element load vector
%
%   4ii. Loop over all Gauss points
%
%        4ii.1. Compute integration parameters
%
%        4ii.2. Compute the basis functions and its derivatives at the parametric location
%
%        4ii.3. Compute the covariant, and contravariant basis as well as the normal to the surface vector and the curvature change vector with respect to the reference configuration
%
%        4ii.4. Compute the the normal and the tangent to the surface boundary vectors
%
%        4ii.5. Compute the B-operator matrix for the rotations in the linear setting
%
%        4ii.6. Form the basis functions matrix
%
%        4ii.7. Compute the base vectors of the current configuration
%
%        4ii.8.. Compute the local force vector at the Gauss point
%
%        4ii.9. Add contribution to the local force vector
%
%  4iii. On the assembly phase
%
%% Function main body

%% 0. Read input

% Length of the knot vectors
mu = length(U);
mv = length(V);

% Number of control points in u,v-direction
nu = length(CP(:,1,1));
nv = length(CP(1,:,1));

% Number of local Control Points
nNode_loc = (p+1)*(q+1);
                    
% Number of local degrees of freedom
nDoF_loc = 3*nNode_loc;

% Number of global Control Points
nNode = nu*nv;

% Number of global degrees of freedom
nDoF = 3*nNode;

% if the fload parameter is numeric then assign directly its value
if isnumeric(mload)==1  
    m = mload;  
end

%% 1. On the integration

% Issue the Gauss points for the selected integration scheme:
if int.type == 0
   % Default scheme is the full gaussian quadrature element-wise (FGI)
   ngauss(1) = ceil((p+1)/2) ;
   ngauss(2) = ceil((q+1)/2) ;
elseif int.type == 1
    % Manual choice of the gauss points
    ngauss = int.ngaussLoad;
end

%% 2. On the load application area

% Find the span in u-direction where to apply the load
i1 = findKnotSpan(ub(1),U,nu);
if (isscalar(ub))
    % If its a scalar one Gauss point sufficies
    u = ub(1);    i2 = i1;
    ugauss = 1;   gwu = 1;
    mapu = 1;     is_on_u = 0;
else
    % If its not a scalar adjust the quadrature respectively
    i2 = findKnotSpan(ub(2),U,nu);
    if ub(2)~=U(mu)    
        i2=i2-1;    
    end
    ugauss = ngauss(1);
    is_on_u = 1;
end

% Find the span in v-direction where to apply the load
j1 = findKnotSpan(vb(1),V,nv);
if (isscalar(vb))
    % If its a scalar one Gauss point sufficies
    v = vb(1);    j2 = j1;
    vgauss = 1;   gwv = 1;
    mapv = 1;     is_on_u = 1;
else
    % If its not a scalar adjust the quadrature respectively
    j2 = findKnotSpan(vb(2),V,nv);
    if vb(2)~=V(mv)    
        j2=j2-1;    
    end
    vgauss = ngauss(2);
    is_on_u = 0;
end  

% Initialize global load vector
V_nc = zeros(nDoF,nDoF);


%% 3. On the DoF numbering

% dof array assigns two DoF (x,y) to every CP for the membrane:
% numbering follows CP: CP1->dof1,dof2 CP2->dof3,dof4
% Initialize the array of the degrees of freedom
dofs = zeros(nu,nv,3);

% Initialize counter
r=1;

% Loop over all the Control points
for cpj = 1:nv
    for cpi = 1:nu
        dofs(cpi,cpj,1)=r;
        dofs(cpi,cpj,2)=r+1;
        dofs(cpi,cpj,3)=r+2;
        
        % Update counter
        r=r+3;
    end
end

% Issue the Gauss points for the correct integration
[GPu,GWu] = gauss(ngauss(1));
[GPv,GWv] = gauss(ngauss(2));
            
%% 4. Loop over all elements on which the load is applied
for j = j1:j2
    for i = i1:i2
        % check if we are in a non-zero knot span
        if (U(i+1)~=U(i) && V(j+1)~=V(j))
            
            %% 4i. Initialize element load vector
            V_nc_el = zeros(nDoF_loc,nDoF_loc); 
            
            %% 4ii. Loop over all Gauss points
            for kv = 1:vgauss
                for ku = 1:ugauss
                    %% 4ii.1. Compute integration parameters
                    
                    % load extension in u-direction
                    if (isscalar(ub)==0)
                        % map the quadrature point into the knot span in
                        % u-direction
                        u = ( U(i+1)+U(i) + GPu(ku)*(U(i+1)-U(i)) )/2;
                        % compute the respective Jacobian determinant
                        mapu = (U(i+1)-U(i))/2;
                        % issue quadrature weight in u-direction
                        gwu = GWu(ku);
                    end
                    % load extension in v-direction
                    if (isscalar(vb)==0)
                        % map the quadrature point into the knot span in
                        % v-direction
                        v = ( V(j+1)+V(j) + GPv(kv)*(V(j+1)-V(j)) )/2;
                        % compute the respective Jacobian determinant
                        mapv = (V(j+1)-V(j))/2;
                        % issue quadrature weight in v-direction
                        gwv = GWv(kv);
                    end
                    % Compute the product of the Jacobian determinants in
                    % u,v-directions
                    map = mapu*mapv;
                    
                    % Compute the product of the quadrature weights 
                    gw = gwu*gwv;
                    
                    if isnumeric(mload)==0  
                        m = mload(p,i,u,U,q,j,v,V,CP);  
                    end

                    %% 4ii.2. Compute the basis functions and its derivatives at the parametric location
                    [R,dR] = computeNurbsBasisFunctionsAndFirstDerivatives2D(i,p,u,U,j,q,v,V,CP);
                    
                    %% 4ii.3. Compute the covariant, and contravariant basis as well as the normal to the surface vector and the curvature change vector with respect to the reference configuration
                    [g_covariant,g3,dA,H,~,~,g_contravariant,~,Bv] = computeCovariantBaseVectorsDerivativesAndMetrics2D(i,p,u,U,j,q,v,V,CP);
                   
                    %% 4ii.4. Compute the the normal and the tangent to the surface boundary vectors
                    [n,t] = computeNormalAndTangentVectorsToBSplineSurfaceBoundary(i,p,u,U,j,q,v,V,CP,is_on_u);           
                    
                    %% 4ii.5. Compute the B-operator matrix for the rotations in the linear setting
                    B_rotations = computeBOperatorMatrixForRotationsLinear(R,dR,g3,dA,g_covariant,H,g_contravariant,Bv,nNode_loc,n,t);
                        
                    % The applied moment tractions on the boundary that is a 2-dimensional entity on the curvilinear coordinate system
                    moment = zeros(2,1);
                    
                    if dir==1
                       moment(1,1) = m;
                       moment(2,1) = 0;
                    elseif dir==2
                       moment(1,1) = 0;
                       moment(2,1) = m;
                    end
                    
                    %% 4ii.6. Form the basis functions matrix
                    
                    % Initialize matrix
                    R_s = zeros(3,nNode_loc);
                    
                    % Assign the entries recursively
                    for iRs=1:nNode_loc
                        R_s(1,3*iRs-2) = R(iRs);
                        R_s(2,3*iRs-1) = R(iRs);
                        R_s(3,3*iRs) = R(iRs);
                    end
                   
                    
                    %% 4ii.7. Compute the base vectors of the current configuration
                    [g,~,~,~,~,~,~] = computeMetricsForKirchhoffLoveShellNonLinearCurrent(p,i,u,U,q,j,v,V,CPd);
                    
                    %% 4ii.8.. Compute the local force vector at the Gauss point
                    if is_on_u==1
                        V_nc_el_GP = B_rotations*moment*g(:,2)'*R_s*gw*map/norm(g(:,2));
                    else
                        V_nc_el_GP = B_rotations*moment*g(:,1)'*R_s*gw*map/norm(g(:,1));
                    end
                    
                    %% 4ii.9. Add contribution to the local force vector
                    V_nc_el = V_nc_el + V_nc_el_GP;
                end  
            end
            
            %% 4iii. On the assembly phase
    
            % Assemble the element coupling matrix to the global one via element
            % freedom tables

            % Element freedom table for the membrane part:
            % Initialize element freedom table
            dof_s = zeros(1,nDoF_loc);

            % Initialization of the counter for the element freedom table on the
            % structural side
            r = 1;

            % Relation global-local DoFs
            for cpj = j-q:j
                for cpi = i-p:i
                    dof_s(r)   = dofs(cpi,cpj,1);
                    dof_s(r+1) = dofs(cpi,cpj,2);
                    dof_s(r+2) = dofs(cpi,cpj,3);

                    % update counter
                    r = r + 3;
                end
            end

            % insert Me into M according to the algorithm under use:
            for mei = 1:nDoF_loc
                V_nc(dof_s(mei),1) = V_nc(dof_s(mei),1) + V_nc_el(mei,1);
            end  
            
        end
    end
end 


end

