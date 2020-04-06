function [Xir, Etar, CPr] = knotRefineBSplineSurface ...
    (p, Xi, q, Eta, CP, Rxi, Reta, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Perform knot insertion in a given NURBS surface in the 3D space. Source 
% reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag, Berlin 1995
% p. 72.
%
%    Input :
%      p,q : the polynomial degrees
%   Xi,Eta : the knot vectors
%       CP : control point coordinates and weights
% Rxi,Reta : vectors with knots to be inserted in U and V
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
%   Output :
%      CPr : the new set of Control Point coordinates and weights
% Xir,Etar : the new knot vectors
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the projective control points Pw
%
% 2. Knot refine the surface in the xi-direction
%
% 3. Knot refine the surface in the eta-direction
%
% 4. Transform the projective Control Points Qw to CPr
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Knot insertion for a B-Spline surface has been initiated \n\n');
    fprintf('Number of knots before knot insertion in xi-direction nxi = %d\n', length(Xi));
    fprintf('Number of knots after knot insertion in xi-direction nxi = %d\n', length(Xi) + length(Rxi));
    fprintf('Number of knots before knot insertion in eta-direction neta = %d\n', length(Eta));
    fprintf('Number of knots after knot insertion in eta-direction neta = %d\n', length(Eta) + length(Reta));
    fprintf('______________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Number of Control Points in xi-direction
nxi = length(CP(:, 1, 1));

% Number of Control Points in eta-direction
neta = length(CP(1, :, 1));

% Number of knots of the unrefined knot vector in xi-direction
mxi = nxi + p + 1;

% Number of knots of the unrefined knot vector in eta-direction
meta = neta + q + 1;

% Number of Control Points of the refined geometry in xi-direction
rxi = length(Rxi);
nxir = nxi + rxi;

% Number of Control Points of the refined geometry in eta-direction
reta = length(Reta);
netar = neta + reta;

% Initialize output arrays
Xir = zeros(1, mxi + rxi);
Etar = zeros(1, meta + reta);
CPr = zeros(nxir, netar, length(CP(1, 1, :)));
Pw = zeros(nxi, neta, length(CP(1, 1, :)));
Qw = zeros(nxir, netar, length(CP(1, 1, :)));

%% 1. Compute the projective control points Pw
for j = 1:neta
    for i = 1:nxi
        Pw(i, j, 1:3) = CP(i, j, 1:3)*CP(i, j, 4);
        Pw(i, j, 4)   = CP(i, j, 4);
    end
end

%% 2. Knot refine the surface in the xi-direction

% If the given vector of knots is not empty
if isempty(Rxi)
    Xir = Xi;
    Qw = Pw;
else
    a = findKnotSpan(Rxi(1), Xi, nxi);
    b = findKnotSpan(Rxi(rxi), Xi, nxi) + 1;

    for col = 1:neta
        for j = 1:a - p
            Qw(j, col, :) = Pw(j, col, :);   
        end
        for j = b - 1:nxi
            Qw(j + rxi, col, :) = Pw(j, col, :);
        end
    end
    
    for j = 1:a
        Xir(j) = Xi(j);
    end
    
    for j = b + p:mxi
        Xir(j+rxi) = Xi(j);
    end
    
    i = b + p - 1;
    k = i + rxi;
  
    for  j = rxi:-1:1
        while Rxi(j) <= Xi(i) && i > a
            for col = 1:neta
                Qw(k - p - 1, col, :) = Pw(i - p - 1, col, :);
            end
            Xir(k) = Xi(i);
            k = k-1;
            i = i-1;
        end
    
        for col = 1:neta
            Qw(k - p - 1, col, :) = Qw(k - p, col, :);
        end
        
        for l = 1:p
            ind = k - p + l;
            alpha = (Rxi(j) - Xir(k + l))/(Xi(i - p + l) - Xir(k + l));
            for col = 1:neta
                Qw(ind - 1, col, :) = alpha*Qw(ind - 1, col, :) + (1 - alpha)*Qw(ind, col, :);
            end
        end
        
        Xir(k) = Rxi(j);
        k = k - 1;
        
    end
    Pw = Qw;
end

%% 3. Knot refine the surface in the eta-direction

if isempty(Reta)
    Etar = Eta;
else
    a = findKnotSpan(Reta(1), Eta, neta);
    b = findKnotSpan(Reta(reta), Eta, neta) + 1;
  
    for row = 1:nxir
        for j = 1:a - q
            Qw(row, j, :) = Pw(row, j, :);
        end
        
        for j = b - 1:neta
            Qw(row, j + reta, :) = Pw(row, j, :);
        end
    end 
    for j = 1:a
        Etar(j) = Eta(j);  
    end
    
    for j = b + q:meta
        Etar(j + reta) = Eta(j);
    end

    i = b + q - 1;
    k = i + reta;
  
    for j = reta:-1:1
        while (Reta(j) <= Eta(i)) && (i > a)
            for row = 1:nxir
                Qw(row, k - q - 1, :) = Pw(row, i - q - 1, :);
            end
            
            Etar(k) = Eta(i);
            k = k - 1;
            i = i - 1;
        end

        for row = 1:nxir
            Qw(row, k - q - 1, :) = Qw(row, k - q, :);
        end
        for l = 1:q
            ind = k - q + l;
            alpha = (Reta(j) - Etar(k + l))/(Eta(i - q + l) - Etar(k + l));
            for row = 1:nxir
                Qw(row, ind - 1, :) = alpha*Qw(row, ind - 1, :) + ....
                    (1 - alpha)*Qw(row, ind, :);
            end
        end
        Etar(k) = Reta(j);
        k = k - 1;
    end
end

%% 4. Transform the projective Control Points Qw to CPr
for j = 1:length(Qw(1, :, 1))
    for i = 1:length(Qw(:, 1, 1))
        CPr(i, j, 1:3) = Qw(i, j, 1:3)/Qw(i, j, 4);
        CPr(i, j, 4) = Qw(i, j, 4);
    end
end

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Knot insertion took %d seconds \n\n', computationalTime);
    fprintf('_____________________Knot Insertion Ended_____________________\n');
    fprintf('##############################################################\n\n\n');
end

end
