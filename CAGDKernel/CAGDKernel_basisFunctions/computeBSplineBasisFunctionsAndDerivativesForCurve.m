function dN = computeBSplineBasisFunctionsAndDerivativesForCurve(knotSpanIndex,p,xi,Xi,nDeriv)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Evaluates the B-Spline basis functions and its up to n-th derivatives 
% using the Cox De-Boor recursion formula at parameter xi for the 
% parametrization of a curve in 2D/3D. Source reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag, Berlin 1995
% p. 72.
%
%         Input :
% knotSpanIndex : The span index where xi is contained
%             p : The polynomial degree of the curve
%            xi : The coordinate in the unit interval where base vectors are to be
%                 evaluate
%            Xi : The knot vector of the B-Spline basis functions
%        nDeriv : The number of derivatives to be computed
%
%        Output :
%            dN : matrix for which the element dN(i,j) contains the (j-1)-
%                 th derivatives of the i-th non identically zero basis 
%                 function of the knot span with index knotSpanIndex. The 
%                 range of i and j indices is i=1,...,n+1 and j=1,...,p+1
%
% Function layout :
%
% 0. Read input
%
% 1. Save necessary values
%
% 2. Load the basis functions and compute the derivatives
%
% 3. Multiply through by the correct factors and compute the derivatives
%
%% Function main body

%% 0. Read input

left = zeros(1,p+1);
right = zeros(1,p+1);
ndxi = zeros(p+1,p+1);
dN = zeros(p+1,nDeriv+1);
a = zeros(2,p+1);

%% 1. Save necessary values

ndxi(1,1) = 1;
for j = 1:p
    left(j+1) = xi - Xi(knotSpanIndex+1-j);
    right(j+1) = Xi(knotSpanIndex+j) - xi;
    saved = 0;
    
    for r = 0:j-1
        % Lower triangle
        ndxi(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = ndxi(r+1,j)/ndxi(j+1,r+1);
        
        % Upper triangle
        ndxi(r+1,j+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    ndxi(j+1,j+1) = saved;
end

%% 2. Load the basis functions and compute the derivatives

for j = 0:p
    dN(j+1,1) = ndxi(j+1,p+1);
end

% compute derivatives
for r = 0:p          
    % loop over function index
    s1 = 0;
    
    % alternate rows in array a
    s2 = 1;                
    a(1,1) = 1;
    
    % loop to compute 1st derivative                      
    for k = 1:nDeriv
        d = 0;
        rk = r-k;
        pk = p-k;
        if r >= k
            a(s2+1,1) = a(s1+1,1)/ndxi(pk+2,rk+1);
            d = a(s2+1,1)*ndxi(rk+1,pk+1);
        end
        if (rk >= -1)
            j1 = 1;
        else 
            j1 = -rk;
        end
        if (r-1) <= pk
            j2 = k-1;
        else 
            j2 = p-r;
        end
        for j = j1:j2
            a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndxi(pk+2,rk+j+1);
            d = d + a(s2+1,j+1)*ndxi(rk+j+1,pk+1);
        end
        if (r <= pk)
            a(s2+1,k+1) = -a(s1+1,k)/ndxi(pk+2,r+1);
            d = d + a(s2+1,k+1)*ndxi(r+1,pk+1);
        end
        dN(r+1,k+1) = d;
        j = s1;
        s1 = s2;
        
        % switch rows
        s2 = j;             
     end
end
      
%% 3. Multiply through by the correct factors and compute the derivatives
r = p;
for k = 1:nDeriv
    for j = 0:p
        dN(j+1,k+1) = dN(j+1,k+1)*r;
    end
    r = r*(p-k);
end

%% 4. Appendix

end
