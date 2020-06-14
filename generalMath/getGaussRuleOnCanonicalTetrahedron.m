function [GP,GW] = getGaussRuleOnCanonicalTetrahedron(n)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Gauss Points and Weights for the exact integration of a
% polynomial of degree n = p + q, where p is the polynomial degree in the
% u-direction and q is the polynomial degree in the v-direction. The 
% function works up to n = 8. The Gaussian integration is exact over 
% the canonical element:
%   
%         (0,1)
%           |\
%           | \
%           |  \
%           |   \
%           |    \
%           |     \
%           |      \
%           ---------
%         (0,0)     (1,0)
%
% Source : Quadrature Formulas in Two Dimensions, 
%          Math 5172 - Finite Element Method, Section 001, Spring 2010
%
%   Input :
%       n : The polynomial degree of the integrand
%
%  Output :
%      GP : The set of Gauss Points with their three linear dependent 
%           coordinates GPi = [1-lambda1 - lambda2 lambda1 lambda2]
%      GW : The respective set of Gauss weights 
%
% Function layout :
%
% 1. Assign the linear independent Gauss Point coordinates and their respective weights
%
% 2. Add the linear dependent coordinate to the Gauss Point array
%
%% Function main body

%% 1. Assign the linear independent Gauss Point coordinates and their respective weights
if n==1
    GPLI = [0.250000000000000 0.250000000000000 0.250000000000000];
    GW = 1.00000000000000;
elseif n==2
    warning('Gauss rule for 2 Gauss points needs to be checked')
    GPLI = [0.138196601125011 0.138196601125011 0.138196601125011
            0.138196601125011 0.585410196624969 0.138196601125011
            0.585410196624969 0.138196601125011 0.138196601125011
            0.138196601125011 0.138196601125011 0.585410196624969];
    GW = [0.250000000000000
          0.250000000000000
          0.250000000000000
          0.250000000000000];
elseif n==3
    GPLI  =  [0.250000000000000 0.250000000000000 0.250000000000000
              0.166666666666667 0.166666666666667 0.166666666666667
              0.166666666666667 0.500000000000000 0.166666666666667
              0.500000000000000 0.166666666666667 0.166666666666667
              0.166666666666667 0.166666666666667 0.450000000000000];
    GW = [-0.800000000000000
           0.450000000000000
           0.450000000000000
           0.450000000000000
           0.450000000000000];
elseif n==4
    error('gauss rule not implemented!')
elseif n==5
    error('gauss rule not implemented!')
elseif n==6
    error('gauss rule not implemented!')
elseif n==7
    error('gauss rule not implemented!')
elseif n==8
    error('gauss rule not implemented!')
else
    error('Gauss Points and weights for the integration of polynomial integrands of order higher than 8 over triangular domain have not been yet implemented');
end

%% 2. Add the linear dependent coordinate to the Gauss Point array
noGP = length(GPLI(:,1));
GP = zeros(noGP,3);
GP(:,2:4) = GPLI;
GP(:,1) = ones(noGP,1) - GP(:,2) - GP(:,3) - GP(:,4);

end
