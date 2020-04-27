function [Xir, CPr] = knotRefineUniformlyBSplineCurve(n, p, Xi, CP, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Control Points and the knot vector after equidistant knot 
% insertion into the given B-Spline curve
%
%  Input :
%      n : Number of knots to be inserted uniformly
%      p : The polynomial degree of the curve
%     Xi : The knot vector of the B-Spline curve
%     CP : The Control Points of the B-Spline curve
% outMsg : Whether or not to output message on refinement progress
%          'outputEnabled' : enables output information
%
% Output :
%    Xir : The knot vector after the refinement
%    CPr : The set of Control points after the refinement
%
% Function layout :
%
% 1. Knot refine uniformly
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________________\n');
    fprintf('##############################################################\n');
    fprintf('Uniform Knot insertion for a B-Spline curve has been initiated \n\n');
    fprintf('Number of knots before knot insertion nxi = %d \n',length(Xi));
    fprintf('Number of knots after knot insertion nxi = %d \n',length(Xi) + n);
    fprintf('______________________________________________________________\n\n');
    tic;
end

%% 1. Knot refine uniformly

% Initialize the knot vector containing the additional knots 
Ru = zeros(n-1);

% Loop over all the knots to be inserted
for j = 1:n-1  
    Ru(j) = j/n*(Xi(length(Xi))-Xi(1)) ;  
end

% Apply knot insertion onto the curve
[Xir,CPr] = knotRefineBSplineCurve(p,Xi,CP,Ru,'');

%% 2. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Degree elevation took %d seconds \n\n', computationalTime);
    fprintf('_____________________Knot Insertion Ended_____________________\n');
    fprintf('##############################################################\n\n\n');
end

end

