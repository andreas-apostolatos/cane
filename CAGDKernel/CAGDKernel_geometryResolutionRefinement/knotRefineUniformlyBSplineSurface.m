function [Xir, Etar, CPr] = knotRefineUniformlyBSplineSurface ...
    (p, Xi, q, Eta, CP, noXi, noEta, outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Computes the CP's and U after equidistant knot insertion into the given
% curve
%
%      Input :
%   nxi,neta : Number of knots to be inserted equidistantly in xi-,eta 
%              direction
%        p,q : The polynomial degrees of the surface
% noXi,noEta : Number of Control Points along xi- and eta- parametric
%              directions
%     Xi,Eta : The knot vectors of the NURBS surface in xi,eta-direction
%         CP : The Control Points of the NURBS surface
%     outMsg : Whether or not to output message on refinement progress
%              'outputEnabled' : enables output information
%
%     Output :
%        Xir : The knot vector in xi-direction after the refinement
%       Etar : The knot vector in eta-direction after the refinement
%        CPr : The set of Control points after the refinement
%
% Function layout :
%
% 0. Read input
%
% 1. Create the vectors with the knots to be added in both parametric directions
%
% 2. Appendix
%
%% Function main body
if strcmp(outMsg, 'outputEnabled')
    fprintf('________________________________________________________________\n');
    fprintf('################################################################\n');
    fprintf('Uniform Knot insertion for a B-Spline curve has been initiated \n\n');
    fprintf('Number of knots before knot insertion in xi-direction nxi = %d\n', length(Xi));
    fprintf('Number of knots after knot insertion in xi-direction nxi = %d\n', length(Xi) + noXi);
    fprintf('Number of knots before knot insertion in eta-direction neta = %d\n', length(Eta));
    fprintf('Number of knots after knot insertion in eta-direction neta = %d\n', length(Eta) + noEta);
    fprintf('________________________________________________________________\n\n');
    tic;
end

%% 0. Read input

% Initialize the vectors containing the additional knots 
Rxi = [];
Reta = [];

%% 1. Create the vectors with the knots to be added in both parametric directions

% Loop and add knots into the vector in xi-direction
counterXi = 1;
for i = 1:noXi - 1
    xiKnot = i/noXi*(Xi(length(Xi)) - Xi(1));
    index = find(xiKnot == Xi, 1);
    isKnot2Add = true;
    if ~isempty(index)
        noMultipleKnots = length(index);
        if noMultipleKnots == p
            if strcmp(outMsg,'outputEnabled')
                warning('Knot vector Xi has already a %d multiplicity at knot %d creating discontinuity, no knot will be added', p, xiKnot);
            end
            isKnot2Add = false;
        elseif noMultipleKnots == p - 1
            if strcmp(outMsg,'outputEnabled')
                warning('Knot vector Xi has already a %d multiplicity at knot %d creating interpolation, no knot will be added', p, xiKnot);
            end
            isKnot2Add = false;
        end
    end
    if isKnot2Add
        Rxi(counterXi) = xiKnot;
        counterXi = counterXi + 1;
    end
end

% Loop and add knots into the vector in eta-direction
counterEta = 1;
for i = 1:noEta - 1
    etaKnot = i/noEta*(Eta(length(Eta)) - Eta(1));
    index = find(etaKnot == Eta, 1);
    isKnot2Add = true;
    if ~isempty(index)
        noMultipleKnots = length(index);
        if noMultipleKnots == q
            if strcmp(outMsg,'outputEnabled')
                warning('Knot vector Eta has already a %d multiplicity at knot %d creating discontinuity, no knot will be added', p, etaKnot);
            end
            isKnot2Add = false;
        elseif noMultipleKnots == q - 1
            if strcmp(outMsg, 'outputEnabled')
                warning('Knot vector Eta has already a %d multiplicity at knot %d creating interpolation, no knot will be added', p, etaKnot);
            end
            isKnot2Add = false;
        end
    end
    if isKnot2Add
        Reta(counterEta) = etaKnot;
        counterEta = counterEta + 1;
    end
end

% Apply knot insertion onto the surface
[Xir, Etar, CPr] = knotRefineBSplineSurface ...
    (p, Xi, q, Eta, CP, Rxi, Reta, '');

%% 2. Appendix
if strcmp(outMsg,'outputEnabled')
    computationalTime = toc;
    fprintf('Knot insertion took %.2d seconds \n\n',computationalTime);
    fprintf('______________________Knot Insertion Ended______________________\n');
    fprintf('################################################################\n\n\n');
end

end
