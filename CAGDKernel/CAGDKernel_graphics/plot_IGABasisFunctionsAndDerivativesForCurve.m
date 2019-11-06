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
function index = plot_IGABasisFunctionsAndDerivativesForCurve(p,Xi,CP,numEval,isNURBS,graph,outMsg)
%% Function documentation
%
% Plots the NURBS basis functions and their derivatives up to the desired
% derivative for the parametrization of a curve.
%
%    Input :
%        p : The polynomial degree of the curve
%       Xi : The knot vector of the NURBS curve
%       CP : The Control Points of the NURBS curve
%  numEval : How many points to use for the plot
%  isNURBS : Flag on whether the geometrical basis is a NURBS or a B-Spline
%    graph : Structure containing all information for the graphs
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
%   Output : 
%   index  : The new index of the graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the polynomial degrees
%
%    1i. Loop over all the parametric locations
%
%        1ii.1. Find the knot span
%
%        1ii.2. Compute the NURBS basis functions and their derivatives
%
%        1ii.3. Save the functions and their derivatives
%
% 2. Plot all the graphs
%
% 3. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('___________________________________________________________\n');
    fprintf('###########################################################\n');
    if isNURBS
        fprintf('Plotting the NURBS basis functions and their derivatives \n\n');
    else
        fprintf('Plotting the B-Spline basis functions and their derivatives \n\n');
    end
    if graph.plotBasisFunctionsAndDerivs == 0
        fprintf('No output plot has been chosen to be displayed\n');
    elseif graph.plotBasisFunctionsAndDerivs == 1
        fprintf('Only the basis functions themselves were selected to be\n');
        fprintf('displayed\n');
    else
        fprintf('The basis functions and their derivatives up to order %d \n',graph.plotBasisFunctionsAndDerivs-1);
        fprintf('were chosen to be displayed \n');
    end
    fprintf('___________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Number of derivatives to be displayed
nDeriv = graph.plotBasisFunctionsAndDerivs;

% Increment
step = (Xi(length(Xi))-Xi(1))/numEval;

% values at which the basis functions are calculated
xi = Xi(1):step:Xi(length(Xi));

% Compute the number of basis functions
nBF = length(Xi)-p-1;

% Initialize matrix containing the values of the base functions at xi
dRi = zeros(nBF,length(xi),graph.plotBasisFunctionsAndDerivs);

%% 1. Loop over all the polynomial degrees
for i = 1:nDeriv
    %% 1i. Loop over all the parametric locations
    for j=1:length(xi)
        %% 1ii.1. Find the knot span
        knotSpan = findKnotSpan(xi(j),Xi,nBF);
        
        %% 1ii.2. Compute the B-Spline/NURBS basis functions and their derivatives
        dR = computeIGABasisFunctionsAndDerivativesForCurve(knotSpan,p,xi(j),Xi,CP,isNURBS,nDeriv);
        
        %% 1ii.3. Save the functions and their derivatives
        dRi(findKnotSpan(xi(j),Xi,nBF)-p:findKnotSpan(xi(j),Xi,nBF),j,i) = dR(:,i);
    end
end

%% 2. Plot all the graphs
for i=1:nDeriv
    figure (graph.index)
    plot(xi(:),dRi(:,:,i));
    graph.index = graph.index + 1;
    xlabel('\xi \in \Xi');
    if i==1
        if isNURBS
            title('The NURBS basis functions'); 
        else
            title('The B-Spline basis functions'); 
        end
        ylabel('N_i(\xi)');
    elseif i==2
        if isNURBS
            title(sprintf('The first derivatives of the NURBS basis functions'));
        else
            title(sprintf('The first derivatives of the B-Spline basis functions'));
        end
        ylabel('$\displaystyle\frac{dN_i}{d\xi}$','interpreter','latex');
    elseif i==3
        if isNURBS
            title(sprintf('The second derivatives of the NURBS basis functions'));
        else
            title(sprintf('The second derivatives of the B-Spline basis functions'));
        end
        ylabel('$\displaystyle\frac{d^2N_i}{d^2\xi}$','interpreter','latex');
    elseif i==4
        if isNURBS
            title(sprintf('The third derivatives of the NURBS basis functions'));
        else
            title(sprintf('The third derivatives of the B-Spline basis functions'));
        end
        ylabel('$\displaystyle\frac{d^3N_i}{d^3\xi}$','interpreter','latex');
    else
        if isNURBS
            title(sprintf('The %d-th derivatives of the NURBS basis functions',i-1));
        else
            title(sprintf('The %d-th derivatives of the B-Spline basis functions',i-1));
        end
        ylabel('$\displaystyle\frac{d^{i}N_i}{d^{i}\xi}$','interpreter','latex');
    end
end

% Save the graphics index
index = graph.index;

%% 3. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Plotting the basis functions took %d seconds \n\n',computationalTime);
    fprintf('_______________Plotting Basis Functions Ended______________\n');
    fprintf('###########################################################\n\n\n');
end

end