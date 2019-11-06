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
function plot_knotsForBSplineSurfaceOnPArametricSpace...
    (Xi,Eta,xiGrid,etaGrid)
%% Function documentation
%
% Draws the element edges for the NURBS surface, i.e the knots on the
% parametric space
%
%         Input :
%        Xi,Eta : Knot vectors in xi,eta-direction
%        xiGrid : Points to use in xi-direction
%       etaGrid : Points to use in eta-direction
%
%     Output : 
%              Graphics
%
% Function layout :
%
% 0. Read input
%
% 1. Plot the edges in xi-direction
%
% 2. Plot the edges in eta-direction
%
% 3. Plot all the element edges
%
%% Function main body

%% 0. Read input

% Number of knots in u,v-direction
mxi = length(Xi);
meta = length(Eta);

% Make the knot vectors unique
XiUnique = unique(Xi);
EtaUnique = unique(Eta);

% Color of the edge
color_edge = 'black';
% color_edge = 'red';
% color_edge = 'none';

% Assign a tolerance value
eps = 1e-9;

% Initialize counter
l = 1;

% Compute step size in xi-direction
dxi = (Xi(mxi)-Xi(1))/xiGrid;

% Compute step size in eta-direction
deta = (Eta(meta)-Eta(1))/etaGrid; 

% Initialize plotting array
P = zeros(xiGrid,etaGrid,2);

%% 1. Plot the edges in xi-direction
for j2 = 1:length(EtaUnique)
    eta = EtaUnique(j2);
    xi = XiUnique(1);
    k = 1;
    while xi <= XiUnique(end) + eps 
        P(k,l,1:2) = [xi eta];
        k = k + 1;
        xi = xi + dxi;
    end
    l = l + 1;
end

%% 2. Plot the edges in eta-direction
for i2 = 1:length(XiUnique)
    xi = XiUnique(i2);
    eta = EtaUnique(1);
    k = 1;
    while eta <= EtaUnique(end) + eps
        P(k,l,1:2) = [xi eta];
        k = k + 1;
        eta = eta + deta;
    end
    l = l + 1;
end

%% 3. Plot all the element edges
plot(P(:,:,1),P(:,:,2),'Color',color_edge,'LineWidth',.01);
axis equal;
grid on;
xlabel('x','FontSize',18);
ylabel('y','FontSize',18);
zlabel('z','FontSize',18);
    
end