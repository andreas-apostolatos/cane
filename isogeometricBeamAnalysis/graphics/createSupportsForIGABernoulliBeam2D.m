function [xs, ys, zs] = createSupportsForIGABernoulliBeam2D ...
    (CP, rb)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns coordinates for the triangles at the support locations for the
% isogeometric Bernoulli beam in 2D
%
%   Input : 
%      CP : Control Point coordinates
%      rd : Array containing information on the supports
%
%  Output :
%      xs : x-coordinates of the support triangle vertices
%      ys : y-coordinates of the support triangle vertices
%      zs : z-coordinates of the support triangle vertices
%
%% Function main body

% Total number of Control Points
nu = length(CP(:,1));

% scaling factors for the support triangles
maximum = max(max(max(max(CP))));
minimum = min(min(min(min(CP))));

% Average the factor with respect to the maximum and minimum values 
factor = (maximum-minimum)/5;

% Initialize the output arrays
xs = zeros(length(rb),4);
ys = zeros(length(rb),4);
zs = zeros(length(rb),4);

for k = 1:length(rb)
    % Get the corresponding Control Point number p and indices CP(i,j)
    h=rb(k)/2;
    p=ceil(h);
    j=ceil(p/nu);
    i=p-(j-1)*nu;
    
    %(rb is odd -> horizontal support)
    if (p~=h)   
        xs(k,1)=CP(i,1);
        xs(k,2)=CP(i,1)-0.1732*factor;
        xs(k,3)=CP(i,1)-0.1732*factor;
        xs(k,4)=xs(k,1);
        ys(k,1)=CP(i,2);
        ys(k,2)=CP(i,2)+0.1*factor;
        ys(k,3)=CP(i,2)-0.1*factor;
        ys(k,4)=ys(k,1);
        zs(k,1:4)=CP(i,3);
    %(rb is even -> vertical support)
    else        
        xs(k,1)=CP(i,1);
        xs(k,2)=CP(i,1)-0.1*factor;
        xs(k,3)=CP(i,1)+0.1*factor;
        xs(k,4)=xs(k,1);
        ys(k,1)=CP(i,2);
        ys(k,2)=CP(i,2)-0.1732*factor;
        ys(k,3)=CP(i,2)-0.1732*factor;
        ys(k,4)=ys(k,1);
        zs(k,1:4)=CP(i,3);
    end
end

end
