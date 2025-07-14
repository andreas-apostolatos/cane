function [xs, ys, zs] = createSupports5D(CP, rb)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation 
%
% Returns the coordinates of the triangles needed for drawing the supports
% Modified version for 5 DOFs per control point.
%
%       Input : 
%          CP : The set of Control Points and coordinates
%          rb : Global numbering of the suported DoFs
%
%
%      Output :
%    xs,ys,zs : The coordinates of the triangles
%
%% Function main body

nu = length(CP(:,1,1));

% scaling factor
up = max(max(max(CP)));
lo = min(min(min(CP)));
fac = (up-lo)/50;

% Initialize the output arrays
xs = zeros(length(rb),4);
ys = zeros(length(rb),4);
zs = zeros(length(rb),4);

for l = 1:length(rb)
    % get the corresponding Control Point number p and indices CP(i,j) for 5 DOFs
    p=ceil(rb(l)/5);  % Control point number assuming 5 DOFs per CP
    j=ceil(p/nu);     % Get j index 
    i=p-(j-1)*nu;     % Get i index
    dir=rb(l)-((j-1)*nu+i-1)*5;  % Get direction (1,2,3,4,5)
    
    % Only plot supports for the first 3 DOFs (translational)
    if (dir==1)              %(x-support)
        xs(l,1)=CP(i,j,1);
        xs(l,2)=CP(i,j,1)-1.732*fac;
        xs(l,3)=CP(i,j,1)-1.732*fac;
        xs(l,4)=xs(l,1);
        ys(l,1)=CP(i,j,2);
        ys(l,2)=CP(i,j,2)+fac;
        ys(l,3)=CP(i,j,2)-fac;
        ys(l,4)=ys(l,1);
        zs(l,1:4)=CP(i,j,3);
    elseif (dir==2)          %(y-support)
        xs(l,1)=CP(i,j,1);
        xs(l,2)=CP(i,j,1)-fac;
        xs(l,3)=CP(i,j,1)+fac;
        xs(l,4)=xs(l,1);
        ys(l,1)=CP(i,j,2);
        ys(l,2)=CP(i,j,2)-1.732*fac;
        ys(l,3)=CP(i,j,2)-1.732*fac;
        ys(l,4)=ys(l,1);
        zs(l,1:4)=CP(i,j,3);
    elseif (dir==3)          %(z-support)
        xs(l,1)=CP(i,j,1);
        xs(l,2)=CP(i,j,1)-fac;
        xs(l,3)=CP(i,j,1)+fac;
        xs(l,4)=xs(l,1);
        ys(l,1:4)=CP(i,j,2);
        zs(l,1)=CP(i,j,3);
        zs(l,2)=CP(i,j,3)-1.732*fac;
        zs(l,3)=CP(i,j,3)-1.732*fac;
        zs(l,4)=zs(l,1);
    elseif (dir==4 || dir==5)  %(rotational DOFs - no visual support for now)
        % Set to zero or NaN to indicate no visual representation
        xs(l,1:4) = NaN;
        ys(l,1:4) = NaN;
        zs(l,1:4) = NaN;
    end
end

end