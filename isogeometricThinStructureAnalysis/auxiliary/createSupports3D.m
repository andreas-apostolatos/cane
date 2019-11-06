function [xs,ys,zs] = createSupports3D(CP,rb)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
% Returns the coordinates of the triangles needed for drawing the supports
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
% fac=0.1;

% Initialize the output arrays
xs = zeros(length(rb),4);
ys = zeros(length(rb),4);
zs = zeros(length(rb),4);

for l = 1:length(rb)
    % get the corresponding Control Point number p and indices CP(i,j)
    p=ceil(rb(l)/3);
    j=ceil(p/nu);
    i=p-(j-1)*nu;
    dir=rb(l)-((j-1)*nu+i-1)*3;
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
    end
end

end
