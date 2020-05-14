function  plot_EdgesTriangularMesh2D(mesh)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Plots the edges of the given triangular mesh
%
%   Input :
%    mesh : The nodes, the edges and the coonectivity of the triabgular
%           mesh
%
%  Output : 
%           Graphics
% 
%% Function main body

% Plot the edges of the mesh
cla,patch('vertices',mesh.nodes(:,2:end),'faces',mesh.elements(:,2:end),'edgecol','k','facecol','none');
% cla,patch('vertices',mesh.nodes,'faces',mesh.elements,'edgecol','k','facecol',[217 218 219]/255);

end
