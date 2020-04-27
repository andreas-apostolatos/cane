function [ result ] = phorzcat(varargin)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Concatenates the 3D matrix arrays horizontally (means the third dimension): 
% in 2-d equal to [input_1, input_2, ... , input_n]
%
%    Input :
% varargin : A variable number if 3-d input arrays that should be 
%            concatenated horizontally. The size in the first two dimensions
%            must match.
%
% Output :
% result : The concatenated 3-d aray
%
%% Function main body
result = cat(3, varargin{:});

end
