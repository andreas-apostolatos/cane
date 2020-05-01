function X = normr(Y)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Normalizes the rows of the given array
%
%	Input :
%		Y : Given array with not normalized rows
%
%  Output :
%	 	X : The array with the normalized rows
%
%% Function main body
X = Y./repmat(sqrt(sum(Y.*Y, 2)),1, size(Y, 2));

end
