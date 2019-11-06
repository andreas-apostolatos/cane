function [ determinant ] = p3Ddet( M )
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Computes determinants of 3x3 matrices pagewise (by rule of Sarrus).
%
%  Input :
%      M : The 3-d array of matrices. The size( M ) should be [n, 3, 3].
% 
% Output :
%    det : The vector of determinants. The size( det ) is [n, 1, 1].
% 
%% Function main body

if size(M, 2) ~= 3 || size(M, 3) ~= 3 || length(size(M)) ~= 3
    error('p3Ddet: error in matrix dimensions!\n');
end

determinant = M(:, 1, 1) .* M(:, 2, 2) .* M(:, 3, 3) ...
            + M(:, 1, 2) .* M(:, 2, 3) .* M(:, 3, 1) ...
            + M(:, 1, 3) .* M(:, 2, 1) .* M(:, 3, 2) ...
            - M(:, 3, 1) .* M(:, 2, 2) .* M(:, 1, 3) ...
            - M(:, 3, 2) .* M(:, 2, 3) .* M(:, 1, 1) ...
            - M(:, 3, 3) .* M(:, 2, 1) .* M(:, 1, 2);

end
