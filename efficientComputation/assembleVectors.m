function [varargout] = assembleVectors(EFT, numDOFs, numDOFsEl, varargin)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Given an arbitrary number of vectors which are storing the individual
% element level vectors pagewisely and the element freedom tables also
% pagewise represented, the function returns the same number of assembled
% vectors globally.
%
%       Input :
%         EFT : The element freedom tables in a pagewise representation
%     numDOFs : The number of DOFs globally
%   numDOFsEl : The number of DOFs in the element level
%    varargin : An arbitrary number of element level vectors stored
%               pagewisely
%
%      Output :
%   varargout : The assembled global vectors
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the input vectors
% ->
%    1i. Get the dimensions of the input vector
%   
%   1ii. Check the input vector for inconsistency in its dimensions#
%
%  1iii. Compute the size of the vectors and initialize output vector
%
%   1iv. Reshape indices to a vector for the element freedom tables
%
%    1v. Reshape the values of the element vectors into a vector
%
%    1vi. Assemble the contributions from each element vector to the global one
% <-
%
%% Function main body

%% 0. Read input

% number of elements is the second dimension of the element freedom table
noElmnts = size(EFT,2);

% initialize output cell array
varargout = cell(nargin-3,1);

%% 1. Loop over all the input vectors
for matID = 1:nargin - 3
    %% 1i. Get the dimensions of the input vector
    sizeOfInput = size(varargin{matID});
    if length(sizeOfInput) == 2
        sizeOfInput(3) = 1;
    end
    
    %% 1ii. Check the input vector for inconsistency in its dimensions
    if ~any(sizeOfInput(2:3) == 1) || ~prod(sizeOfInput(2:3)) == numDOFsEl
        error('Assembly of vectors has failed due to inconsistent dimensions');
    end
    
    %% 1iii. Compute the size of the vectors and initialize output vector
    sizeOfOutput = ( sizeOfInput(2:3) ~= 1 ) * (numDOFs - 1) + 1;
    sizeOfReshapedVectors = ( sizeOfInput(2:3) ~= 1 ) * (noElmnts * numDOFsEl - 1) + 1;
    varargout{matID} = zeros(sizeOfOutput(1), sizeOfOutput(2));
    
    %% 1iv. Reshape indices to a vector for the element freedom tables
    i = reshape(transpose(EFT), sizeOfReshapedVectors(1), sizeOfReshapedVectors(2));
    
    %% 1v. Reshape the values of the element vectors into a vector
    s = reshape(varargin{matID}, sizeOfReshapedVectors(1), sizeOfReshapedVectors(2));
    
    %% 1vi. Assemble the contributions from each element vector to the global one
    for j = 1:length(i)
        index = i(j);
        varargout{matID}(index) = varargout{matID}(index) + s(j);
    end
end

end
