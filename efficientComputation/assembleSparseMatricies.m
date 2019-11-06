function [varargout] = assembleSparseMatricies(EFT,noDOFs,noDOFsEl,varargin)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Given an arbitrary number of matrices which are storing the individual
% element level matrices pagewisely and the element freedom tables also
% pagewise represented, the function returns the same number of assembled
% matrices globally and in a sparse form.
%
%            Input :
%              EFT : The element freedom tables of all elements
%           noDOFs : The total number of degrees of freedom 
%         noDOFsEl : The number of dofs per element
%         varargin : The input element matrices
%
%           Output :
%        varargout : The assembled global matrices
%
% Function layout :
%
% 0. Read input
%
% 1. loop over all input matrices
% ->
%    1i. Check input matrix
%
%   1ii. Get the dimensions of the input matrix
%
%  1iii. Check the input matrix for inconsistency in its dimensions
%
%   1iv. Compute first indices of non-zero matrix entries
%
%    1v. Compute second indices of non-zero matrix entries
%
%   1iv. Compute values of non-zero matrix entries
%
%  1vii. Create the sparse matrix from the i, j and s arrays
% <-
%
%% Function main body

%% 0. Read input

% number of elements is the second dimension of the element freedom table
noElem = size(EFT,2);

% transform the element freedom table in the right shape
EFT = permute(EFT,[3,1,2]);

% initialize output cell array
varargout = cell(nargin - 3,1);

%% 1. loop over all input matrices
for matID = 1:nargin-3
    %% 1i. Check input matrix
    if ischar(varargin{matID})
        varargout{matID} = 'undefined';
        continue;
    end
    
    %% 1ii. Get the dimensions of the input matrix
    sizeOfInput = size(varargin{matID});
    if length(sizeOfInput) == 2
        sizeOfInput(3) = 1;
    end
    
    %% 1iii. Check the input matrix for inconsistency in its dimensions
    if sizeOfInput(2) ~= sizeOfInput(3) || sizeOfInput(2) ~= noDOFsEl
        error('Assembly of sparse matrices has failed due to inconsistent dimensions')
    end
    
    %% 1iv. Compute first indices of non-zero matrix entries
    i = permute(reshape(reshape(permute(EFT(ones(1,noDOFsEl),:,:),[2,1,3]), ...
        [noDOFsEl,1,noDOFsEl*noElem]),[1,1,noDOFsEl^2*noElem]),[3,2,1]);

    %% 1v. Compute second indices of non-zero matrix entries
    j = permute(reshape(reshape(EFT(ones(1,noDOFsEl),:,:),[noDOFsEl,1, ...
        noDOFsEl*noElem]),[1,1,noDOFsEl^2*noElem]),[3,2,1]);

    %% 1vi. Compute values of non-zero matrix entries
    s = permute(reshape(reshape(permute(varargin{matID},[2,3,1]),...
        [noDOFsEl,1,noDOFsEl*noElem]),[1,1,noDOFsEl^2*noElem]),[3,2,1]);

    %% 1vii. Create the sparse matrix from the i, j and s arrays
    varargout{matID} = sparse(i,j,s,noDOFs,noDOFs);
end
    
end
