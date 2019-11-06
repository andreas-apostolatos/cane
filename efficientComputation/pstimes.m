function [ mat ] = pstimes( mat, scalar )
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Performs pagewise computation of matrix times scalar.
%
%  Input :
%    mat : The 3-d array of first operand matrices. The size( mat1 ) 
%          should be [nMat, m, n], where nMat is the number of matrices to
%          multiply and m, n are the dimensions of each matrix. There is
%          no restriction on m and n
% scalar : The array of scalars. The size( scalar ) should be [nMat, 1].
%
%          (the number of entries in the first dimension has to be equal 
%          for both input parameters)
%
% Output :
%    mat : The pagewise product of each matrix in mat with the corresponding 
%          scalar in scalar. The size of the result is the same as the size
%          of the input matrix
% 
%% Function main body

scalarSize = size(scalar);

if length(scalarSize) > 2 || prod(scalarSize) == 0        
    error('page-wise scalar multiplication: error in scalar input dimensions!\n')
end

if( scalarSize(2) > 1 )
    if( scalarSize(1) > 1 )
         error('page-wise multiplication: error in scalar input dimensions!\n')
    end
    scalar = scalar';
end

mat = bsxfun(@times,mat,scalar);

end
