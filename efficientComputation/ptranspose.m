function MTranspose = ptranspose(M)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Transposes the pagewise matrix array by switching the second and third
% dimensions (as the first dimension loops over all pages).
%
%  Input :
%      M : The 3-d array of matrices to transpose. The size( M ) 
%          should be [nMat, m, n], where nMat is the number of matrices to
%          multiply and m, n are the dimensions of each matrix
%
% Output :
% result : The transpose of M
% 
%% Function main body
MTranspose = permute(M,[1,3,2]);

end
