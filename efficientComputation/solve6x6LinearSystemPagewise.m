function x = solve6x6LinearSystemPagewise(A, b)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the solution to a 6x6 linear equation system of the special form
%
% |A11 A12| . |x11 x12 x13| = |b11 b12 b13|
% |A21 A22|   |x21 x22 x23|   |b21 b22 b23|
%
% when matrices A and b are pagewise formed matrices.
%
%   Input :
%       A : Matrix of dimensions (noPages,2,2)
%       b : Matrix of dimensions (noPages,2,3)
%   
%  Output :
%       x : Solution vector of dimensions (noPages,2,3)
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the pagewise solution vector
%
%% Function main body

%% 0. Read input

if length(size(A)) ~= 3
    error('Input matrix A has to be a pagewise matrix');
end
if length(size(b)) ~= 3
    error('Input matrix b has to be a pagewise matrix');
end
if size(A,2) ~= 2 || size(A,3) ~= 2
    error('Input matix A has to have dimensions (noPages,2,2)');
end
if size(b,2) ~= 2 || size(b,3) ~= 3
    error('Input matix A has to have dimensions (noPages,2,3)');
end

%% 1. Compute the pagewise solution vector
%
% x11 = (A(:,2,2).*b(:,1,1) - A(:,1,2).*b(:,2,1))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1));
% x12 = (A(:,2,2).*b(:,1,2) - A(:,1,2).*b(:,2,2))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1));
% x13 = (A(:,2,2).*b(:,1,3) - A(:,1,2).*b(:,2,3))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1));
% x21 = (A(:,1,1).*b(:,2,1) - A(:,2,1).*b(:,1,1))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1));
% x22 = (A(:,1,1).*b(:,2,2) - A(:,2,1).*b(:,1,2))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1));
% x23 = (A(:,1,1).*b(:,2,3) - A(:,2,1).*b(:,1,3))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1));
%
x = cat(3,cat(2,(A(:,2,2).*b(:,1,1) - A(:,1,2).*b(:,2,1))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1)),...
    (A(:,1,1).*b(:,2,1) - A(:,2,1).*b(:,1,1))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1))),...
    cat(2,(A(:,2,2).*b(:,1,2) - A(:,1,2).*b(:,2,2))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1)),...
    (A(:,1,1).*b(:,2,2) - A(:,2,1).*b(:,1,2))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1))),...
    cat(2,(A(:,2,2).*b(:,1,3) - A(:,1,2).*b(:,2,3))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1)),...
    (A(:,1,1).*b(:,2,3) - A(:,2,1).*b(:,1,3))./(A(:,1,1).*A(:,2,2) - A(:,1,2).*A(:,2,1))));

end
