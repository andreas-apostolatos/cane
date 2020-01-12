function c = mergesorted(a,b)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
%% Function documentation
%
% Merge and sorts two given vectors a and b into c
%
% Input :
%  a, b : Two vectors to be merged and sorted
%
% Output :
%      c : The merged and sorted vector
%
%% Function main body

lena = length(a);
lenb = length(b);
c = zeros(1,lena+lenb);

% index to move along vector 'a'
inda = 1;

% index to move along vector 'b'
indb = 1;

% index to move along vector 'c'
indc = 1;

while ((inda <= lena) && (indb <= lenb))
 if a(inda) < b(indb)
    c(indc) = a(inda);
    inda = inda + 1;
 else
    c(indc) = b(indb);
    indb = indb + 1;
 end
 indc = indc + 1;
end

% Copy any remaining elements of the 'a' into 'c'
while (inda <= lena)
  c(indc) = a(inda);
  indc = indc + 1;
  inda = inda + 1;
end

% Copy any remaining elements of the 'b' into 'c'
while (indb <= lenb)
  c(indc) = b(indb);
  indc = indc + 1;
  indb = indb + 1;
end

end