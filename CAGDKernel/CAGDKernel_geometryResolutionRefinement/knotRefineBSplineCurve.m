%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%   Technische Universität München                                        %
%   Lehrstuhl für Statik, Prof. Dr.-Ing. Kai-Uwe Bletzinger               %
%   _______________________________________________________               %
%   _______________________________________________________               %
%                                                                         %
%                                                                         %
%   Authors                                                               %
%   _______________________________________________________________       %
%                                                                         %
%   Dipl.-Math. Andreas Apostolatos    (andreas.apostolatos@tum.de)       %
%   Dr.-Ing. Roland Wüchner            (wuechner@tum.de)                  %
%   Prof. Dr.-Ing. Kai-Uwe Bletzinger  (kub@tum.de)                       %
%   _______________________________________________________________       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Xir,CPr] = knotRefineBSplineCurve(p,Xi,CP,Rxi,outMsg)
%% Function documentation
%
% Perform knot insertion for a B-Spline curve in 3D space. Source 
% reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag: Berlin 1995
% p. 72.
%
%   Input :
%       p : The polynomial degree of the curve
%      Xi : The knot vector of the NURBS curve
%      CP : The Control Points of the NURBS curve in 3D
%     Rxi : The set of knots to be inserted into the knot vector Xi
%  outMsg : Whether or not to output message on refinement progress
%           'outputEnabled' : enables output information
%
%  Output :
%     Xir : The knot vector after the refinement
%     CPr : The set of Control points after the refinement
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the projective control points Pw
%
% 2. Knot refine the curve
%
% 3. Transform the projective Control Points Qw to CPr
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('______________________________________________________\n');
    fprintf('######################################################\n');
    fprintf('Knot insertion for a B-Spline curve has been initiated \n\n');
    fprintf('Number of knots before knot insertion nxi = %d\n',length(Xi));
    fprintf('Number of knots after knot insertion nxi = %d\n',length(Xi)+length(Rxi));
    fprintf('______________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Compute the number of Control Points of the unrefined geometry
nxi = length(CP(:,1));

% Compute the number of knots of the unrefined geometry
mxi = nxi + p + 1;

% Compute the number of knots of the refined knot vector
mxir = mxi + length(Rxi);

% Compute the number of Control Points of the refined geometry
nxir = mxir - p - 1;

% Initialize arrays
Pw = zeros(nxi,length(CP(1,:)));
Xir = zeros(1,mxir);
Qw = zeros(nxir,length(CP(1,:)));
CPr = zeros(nxir,length(CP(1,:)));

%% 1. Compute the projective control points Pw
for i = 1:nxi
    Pw(i,1:3) = CP(i,1:3)*CP(i,4);
    Pw(i,4)   = CP(i,4);
end

%% 2. Knot refine the curve
if isempty(Rxi)
    Xir = Xi;
    Qw = Pw;
else
    r = length(Rxi); 
    a = findKnotSpan(Rxi(1),Xi,nxi);
    b = findKnotSpan(Rxi(r),Xi,nxi)+1;
 
    for j = 1:a-p 
        Qw(j,:) = Pw(j,:);   
    end
    
    for j = b-1:nxi
        Qw(j+r,:) = Pw(j,:);   
    end
 
    for j = 1:a      
        Xir(j) = Xi(j); 
    end
    
    for j = b+p:mxi   
        Xir(j+r) = Xi(j);  
    end
    
    i = b + p - 1;   
    k = i + r;
    
    for  j = r:-1:1
        while Rxi(j)<=Xi(i) && i>a
            Qw(k-p-1,:) = Pw(i-p-1,:);
            Xir(k) = Xi(i);
            k = k-1;   i = i-1;
        end
        
        Qw(k-p-1,:) = Qw(k-p,:);
    
        for l = 1:p
            ind = k-p+l;
            alpha = (Rxi(j)-Xir(k+l)) / (Xi(i-p+l)-Xir(k+l));
        
            Qw(ind-1,:) = alpha*Qw(ind-1,:) + (1-alpha)*Qw(ind,:);
     
        end
        
        Xir(k) = Rxi(j);
        k = k - 1;
    end
end

%% 3. Transform the projective Control Points Qw to CPr
for i = 1:length(Qw(:,1))
    CPr(i,1:3) = Qw(i,1:3)/Qw(i,4);
    CPr(i,4)   = Qw(i,4);
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('\t Knot insertion took %d seconds \n\n',computationalTime);
    fprintf('_________________Knot Insertion Ended_________________\n');
    fprintf('######################################################\n\n\n');
end

end