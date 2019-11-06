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
function [Xir,CPr,pr] = degreeElevateBSplineCurve(p,Xi,CP,tp,outMsg)
%% Function documentation
%
% elevates B-Spline curve from polynomial degree p  by tp. This function
% handles B-Spline curves in 3D. Source reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag, Berlin 1995
% p. 72.
%
%  Input :
%      p : The polynomial degree of the curve
%     Xi : The knot vector of the NURBS curve
%     CP : The Control Points of the NURBS curve in 3D
%     tp : polynomial degree by which to elevate the curve
% outMsg : Whether or not to output message on refinement progress
%          'outputEnabled' : enables output information
%
% Output :
%    Xir : The knot vector after refinement
%    CPr : The set of Control points after refinement
%     pr : pr = p + tp, the degree of the curve after refinement
%
% Function layout :
%
% 0. Read input
%
% 1. Compute the projective Control Points Pw
%
% 2. Degree elevate to polynomial degree p + tp
%
% 3. Retransform the projective Control Points from Qw to CPr
%
% 4. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('________________________________________________________\n');
    fprintf('########################################################\n');
    fprintf('Degree elevation for a B-Spline curve has been initiated \n\n');
    fprintf('Polynomial degree before elevation p = %d\n',p);
    fprintf('Polynomial degree after elevation p = %d\n',p + tp);
    fprintf('________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read input

% Compute the number of Control Points
nxi = length(CP(:,1));

% Number of knots of the unrefined geometry
mxi = nxi + p + 1;

% Polynomial degree of the curve after refinement
pr = p + tp;

% Initialize arrays
bezalfs = zeros(pr+1,p+1);
bpts = zeros(p+1,length(CP(1,:)));
ebpts = zeros(pr+1,length(CP(1,:)));
Nextbpts = zeros(p-1,length(CP(1,:)));
Pw = zeros(nxi,length(CP(1,:)));

%% 1. Compute the projective Control Points Pw
for i = 1:nxi
    Pw(i,1:3) = CP(i,1:3)*CP(i,4);
    Pw(i,4)   = CP(i,4);
end

%% 2. Degree elevate to polynomial degree p + tp
if tp>0
    % preallocate Ur
    ku = 1;
    
    for i=1:mxi-1
        if Xi(i)~=Xi(i+1)  
            ku = ku + 1;   
        end
    end 
    Xir = zeros(1,mxi+tp*ku);

    % compute Bezier degree elevation coefficients
    bezalfs(1,1) = 1.0;
    bezalfs(pr+1,p+1) = 1.0;

    % bezalfs are symmetric to row pr/2
    for i = 1:fix(pr/2)
        inv = 1/getBinomialCoefficients(pr,i);
        mpi = min(p,i);
        
        % alphas for degree elevation
        for j = max(0,i-tp):mpi
            bezalfs(i+1,j+1) = inv*getBinomialCoefficients(p,j)*getBinomialCoefficients(tp,i-j);
            bezalfs(pr-i+1,p-j+1) = bezalfs(i+1,j+1);    % symmetry !
        end
    end

    mh = pr;
    kind = pr+1;
    r = -1;
    a = p + 1;
    b = p + 2;
    cind = 1;
    ua = Xi(a);
    Qw(1,:) = Pw(1,:);

    % left end of Ur
    for i = 1:pr+1
        Xir(i) = ua;
    end

    % initialize first Bezier seg
    for i = 1:p+1
        bpts(i,:) = Pw(i,:);
    end

    % big loop through knot vector
    while b<mxi
	i = b;
        
        % make b the rightmost ocurrence of ub
        while b<mxi && Xi(b)==Xi(b+1)
            b = b + 1;
        end
        
        % multiplicity of ub
        mult = b - i + 1;
        mh = mh + mult + tp;
        ub = Xi(b);
        
        % r from last segment
        oldr = r;
        
        % ub to be inserted r times
        r = p - mult;
  
        if oldr>0    
            lbz = fix((oldr+2)/2);
        else
            lbz = 1;
        end
  
        if r>0
            rbz = pr-fix((r+1)/2);
        else
            rbz = pr;
        end
        
        % insert knots to get Bezier segment
        alphas = zeros(1,p-mult);
        if r>0
            numer = ub - ua;
        
            % alpha for knot insertion
            for k = p:-1:mult+1
                alphas(k-mult) = numer/(Xi(a+k+1)-ua);
            end
        
            % r times knot insertion
            for j = 1:r
                save = r - j + 1;
                s = mult + j;
            
                % new CP due to knot insertion
                for k = p+1:-1:s+1
                    bpts(k,:) = alphas(k-s)*bpts(k,:)+(1-alphas(k-s))*bpts(k-1,:);
                end
                Nextbpts(save,:) = bpts(p+1,:);
            end
        end
        % end of inserting knots
  
        % degree elevate Bezier
        for i = lbz:pr
            % only points lbz,..,pr are used below
            ebpts(i+1,1:4) = 0;
            mpi = min(p,i);
      
            % new CP due to degree elevation
            for j = max(0,i-tp):mpi
                ebpts(i+1,:) = ebpts(i+1,:) + bezalfs(i+1,j+1)*bpts(j+1,:);
            end
        end
        % end of degree elevating Bezier
  
        % knot removal ua oldr times
        if oldr>1
            first = kind-2;
            last = kind;
            den = ub-ua;
            bet = (ub-Xir(kind))/den;
      
            for tr = 1:oldr-1
                i = first;
                j = last;
                kj = j-kind+1;
            
                while j-i>tr
                    % loop and compute the new CP for one removal step
                    if i<cind
                        alpha = (ub-Xir(i+1))/(ua-Xir(i+1));
                        Qw(i+1,:) = (alpha*Qw(i+1,:)+(1.0-alpha)*Qw(i,:));
                    end
                    if j>=lbz
                        if ((j-tr)<=(kind-pr+oldr))
                            gam = (ub-Xir(j-tr+1))/den;
                            ebpts(kj+1,:) = gam*ebpts(kj+1,:) + (1.0-gam)*ebpts(kj+2,:);
                        else
                            ebpts(kj+1,:) = bet*ebpts(kj+1,:) + (1.0-bet)*ebpts(kj+2,:);
                        end
                    end
                    i = i + 1;
                    j = j - 1;
                    kj = kj - 1;
                end
                
                first = first - 1;
                last = last + 1;
            end
        end
        % end of removing knot, u=U(a)
     
        % load the knot ua
        if a~=p+1
            for i = 0:(pr-oldr-1)
                Xir(kind+1) = ua;
                kind = kind+1;
            end
        end
  
        % load CPs into Qw
        for j = lbz:rbz
            Qw(cind+1,:) =  ebpts(j+1,:);
            cind = cind +1;
        end
  
        % setup for the next pass through loop
        if b<mxi
            for j = 0:r-1
                bpts(j+1,:) = Nextbpts(j+1,:);
            end
            for j = r:p
                bpts(j+1,:) = Pw(b-p+j,:);
            end
            a = b;
            b = b+1;
            ua = ub;
            % end knots
        else
            for i = 0:pr
                Xir(kind+i+1) = ub;
            end
        end
    end
    % end of big loop through knot vector
elseif (tp==0)
    Xir = Xi;  
    pr = p;
    Qw = Pw;
end

%% 3. Retransform the projective Control Points from Qw to CPr

% Initialize output array
CPr = zeros(length(Qw(:,1)),length(CP(1,:)));

% Loop over all the projective Control Points
for i = 1:length(Qw(:,1))  
    CPr(i,1:3) = Qw(i,1:3)/Qw(i,4);
    CPr(i,4)   = Qw(i,4); 
end

%% 4. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Degree elevation took %d seconds \n\n',computationalTime);
    fprintf('_________________Degree Elevation Ended_________________\n');
    fprintf('########################################################\n\n\n');
end

end
