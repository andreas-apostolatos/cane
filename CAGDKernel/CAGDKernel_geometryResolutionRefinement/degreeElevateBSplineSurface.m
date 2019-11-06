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
function [Xir,Etar,CPr,pr,qr] = degreeElevateBSplineSurface...
    (p,q,Xi,Eta,CP,tp,tq,outMsg)
%% Function documentation
%
% Elevates polynomial degrees p and q by tp and tq to the given NURBS
% curve. Source reference:
%
% Les Piegl and Wayne Tiller, The NURBS Book. Springer-Verlag, Berlin 1995
% p. 72.
%
%    Input :
%      p,q : the polynomial degrees
%   Xi,Eta : the knot vectors in xi-, eta-directions
%       CP : control point coordinates and weights
%       tp : integer by which to order elevate in xi-direction
%       tq : integer by which to order elevate in eta-direction
%   outMsg : Whether or not to output message on refinement progress
%            'outputEnabled' : enables output information
%
%   Output :
% Xir,Etar : the new knot vectors
%      CPr : the new set of control points
%    pr,qr : the update of the polynomial degrees
%
% Function layout :
%
% 0. Read Input
%
% 1.Compute the projective control points Pw 
%
% 2. Degree elevate in xi-direction
%
% 3. Degree elevate in eta-direction
%
% 4. Retransform from 4th dimensional to the 3D space
%
% 5. Appendix
%
%% Function main body
if strcmp(outMsg,'outputEnabled')
    fprintf('__________________________________________________________\n');
    fprintf('##########################################################\n');
    fprintf('Degree elevation for a B-Spline srface has been initiated \n\n');
    fprintf('Polynomial degree in xi-direction before elevation p = %d\n',p);
    fprintf('Polynomial degree in xi-direction after elevation p = %d\n',p + tp);
	fprintf('Polynomial degree in eta-direction before elevation q = %d\n',q);
    fprintf('Polynomial degree in eta-direction after elevation q = %d\n',q + tq);
    fprintf('__________________________________________________________\n\n');

    % start measuring computational time
    tic;
end

%% 0. Read Input

% Number of Control Points in u- and in v- direction
nxi = length(CP(:,1,1));
neta = length(CP(1,:,1));

% Number of knots of the unrefined geometry
mxi = nxi + p + 1;
meta = neta + q + 1;

% Polynomial degree of the curve after refinement
pr = p + tp;
qr = q + tq;

% Initializations
Pw = zeros(nxi,neta,4);

%% 1.Compute the projective control points Pw 
for i = 1:nxi
    for j = 1:neta
        Pw(i,j,1:3) = CP(i,j,1:3)*CP(i,j,4);
        Pw(i,j,4) = CP(i,j,4);
    end
end

%% 2. Perform degree elevation in xi-direction

% Initializations
bezalfs = zeros(pr+1,p+1);
bpts = zeros(p+1,length(CP(1,:,1)),length(CP(1,1,:)));
ebpts = zeros(pr+1,length(CP(1,:,1)),length(CP(1,1,:)));
Nextbpts = zeros(p-1,length(CP(1,:,1)),length(CP(1,1,:)));

% If the given polynomial degree is a positive integer
if tp>0
    % preallocate Ur
    ku = 1;
    for i=1:mxi-1  
        if Xi(i)~=Xi(i+1)
            ku = ku+1;
        end
    end
    
    Xir = zeros(1,mxi+tp*ku);

    % compute Bezier degree elevation coefficients
    bezalfs(1,1) = 1.0;
    bezalfs(pr+1,p+1) = 1.0;

    for i = 1:fix(pr/2)   % bezalfs are symmetric to row pr/2
        inv = 1/getBinomialCoefficients(pr,i);
        mpi = min(p,i);
        % alfas for degree elevation
        for j = max(0,i-tp):mpi
            bezalfs(i+1,j+1) = inv*getBinomialCoefficients(p,j)*getBinomialCoefficients(tp,i-j);
            % Apply symmetry
            bezalfs(pr-i+1,p-j+1) = bezalfs(i+1,j+1);
        end
    end

    mh = pr;
    kind = pr+1;
    r=-1;
    a=p+1;
    b=p+2;
    cind=1;
    ua = Xi(a);
    Qw(1,:,:) = Pw(1,:,:);

    % left end of Ur
    for i = 1:pr+1
        Xir(i) = ua;
    end

    % initialize first Bezier seg
    for i = 1:p+1
        bpts(i,:,:) = Pw(i,:,:);
    end

    % big loop through knot vector
    while b<mxi
        i = b;
        while b<mxi && Xi(b)==Xi(b+1)
            % make b the rightmost ocurrence of ub
            b = b + 1;      
        end
        
        % multiplicity of ub
        mult = b-i+1;   
        mh = mh + mult + tp;
        ub = Xi(b);
        
        % r from last segment
        oldr = r;
        
        % ub to be inserted r times
        r = p-mult;
  
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
        if (r>0)
            numer = ub - ua;
            for k = p:-1:mult+1
                % alpha for knot insertion
                alphas(k-mult) = numer/(Xi(a+k)-ua);   
            end
            % r times knot insertion
            for j = 1:r   
                save = r - j + 1;
                s = mult+j;
                % new CP due to knot insertion
                for k = p+1:-1:s+1    
                    bpts(k,:,:) = alphas(k-s)*bpts(k,:,:)+(1-alphas(k-s))*bpts(k-1,:,:);
                end
                Nextbpts(save,:,:) = bpts(p+1,:,:);
            end
        end
        % end of inserting knots
  
        % degree elevate Bezier
        for i = lbz:pr
            % only points lbz,..,pr are used below
            ebpts(i+1,1:neta,1:4) = 0;
            mpi = min(p,i);
            % new CP due to degree elevation
            for j = max(0,i-tp):mpi   
                ebpts(i+1,:,:) = ebpts(i+1,:,:) + bezalfs(i+1,j+1)*bpts(j+1,:,:);
            end
        end
        % end of degree elevating Bezier
  
        % knot removal ua oldr times
        if oldr>1
            first = kind-2;
            last = kind;
            den = ub-ua;
            beta = (ub-Xir(kind))/den;
            for tr = 1:oldr-1
                i = first;
                j = last;
                kj = j-kind+1;

                while j-i>tr
                    % loop and compute the new CP for one removal step
                    if i<cind
                        alpha = (ub-Xir(i+1))/(ua-Xir(i+1));
                        Qw(i+1,:,:) = (alpha*Qw(i+1,:,:)+(1.0-alpha)*Qw(i,:,:));
                    end
                    if j>=lbz
                        if j-tr <= kind-pr+oldr
                            gamma = (ub-Xir(j-tr+1))/den;
                            ebpts(kj+1,:,:) = gamma*ebpts(kj+1,:,:) + (1.0-gamma)*ebpts(kj+2,:,:);
                        else
                            ebpts(kj+1,:,:) = beta*ebpts(kj+1,:,:) + (1.0-beta)*ebpts(kj+2,:,:);
                        end
                    end
                    i = i+1;
                    j = j-1;
                    kj = kj-1;
                end
                first = first-1;
                last = last+1;
            end
        end
        % end of removing knot, u=U(a)
     
        if (a~=p+1)   
            % load the knot ua
            for i = 0:(pr-oldr-1)
                Xir(kind+1) = ua;
                kind = kind+1;
            end
        end
  
        for j = lbz:rbz   
            % load CPs into Qw
            Qw(cind+1,:,:) =  ebpts(j+1,:,:);
            cind = cind +1;
        end
  
        if b<mxi   % setup for the next pass through loop
            for j = 0:r-1
                bpts(j+1,:,:) = Nextbpts(j+1,:,:);
            end
            for j = r:p
                bpts(j+1,:,:) = Pw(b-p+j,:,:);
            end
            a = b;
            b = b+1;
            ua = ub;
      else   % end knots
            for i = 0:pr
                Xir(kind+i+1) = ub;
            end
        end
    end
    % end of big loop through knot vector
elseif tp==0
	Xir = Xi;  
    pr = p;
    Qw = Pw;
end

%% 3. Perform degree elevation in eta-direction

% Initializations
bezalfs = zeros(qr+1,q+1);
bpts = zeros(q+1,length(CP(1,:,1)));
ebpts = zeros(qr+1,length(CP(1,:,1)));
Nextbpts = zeros(q-1,length(CP(1,:,1)));

% If the given polynomial degree is a positive integer
if tq>0
    clear ('Pw');
    Pw(:,:,:) = Qw(:,:,:);
    nxi = length(Pw(:,1,1));
    clear ('Qw');  
    clear ('bpts');
    clear ('Nextbpts');
    clear ('ebpts');
    clear ('bezalfs');
    
    % preallocate Vr
    kv = 1;
    for i=1:meta-1  
        if Eta(i)~=Eta(i+1)
            kv = kv+1;
        end
    end
    Etar = zeros(1,meta+tq*kv);

    % compute Bezier degree elevation coefficients
    bezalfs(1,1) = 1.0;
    bezalfs(qr+1,q+1) = 1.0;

    % bezalfs are symmetric to row qr/2
    for i = 1:fix(qr/2)
        inv = 1/getBinomialCoefficients(qr,i);
        mpi = min(q,i);
        
        % alfas for degree elevation
        for j = max(0,i-tq):mpi   
            bezalfs(i+1,j+1) = inv*getBinomialCoefficients(q,j)*getBinomialCoefficients(tq,i-j);
            % Apply the symmetry
            bezalfs(qr-i+1,q-j+1) = bezalfs(i+1,j+1);
        end
    end

    mh = qr;
    kind = qr+1;
    r = -1;
    a = q + 1;
    b = q + 2;
    cind = 1;
    ua = Eta(a);
    Qw(:,1,:) = Pw(:,1,:);

    % left end of Vr
    for i = 1:qr+1
        Etar(i) = ua;
    end

    % initialize first Bezier seg
    for i = 1:q+1
        bpts(:,i,:) = Pw(:,i,:);
    end

    % big loop through knot vector
    while b<meta
        i = b;
        while b<meta && Eta(b)==Eta(b+1)
            % make b the rightmost ocurrence of ub
            b = b + 1;      
        end
        
        % multiplicity of ub
        mult = b - i + 1;   
        mh = mh + mult + tq;
        ub = Eta(b);
        
        % r from last segment
        oldr = r;   
        
        % ub to be inserted r times
        r = q - mult;     
  
        if oldr>0
            lbz = fix((oldr+2)/2);
        else
            lbz = 1;
        end
  
        if r>0
            rbz = qr-fix((r+1)/2);
        else
            rbz = qr;
        end
  
        % insert knots to get Bezier segment
        if r>0
            numer = ub - ua;
            for k = q:-1:mult+1
                % alfa for knot insertion
                alphas(k-mult) = numer/(Eta(a+k)-ua);   
            end
            
            % r times knot insertion
            for j = 1:r   
                save = r-j+1;
                s = mult + j;
            
                % new CP due to knot insertion
                for k = q+1:-1:s+1    
                    bpts(:,k,:) = alphas(k-s)*bpts(:,k,:)+(1-alphas(k-s))*bpts(:,k-1,:);
                end
                Nextbpts(:,save,:) = bpts(:,q+1,:);
            end
        end
        
        % end of inserting knots
  
        % degree elevate Bezier
        for i = lbz:qr
            % only points lbz,..,qr are used below
            ebpts(1:nxi,i+1,1:4) = 0;
            mpi = min(q,i);
            for j = max(0,i-tq):mpi   % new CP due to degree elevation
                ebpts(:,i+1,:) = ebpts(:,i+1,:) + bezalfs(i+1,j+1)*bpts(:,j+1,:);
            end
        end
        % end of degree elevating Bezier
  
        % knot removal ua oldr times
        if oldr>1
            first = kind - 2;
            last = kind;
            den = ub-ua;
            beta = (ub-Etar(kind))/den;
            for tr = 1:oldr-1
                i = first;
                j = last;
                kj = j-kind+1;
                
                % loop and compute the new CP for one removal step
                while j-i>tr
                    if i<cind
                        alpha = (ub-Etar(i+1))/(ua-Etar(i+1));
                        Qw(:,i+1,:) = (alpha*Qw(:,i+1,:)+(1.0-alpha)*Qw(:,i,:));
                    end
                    if j>=lbz
                        if ((j-tr)<=(kind-qr+oldr))
                            gamma = (ub-Etar(j-tr+1))/den;
                            ebpts(:,kj+1,:) = gamma*ebpts(:,kj+1,:) + (1.0-gamma)*ebpts(:,kj+2,:);
                        else
                            ebpts(:,kj+1,:) = beta*ebpts(:,kj+1,:) + (1.0-beta)*ebpts(:,kj+2,:);
                        end
                    end
                    i = i + 1;
                    j = j - 1;
                    kj = kj - 1;
                end
                first = first-1;
                last = last+1;
            end
        end
        % end of removing knot, u=V(a)
     
        % load the knot ua
        if a~=q+1 
            for i = 0:qr-oldr-1
                Etar(kind+1) = ua;
                kind = kind + 1;
            end
        end
  
        % load CPs into Qw
        for j = lbz:rbz   
            Qw(:,cind+1,:) =  ebpts(:,j+1,:);
            cind = cind +1;
        end
  
        % setup for the next pass through loop
        if b<meta
            for j = 0:r-1
                bpts(:,j+1,:) = Nextbpts(:,j+1,:);
            end
            for j = r:q
                bpts(:,j+1,:) = Pw(:,b-q+j,:);
            end
            a = b;
            b = b+1;
            ua = ub;
        
        % end knots
        else   
            for i = 0:qr
                Etar(kind+i+1) = ub;
            end
        end
    end
  % end of big loop through knot vector

elseif tq==0
    Etar = Eta;  
    qr = q;
end

%% 4. Retransform from 4th dimensional to the 3D space

% Initialize output array
CPr = zeros(length(Qw(:,1,1)),length(Qw(1,:,1)),length(CP(1,1,:)));

% Loop over all the projective Control Points
for i = 1:length(Qw(:,1,1))
    for j = 1:length(Qw(1,:,1))
        CPr(i,j,1:3) = Qw(i,j,1:3)/Qw(i,j,4);
        CPr(i,j,4)   = Qw(i,j,4);
    end
end

%% 5. Appendix
if strcmp(outMsg,'outputEnabled')
    % Save computational time
    computationalTime = toc;

    fprintf('Degree elevation took %.2d seconds \n\n',computationalTime);
    fprintf('__________________Degree Elevation Ended__________________\n');
    fprintf('##########################################################\n\n\n');
end

end
