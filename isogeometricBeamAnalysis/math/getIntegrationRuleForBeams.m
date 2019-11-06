function int = getIntegrationRuleForBeams(int,analysis)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Returns the Gauss Point coordinates and weights for the integration of
% the stiffness of the isogeometric Bernoulli or Timoshenko beam. There are
% two types available: Automatic integration, manual choice of the Gauss 
% Point number and the efficient rule.
%
%    Input :
%      int : Number of Gauss points, and other parameters related to the
%            efficient quadrature if chosen
% analysis : Beam analysis type :
%                   'Bernoulli' : isogeometric Bernoulli beam analysis
%                  'Timoshenko' : isogeometric Timoshenko beam analysis
%
%   Output :
%      int : The integration structure updated
%
%% Function main body

% Contains the number of all Gauss points that have been used
int.nGP = 0;

% On the integration using the efficient quadrature rule
if strcmp(int.type,'efficientRule')
    % Find the interrelement regularity
    q = 0;
    isTheSame = true;
    while isTheSame
        q = q + 1;
        isTheSame = (U(p+2)==U(p+2+q));
    end
    % The interrelement regularity of the Splines before differentiating
    q = p-q;
    % Integrate all parts with the same number of quadrature points
    if int.partial == 0
        % The space is S_{q-2}^{2*(p-1)}
        % Due to the presence of 2nd derivatives the continuity decreases
        q = q-2;
        % Maximal polynomial degree 
        m = 2*(p-2);
        [int.xI,int.wI,int.xL,int.wL,int.xR,int.wR] = efficientQuadratureUniform(m,q);
   % Integrate membrane and bending part with appropriate number of
   % quadrature points
    else
        % For the membrane part the space of the basis is S_{q-1}^{2*(p-1)}
        mMembrane = 2*(p-1);
        qMembrane = q-1;
        [int.xI_mem,int.wI_mem,int.xL_mem,int.wL_mem,int.xR_mem,int.wR_mem] =...
            efficientQuadratureUniform(mMembrane,qMembrane);
        % For the bending part the space of the basis is S_{q-2}^{2*(p-2)}
        mBending = 2*(p-1);
        qBending = q-2;
        [int.xI_ben,int.wI_ben,int.xL_ben,int.wL_ben,int.xR_ben,int.wR_ben] =...
            efficientQuadratureUniform(mBending,qBending);
        if strcmp(analysis.type,'Timoshenko')
            % For the shear part the space of the basis is S_{q-1}^{2*p}
            mShear = 2*p;
            qShear = q-1;
            [int.xI_she,int.wI_she,int.xL_she,int.wL_she,int.xR_she,int.wR_she] =...
                quadrature_database_uniform(mShear,qShear);
        end
    end
end

