function outputPostProc = computePostProc(FComplete, parameters, postProc, domainName, noDimensions)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Marko Leskovar
%                  Andreas Apostolatos
%
%% Function documentation
%
% function that handle post processing defined in git
%
%               Input :
%           FComplete : The complete force vector
%            postProc : Post-processing properties
%          domainName : name of the domain of interest within out problem
%          parameters : Flow parameters
%        noDimensions : dimensionality of the problem 
%                       2 for 2D problems, 3 for 3D problems
%   
%              Output :
%      outputPostProc : variable output defined by function handle when we
%                       setup the problem in GiD
%      
%% Function main body

% prealocate a logical array
indexArray = false(1, length(postProc.nameDomain));

% find the function handle and domains index
for k = 1:length(postProc.nameDomain)
    indexArray(k) = (domainName == postProc.nameDomain(k));
end

% get the nodes and function handle
nodesDomain = postProc.nodesDomain{indexArray};
functionHandle = postProc.computePostProc{indexArray};

% compute forces
if(strcmp(functionHandle,'computeForces'))

    % get the name of the function and compute output
    outputFunction = str2func(functionHandle);
    outputPostProc = outputFunction(FComplete, parameters, nodesDomain, noDimensions);
end

% other functions can be added here in the same way as "compute forces"

end




