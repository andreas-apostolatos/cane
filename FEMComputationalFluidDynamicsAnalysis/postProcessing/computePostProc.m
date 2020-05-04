function postProc = computePostProc ...
    (FComplete, propAnalysis, propParameters, postProc)
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
% Function that handle post processing defined in GiD and adds a new field
% to the struct postProc
%
%               Input :
%           FComplete : The complete force vector
%        propAnalysis : Structure on general analysis properties,
%                           .type : Analysis type
%      propParameters : Structure on the fluid parameters,
%                            .nue : Dynamic viscosity
%            postProc : Structure on postprocessing properties,
%                           .nameDomain : Name of domain of interest
%                          .nodesDomain : ID of nodes on domain of interest
%                      .computePostProc : Function handle on the
%                                         computation of the desirable 
%                                         quantity of interest
%   
%              Output :
%            postProc : Updated structure on postprocessing properties,
%                        .valuePostProc : variable output defined by 
%                                         function handle when we setup the 
%                                         problem in GiD
%
% Function layout :
%
% 0. Read input
%
% 1. Loop over all the domains of interest
% ->
%    1i. Get the nodes and function handle
%
%   1ii. Get the name of the function handle
%
%  1iii. Compute the desirable output
% <-
%                       
%% Function main body

%% 0. Read input
if ~isstruct(postProc)
    error('postProc must be a structure defining the postprocessing properties');
else
    if ~isfield(postProc, 'nameDomain')
        error('postProc must define variable nameDomain');
    end
    if ~isfield(postProc, 'nodesDomain')
        error('postProc must define variable nodesDomain');
    end
    if ~isfield(postProc, 'computePostProc')
        error('postProc must define variable computePostProc');
    end
end

%% 1. Loop over all the domains of interest
for iDomains = 1:length(postProc.nameDomain)
    %% 1i. Get the nodes and function handle
    nodesDomain = postProc.nodesDomain{iDomains};
    functionHandle = postProc.computePostProc{iDomains};
    
    %% 1ii. Get the name of the function handle
    outputFunction = str2func(functionHandle);
    
    %% 1iii. Compute the desirable output
    postProc.valuePostProc{iDomains} = outputFunction ...
        (propAnalysis, nodesDomain, FComplete, propParameters);                                
end

end
