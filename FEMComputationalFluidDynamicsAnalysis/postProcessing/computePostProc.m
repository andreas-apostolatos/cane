function postProc = computePostProc(FComplete, analysis, parameters, postProc)
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
%            analysis : Analysis type and number of dimensions
%          parameters : Flow parameters
%            postProc : Post-processing properties
%   
%                     Output :
%     postProc.valuePostProc : variable output defined by function handle 
%                              when we setup the problem in GiD
%                       
%% Function main body

% loop through the domains
for k = 1:length(postProc.nameDomain)
    
    % get the nodes and function handle
    nodesDomain = postProc.nodesDomain{k};
    functionHandle = postProc.computePostProc{k};
    
    % get the name of the function and compute output
    outputFunction = str2func(functionHandle);
    postProc.valuePostProc{k} = outputFunction(FComplete, analysis,     ...
                                               parameters, nodesDomain);
                                    
end




