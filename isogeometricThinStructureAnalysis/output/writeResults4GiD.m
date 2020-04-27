function writeResults4GiD ...
    (analysis, BSplinePatches, u, uDot, uDDot, propNLinearAnalysis, ...
    propTransientAnalysis, propPostproc, caseName, pathToOutput, title, ...
    noTimeStep)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Writes out the results of an isogeometric thin-walled structural analysis
% into a GiD result file which can be read by GiD given that it is located
% in the same folder (with extension .gid) written out from the function 
% writeOutMultipatchBSplineSurface4GiD.m
%
%                  Input : 
%              analysis : Dummy variable for this function
%        BSplinePatches : The B-Spline patch consisting of its polynomial
%                         orders, knot vectors and Control Point 
%                         coordinates and weights
%                     u : Discrete solution to the primary variable of the
%                         problem
%                  uDot : Dummy variable for this function
%                 uDDot : Dummy variable for this function
%   propNLinearAnalysis : Dummy variable for this function
% propTransientAnalysis : On the properties of the structural dynamics
%                         analysis:
%                                 .scheme : Time integration scheme
%                                 .TStart : Starting time of the simulation
%                                   .TEnd : End time of the simulation
%                            .noTimeSteps : Number of time steps
%                                     .dt : Time step size
%          propPostproc : On the postprocessing        
%                              .resultant : Array of strings containing the 
%                                           name of each resultant to be 
%                                           written out
%                       .computeResultant : Array of strings representing
%                                           the function handle for the
%                                           computation of the desirable
%                                           resultant
%               caseName : The name of the case after which the file will 
%                          be named
%           pathToOutput : Absolute path to the folder containing the 
%                          results file
%                  title : Dummy variable for this function
%             noTimeStep : Number of the current time step
%
% Function layout :
%
% 0. Read input
%
% 1. Make directory to write out the results of the analysis
%
% 2. Write the preamble
%
% 3. Loop over all the resultants
% ->
%    3i. Get the resultants name
%
%   3ii. Loop over the patches in the multipatch geometry
%   ->
%        3ii.1. Write the patch number
%
%        3ii.2. Get the vector of DOFs for the current patch
%
%        3ii.3. Get the vector of DOFs for the current patch containing only displacements
%
%        3ii.4. Loop over the Control Points of the given patch
%   <-
%
%  3iii. Indicate the ending of the values
% <-
%
%% Function main body

%% 0. Read input

% Number of DOFs per Control Point
if strcmp(analysis.type,'isogeometricKirchhoffLoveShellAnalysis') || ...
        strcmp(analysis.type,'isogeometricMembraneAnalysis')
    noDOFsPerCP = 3;
else
    error('Output for the selected analysis type %s is not yet implemented',analysis.type);
end

% Number of patches
noPatches = length(BSplinePatches);

%% 1. Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName,'.gid'),'dir');
if ~isExistent
    error('Directory %s is not existent, function writeOutMultipatchBSplineSurface4GiD has to be called in advanced',strcat(pathToOutput,caseName,'.gid'));
end
isUpdated = isExistent && noTimeStep ~= 0;
if isUpdated
    fileHandle = fopen(strcat(pathToOutput,caseName,'.gid','/',caseName,'.post.res'),'a+');
else
    fileHandle = fopen(strcat(pathToOutput,caseName,'.gid','/',caseName,'.post.res'),'w');
end

%% 2. Write the preamble
if noTimeStep == 0
    fprintf(fileHandle,'GiD Post Results File 1.2 \n\n');
end

%% 3. Loop over all the resultants
for iRes = 1:length(propPostproc.resultant)
    %% 3i. Get the resultants name
    nameResultant = propPostproc.resultant{iRes};
    
    if strcmp(propPostproc.resultant,'displacement')
        fprintf(fileHandle,['Result "' nameResultant '" "time_step" %d Vector OnNurbsSurface\n'],noTimeStep);
    end
    fprintf(fileHandle,'Values\n');

    %% 3ii. Loop over the patches in the multipatch geometry
    if strcmp(nameResultant,'displacement')
        for iPatches = 1:noPatches
            %% 3ii.1. Write the patch number
            fprintf(fileHandle,'%d\n',iPatches);
            
            %% 3ii.2. Get the vector of DOFs for the current patch
            uPatch = u(BSplinePatches{iPatches}.EFTPatches);
            
            %% 3ii.3. Get the vector of DOFs for the current patch containing only displacements
            nxi = length(BSplinePatches{iPatches}.CP(:,1,1));
            neta = length(BSplinePatches{iPatches}.CP(1,:,1));
            noDOFsDisp = noDOFsPerCP*nxi*neta;
            uPatchDisp = uPatch(1:noDOFsDisp);
            
            %% 3ii.4. Loop over the Control Points of the given patch
            for iCPsPatch = 1:BSplinePatches{iPatches}.noCPs
                if noDOFsPerCP == 3
                    for iComp = 1:3
                        if uPatchDisp(3*iCPsPatch - ( 3 - iComp)) ~= 0
                            fprintf(fileHandle,'%2.15f',uPatchDisp(3*iCPsPatch - (3 - iComp)));
                        else
                            fprintf(fileHandle,'%2.15f',uPatchDisp(3*iCPsPatch - (3 - iComp)));
                        end
                        if iComp == 1 || iComp == 2
                            fprintf(fileHandle,' ');
                        else
                            fprintf(fileHandle,'\n');
                        end
                    end
                end
            end
        end
    else
        error('The implementation for the computation of resultant %s is yet missing',nameResultant);
    end
    
    %% 3iii. Indicate the ending of the values
    fprintf(fileHandle,'End Values\n\n');
end

end
