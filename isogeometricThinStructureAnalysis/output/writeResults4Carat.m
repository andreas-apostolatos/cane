function writeResults4Carat ...
    (analysis, BSplinePatches, u, uDot, uDDot, propNLinearAnalysis, ...
    propPostproc, caseName, pathToOutput, title, noTimeStep, GiDResultFile)
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
% into a Carat++ result file which can be read by GiD given that it is
% located in the same folder as the file writen out from the function
% writeOutMultipatchBSplineSurface4Carat.m
%
%              Input : 
%            analysis : Dummy variable for this function
%        BSplinePatch : The B-Spline patch consisting of its polynomial
%                       orders, knot vectors and Control Point coordinates
%                       and weights
%                   u : Discrete solution to the primary variable of the
%                       problem
%                uDot ; Dummy variable for this function
%               uDDot : Dummy variable for this function
% propNLinearAnalysis : Dummy variable for this function
%        propPostproc :        .resultant : Array of strings containing the 
%                                           name of each resultant to be 
%                                           written out
%                       .computeResultant : Array of strings representing
%                                           the function handle for the
%                                           computation of the desirable
%                                           resultant
%            caseName : The name of the case after which the file will be
%                       named
%        pathToOutput : Absolute path to the folder containing the results
%                       file
%               title : Dummy variable for this function
%          noTimeStep : Number of the current time step
%       GiDResultFile : Dummy variable for this function
%
% Function layout :
%
% 1. Make directory to write out the results of the analysis
%
% 2. Write the displacement field for each Control Point
%
%% Function main body

%% 0. Read input

% Initialize counter
counter = 1;

% Check if the geometry is a multipatch NURBS geometry
noPatches = length(BSplinePatches);
if noPatches ~= 1
    warning('This function is can handle only single patch NURBS geometries');
end

%% 1. Make directory to write out the results of the analysis
isExistent = exist(strcat(pathToOutput,caseName),'dir');
if ~isExistent
    mkdir(strcat(pathToOutput,caseName));
end
fileHandle = fopen(strcat(pathToOutput,caseName,'/',caseName,'_GiD.post.res'),'a+');

%% 2. Write the displacement field for each Control Point
%(fid,gtype,ngaus,job)
% write header requered for the post
%
%GiDResultFile = fopen('out.post.res','w');
if noTimeStep == 1
    fprintf(fileHandle,'Rhino Post Results File 1.0 \n');
end
fprintf(fileHandle,'Result "Displacement" "Load Case" %d Vector OnNodes\n',noTimeStep);
fprintf(fileHandle,'Values\n');
for iPatches = 1:noPatches
    for iCPs = 1:BSplinePatches{iPatches}.noCPs
        fprintf(fileHandle,'  %d  %.15f  %.15f  %.15f\n',counter,u(3*counter-2),u(3*counter-1),u(3*counter));
        counter = counter + 1;
    end
end
fprintf(fileHandle,'End Values\n');
%fclose(GiDResultFile);

end
