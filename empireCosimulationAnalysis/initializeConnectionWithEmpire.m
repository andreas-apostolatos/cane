function initializeConnectionWithEmpire(strMatlabXml,tab,outMsg)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
% Initializes connection with empire through MPI.
%
%        Input : 
% strMatlabXml : String containing the name of the xml file needed for
%                establishing connection with Empire
%          tab : Tabulation when printing message in the command window
%       outMsg : Enables message outputting in the command window if chosen
%                as 'outputEnabled'
%
%       Output : 
%                Establishing connection with Empire
%
% Function layout :
%
% 0. Check input
%
% 1. Establish connection with Empire
%
%% Function main body

%% 0. Check input
if strcmp(outMsg,'outputEnabled')
    fprintf(strcat(tab,'Connecting to Empire using file %s.xml\n\n'),strMatlabXml);
end
    
if isempty(strMatlabXml)
    error('No xml file for the Matlab code is provided');
end

%% 1. Establish connection with Empire

% Get the string containing the path to the directory where the
% environment variables from EMPIRE are
empireBaseDir = getenv('EMPIRE_CORE_BASE_DIR');
if isempty(empireBaseDir)
    error('Empire base directory could not be found')
end

% Create a temporary directory containing the path to the Matlab interface
% functions from Empire
tmpDir = '/interface/matlab/binICC';

% Concatenate the strings containing the directory names
mexFuncPath = strcat(empireBaseDir, tmpDir);

% Add the directory to the existing path of included directories
path(path, mexFuncPath);

% Connect to Empire using the provided xml file
EMPIRE_API_Connect(strcat(strMatlabXml,'.xml'));

end
