%% SETUP
% -------------------------------------------------------------------------
% caseName          : GiD .dat file to be used for the calculation
% samplingCall      : desired input sampling method
%                   - 'randomUniform'
%                   - 'randomNormal'
%                   - 'latinHypercube'
%                   - 'quasiMonteCarloHalton'
%
% outputFunction    : function that processes the output of the solver and
%                     computes the quantity of interest
% nSample           : number of samples to run
% -------------------------------------------------------------------------
%% USER INPUT
clear all;
init();
% Case name
GiDFileList = readGiDFiles();
[indx,tf]   = listdlg('PromptString','Select a data file:','SelectionMode','single','ListString',GiDFileList);
caseName = GiDFileList{indx};

%Request user sample call selection
sample_list = {'randomUniform','randomNormal','latinHyperCube', 'quasiMonteCarloHalton'};
[indx,tf] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',sample_list);
% Call different functions depending on user's selection
switch indx
    case 1
        samplingCall = 'randomUniform';   
    case 2
        samplingCall = 'randomNormal';   
    case 3
        samplingCall = 'latinHyperCube';   
    case 4
        samplingCall = 'quasiMonteCarloHalton';
end

% Input function
input_list      = readInputFunctions();
[indx,tf]       = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',input_list);
% Call different functions depending on user's selection
inputFunction   = eval(['@',input_list{indx}]);

% Output function
outputFunction = @(paramStruct, up, FComplete) [                        ...
    drag( paramStruct, up, FComplete ),                                 ...
    lift( paramStruct, up, FComplete )                                  ...
    ];

% Number of samples and input distribution function
prompt = {'Enter input distribution average: ','Enter input distribution deviation: ', 'Enter number of samples: '};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0','1','10'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
avg = str2double(answer(1));
deviation = str2double(answer(2));
nSample = str2double(answer(3));

%% SAMPLING
switch samplingCall
    case 'randomUniform'
        % Uniform distribution with bounds
        distributionFunction = @(varargin) 2*deviation*rand(nSample,1)+avg-deviation;
    case 'randomNormal'
        % Random normal distribution sampling
        distributionFunction = @(varargin) normrnd(avg, deviation, nSample, 1);             
    case 'latinHyperCube'
        % Latin hypercube sampling
        distributionFunction =  @(varargin) lhsnorm(avg, (deviation*deviation), nSample); 
    case 'quasiMonteCarloHalton'
        % Quasi Monte Carlo Halton Sequence sampling
        p = haltonset(1,'Skip',1e3,'Leap',1e2);
        p = scramble(p,'RR2');
        haltonvector = net(p,nSample); % Halton Sequence sampling
        distributionFunction = @(varargin) norminv(haltonvector, avg, deviation); % Invert uniform Halton sequence to normal distribution
    otherwise
        % Error
        error('Invalid input sampling method')
end

% Generate inputs
input = distributionFunction();

%% RUN CALCULATION --------------------------------------------------------
[output, problemStruct] = monteCarlo(                                   ...
    caseName,                                                           ...
    input,                                                              ...
    inputFunction,                                                      ...
    outputFunction                                                      ...
    );

%% POSTPROCESSING ---------------------------------------------------------
% Plot mesh and nodes with boundary conditions ('fldMsh','homDBC','inhomDBC')
inputCaseToPlot     = 1;                                                    % <--- specify which input value should be applied
pstr                = ParameterStructure(problemStruct);                    % <--- create the object to plot
inputFunction(pstr,input(inputCaseToPlot));                                 % <--- apply the input function with the specified input
plotProblem(pstr,'homDBC','inhomDBC');                                      % <--- plot the desired features

% Convert outputs
if isa(output, 'cell') 
    output = cell2mat(output);
    outputVarNum = size(output,2);
end

% Get rid of results with singular matrices
sIndices            = singularIndices(output(:,1));
input(sIndices)     = [];
output(sIndices,:)  = [];

% Input histogram
figure
subplot(2,outputVarNum,1:outputVarNum)
histogram(input)
title('INPUT HISTOGRAM')

% Plot histograms of outputs
for k=1:outputVarNum
    subplot(2, outputVarNum,outputVarNum+k)
    histogram(output(:,k))
    title(['OUTPUT #',num2str(k),' HISTOGRAM'])
end