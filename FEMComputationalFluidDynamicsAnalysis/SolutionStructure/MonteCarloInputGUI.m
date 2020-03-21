function [samplingCall, Umax, iterationLimit, design_penalization, learning_rate, delta_p1, delta_p2, delta_p3, h_limit, w_limit] = MonteCarloInputGUI()
%% Licensing
%
%  License:         BSD License
%                   cane Multiphysics default license: cane/license.txt
%
%  Main authors:    Matthew Keller
%
%% Function documentation
%
%  User input GUI for gradient descent process
%
%              Output :
%        samplingCall : User defined selection of analysis type
%                Umax : Max velocity from research paper
%      iterationLimit : Limit of iteration loops
% design_penalization : Penalization of the gradient descent step:
%                       - Barzilai-Borwein
%       learning_rate : Penalization of the gradient descent step:
%                       - Basic descent
%         delta_p1... : Perturbation of each design parameters
%    h_limit, w_limit : Limit of structure height and width

%% Function main body

% Request user sample call selection
sample_list = {'Height','Width','Height and Width', 'Taper', 'Height and Taper'};
[indx,tf] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',sample_list);
% Call different functions depending on user's selection
switch indx
    case 1
        samplingCall = 'Height';   
    case 2
        samplingCall = 'Width';   
    case 3
        samplingCall = 'Height and Width';   
    case 4
        samplingCall = 'Taper';
    case 5
        samplingCall = 'Height and Taper';
end

% Request user input dependent on sample call selection
switch samplingCall
    case 'Height'
        prompt = {'Enter perturbation of height: ', 'Enter height limit: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'3e-4', '0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p1 = str2double(answer(1));
        delta_p2 = 0;
        delta_p3 = 0;
        h_limit = str2double(answer(2));
        w_limit = 0;
    case 'Width'
        prompt = {'Enter perturbation of width: ', 'Enter width limit: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'2e-5', '0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p2 = str2double(answer(1));
        delta_p1 = 0;
        delta_p3 = 0;
        h_limit = 0;
        w_limit = str2double(answer(2));
    case 'Height and Width'
        prompt = {'Enter perturbation of height: ', 'Enter perturbation of width: ', 'Enter height limit: ', 'Enter width limit: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'3e-4', '2e-5','0.01','0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p1 = str2double(answer(1));
        delta_p2 = str2double(answer(2));
        delta_p3 = 0;
        h_limit = str2double(answer(3));
        w_limit = str2double(answer(4));
    case 'Taper'
        prompt = {'Enter perturbation of taper: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p3 = str2double(answer(1));
        delta_p1 = 0;
        delta_p2 = 0;
        h_limit = 0;
        w_limit = 0;
    case 'Height and Taper'
        prompt = {'Enter perturbation of height: ', 'Enter perturbation of taper: ', 'Enter height limit: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'3e-4', '0.01', '0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p1 = str2double(answer(1));
        delta_p3 = str2double(answer(2));
        delta_p2 = 0;
        h_limit = str2double(answer(3));
        w_limit = 0;
    otherwise
        % Error
        error('Invalid input sampling method')
end

% Request general input definitions from user
prompt = {'Enter max input velocity: ','Enter iteration limit: ', 'Enter design penalization: ', 'Enter learning rate: '};
dlgtitle = 'Input';
dims = [1 35];
definput = {'0.1','20','1e4', '0.01'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Umax = str2double(answer(1)); % Max input velocity defined in the reference paper
iterationLimit = str2double(answer(2)); % Limit search iterations
design_penalization = str2double(answer(3)); % Design penalty used in J calculations
learning_rate = str2double(answer(4));
end