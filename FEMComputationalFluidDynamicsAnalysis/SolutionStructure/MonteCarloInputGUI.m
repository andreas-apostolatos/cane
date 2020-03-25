function [samplingCall,Umax,iterationLimit,learning_rate,...
          delta_p1,delta_p2,delta_p3,convergenceLimit,h_limit,w_limit,t_limit] = MonteCarloInputGUI(analysis_type)
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
%    convergenceLimit : Criteria to exit optimization loop based on the 
%                       minimum value of the djdp finite difference step
%          h_limit... : Limit of structure height and width

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
        t_limit = 0;
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
        t_limit = 0;
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
        t_limit = 0;
    case 'Taper'
        prompt = {'Enter perturbation of upper boundary width: ', 'Enter upper boundary width limit: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'2e-5', '0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p3 = str2double(answer(1));
        delta_p1 = 0;
        delta_p2 = 0;
        h_limit = 0;
        w_limit = 0;
        t_limit = str2double(answer(2));
    case 'Height and Taper'
        prompt = {'Enter perturbation of height: ', 'Enter perturbation of upper boundary width: ', 'Enter height limit: ', 'Enter upper boundary width limit: '};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'3e-4', '2e-5', '0.01', '0.01'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        delta_p1 = str2double(answer(1));
        delta_p3 = str2double(answer(2));
        delta_p2 = 0;
        h_limit = str2double(answer(3));
        w_limit = 0;
        t_limit = str2double(answer(4));
    otherwise
        error('Invalid input sampling method')
end

if strcmp(analysis_type, 'simplified')
    % Request general input definitions from user
    prompt = {'Enter max input velocity: ','Enter iteration limit: ', 'Enter learning rate: ', 'Enter convergence limit: '};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'0.1','20', '0.01','1e-4'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Umax = str2double(answer(1));
    iterationLimit = str2double(answer(2));
    learning_rate = str2double(answer(3));
    convergenceLimit = str2double(answer(4));
    
elseif strcmp(analysis_type, 'full_descent')
    % Request general input definitions from user
    prompt = {'Enter max input velocity: ','Enter iteration limit: ', 'Enter design penalization: ', 'Enter convergence limit: '};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'0.1','20','1e4','1e-4'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    Umax = str2double(answer(1));
    iterationLimit = str2double(answer(2));
    learning_rate = str2double(answer(3));
    convergenceLimit = str2double(answer(4));
end    
end