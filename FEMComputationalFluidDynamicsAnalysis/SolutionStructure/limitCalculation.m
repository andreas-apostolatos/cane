function [limit_flag] = limitCalculation(p1,p2,p3,h_limit,w_limit,t_limit,u_flag)
%% Licensing
%
%  License:         BSD License
%                   cane Multiphysics default license: cane/license.txt
%
%  Main authors:    Matthew Keller
%
%% Function documentation
%
%  Flags limits of design boundaries set by user if violated
%
%               Input :
%            p1,p2,p3 : Internal variables at current step
%          h_limit... : Limit of structure height and width
%              u_flag : Flag for ALE mesh motion 
%
%              Output :
%          limit_flag : Flag of potential violated conditions (true/false)

if u_flag == 1 % Limit analysis dependent on height
    if p1 <= h_limit
        limit_flag = true;
    else
        limit_flag = false;
    end
    
elseif u_flag == 2 % Limit analysis dependent on width
    if p2 <= w_limit
        limit_flag = true;
    else
        limit_flag = false;
    end
     
elseif u_flag == 3 % Limit analysis dependent on height and width
    if p1 <= h_limit
        limit_flag = true;
    elseif p2 <= w_limit
        limit_flag = true;
    else
        limit_flag = false;
    end 
        
elseif u_flag == 4 % Limit analysis dependent on taper ratio
    if p3 <= t_limit
        limit_flag = true;
    else
        limit_flag = false;
    end    
    
elseif u_flag == 5 % Limit analysis dependent on height and taper
    if p1 <= h_limit
        limit_flag = true;
    elseif p3 <= t_limit
        limit_flag = true;
    else
        limit_flag = false;
    end       
end

