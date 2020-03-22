function [denom,gamma] = gradientDescent(propALE,p1_hist,p2_hist,p3_hist,....
                            djd1,djd2,djd3,djdp1,djdp2,djdp3,u_flag)
%% Licensing
%
%  License:         BSD License
%                   cane Multiphysics default license: cane/license.txt
%
%  Main authors:    Matthew Keller
%
%% Function documentation
%
%  Sets gradient descent parameters for Barzilai-Borwein gradient descent
%
%               Input :
%             propALE : Properties regarding the ALE boundary
%                         .nodes : The sequence of the nodal coordinates
%                                  on the ALE boundary
%                     .fcthandle : Function handle to the computation of
%                                  the ALE motion
%                       propUser : Extra user-defined parameters
%          p1_hist... : torage container of design parameter values
%             djd1... : Drag error of perturbed solution
%            djdp1... : Storage container of perturbation error
%              u_flag : Flag for ALE mesh motion 
%
%              Output :
%               denom : Barzilai-Borwein gradient descent parameter
%               gamma : Barzilai-Borwein gradient descent parameter

%% Function main body

if u_flag == 1 % Height
    denom = (djd1 - djdp1(end))^2;
    gamma = ((propALE.propUser.p1 - p1_hist(end)) * (djd1 - djdp1(end)))/denom;
elseif u_flag == 2 % Width
    denom = (djd2 - djdp2(end))^2;
    gamma = ((propALE.propUser.p2 - p2_hist(end)) * (djd2 - djdp2(end)))/denom;
elseif u_flag == 3 % Height and Width
    denom = (djd1 - djdp1(end))^2 + (djd2 - djdp2(end))^2;
    gamma = ((propALE.propUser.p1 - p1_hist(end)) * (djd1 - djdp1(end)) + (propALE.propUser.p2 - p2_hist(end)) * (djd2 - djdp2(end)))/denom;  
elseif u_flag == 4 % Taper
    denom = (djd3 - djdp3(end))^2;    
    gamma = ((propALE.propUser.p3 - p3_hist(end)) * (djd3 - djdp3(end)))/denom;   
elseif u_flag == 5 % Height and Taper
    denom = (djd1 - djdp1(end))^2 + (djd3 - djdp3(end))^2;    
    gamma = ((propALE.propUser.p1 - p1_hist(end)) * (djd1 - djdp1(end)) + (propALE.propUser.p3 - p3_hist(end)) * (djd3 - djdp3(end)))/denom;   
end