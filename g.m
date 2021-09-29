function xplus = g(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Kat Hendrickson
%
% Project: Exponentially Converging Distributed Gradient Descent with 
%          Intermittent Communication via Hybrid Methods
%
% Name: g.m
%
% Description: Jump map
%
% Dependencies: run.m, grad.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global taumax
global taumin
global numAgents

% states
z1 = x(1:numAgents);
%z2 = x(numAgents+1:numAgents*2);
tau = x(end);

%% Jump Map

% Communication and Timer Reset
if tau <= 0
    z1plus = z1;
    z2plus = z1;
    tauplus = taumin + rand(1)*(taumax-taumin);
end

% Assembles state vector, xplus
xplus=[z1plus;z2plus;tauplus];
    
end