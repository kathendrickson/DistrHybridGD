function v  = D(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Kat Hendrickson
%
% Project: Exponentially Converging Distributed Gradient Descent with 
%          Intermittent Communication via Hybrid Methods
%
% Name: D.m
%
% Description: Jump set
%
% Dependencies: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% timer state
tau = x(end);

% jump set
if tau <= 0 % jump condition
    v = 1;  % report jump
else
    v = 0;  % do not report jump
end