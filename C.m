function v  = C(x) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Kat Hendrickson 
%
% Project: Exponentially Converging Distributed Gradient Descent with 
%          Intermittent Communication via Hybrid Methods
%
% Name: C.m
%
% Description: Flow set
%
% Dependencies: 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% timer state
tau = x(end);

% flow set
if tau > 0
    v = 1;  % report flow
else 
    v = 0;  % do not report flow
end