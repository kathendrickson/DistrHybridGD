function xdot = f(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Kat Hendrickson
%
% Project: Exponentially Converging Distributed Gradient Descent with 
%          Intermittent Communication via Hybrid Methods
%
% Name: f.m
%
% Description: Flow map
%
% Dependencies: run.m, grad.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global numAgents
global converged


% states
%z1 = x(1:numAgents);
z2 = x(numAgents+1:numAgents*2);
%tau = x(end);

%% Flow Map
xdot = zeros(size(x));
for i = 1:numAgents
    xidot = -grad(z2,i); 
    xdot(i)=xidot;
    if grad(z2,i) <= 10e-8
        converged(i)=1;
    end
end
xdot(end)=-1;

end