function gradient = grad(z2,agent)
%PRIMAL GRADIENT Returns a gradient given an agent's copy of x.
%   Inputs: (x, agent)
%       x - agent's complete copy of x.
%       agent - agent number; used to retrieve its copy of x_i. 
%   Returns: gradient
%       gradient - scalar gradient for defined problem

%% Global variables:
global A
global b

%% Define gradient for x_i:

gradient= A(agent,:)*z2 + b(agent);

end

