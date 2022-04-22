%--------------------------------------------------------------------------
%
% Project: Exponentially Converging Distributed Gradient Descent with 
%          Intermittent Communication via Hybrid Methods
%
% Author: Kat Hendrickson
%
% Utilizes the HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
%--------------------------------------------------------------------------
% Description: Run file that compares results for 5, 100, 500, 1000, and 
%              5000 agents.
%--------------------------------------------------------------------------

clear all                                                               

%% Global Parameters
global taumax
global taumin
global numAgents
global b
global A
global converged

%% Setup for 5 Agents
% Problem Parameters
numAgents = 5;  %number of agents
[Q,R] = qr(rand(numAgents));
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 2;
v(2) = 4;
A = transpose(Q)*diag(v)*Q;
b=rand(numAgents,1)*5;
beta = min(eig(A));
K = max(eig(A));

% Global Parameters
taumax=beta^2/(3*K^3+1);   %max value for timers
taumin=taumax/2;   %min value for timers
converged=zeros(numAgents,1);

%% Initial Conditions
% Initial conditions must begin within boundaries otherwise flow will never
% occur.
z1o = 2*ones(numAgents,1);
z2o = z1o;
tauo = taumin + rand(1)*(taumax-taumin);
xi0 = [z1o;z2o;tauo];


%% HyEqsolver Inputs
% simulation horizon
TSPAN = [0 16];                                                                 
JSPAN = [0 15000];                                                                 
                                                                        
% rule for jumps                                                        
% rule = 1 -> priority for jumps                                        
% rule = 2 -> priority for flows                                        
% rule = 3 -> no priority, random selection when simultaneous conditions
rule = 1;
                                                                        
%solver tolerances
RelTol = 10e-8;
MaxStep = taumax/10;
options = odeset('RelTol', RelTol, 'MaxStep', MaxStep);

%% Call HyEQsolver.m for 5 Agents
[t_5 j_5 x_5] = HyEQsolver(@f,@f,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

convdiff_5 = ones(size(t_5))*0.5;
for i = 2:size(t_5)
    convdiff_5(i) = norm(x_5(i,numAgents+1:2*numAgents) - x_5(i-1,numAgents+1:2*numAgents),2);
    if convdiff_5(i) == 0
        convdiff_5(i) = convdiff_5(i-1);
    end
end

%% Setup for 100 Agents
% Problem Parameters
numAgents = 100;  %number of agents
[Q,R] = qr(rand(numAgents));
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 2;
v(2) = 4;
A = transpose(Q)*diag(v)*Q;
b=rand(numAgents,1)*5;
beta = min(eig(A));
K = max(eig(A));

% Global Parameters
taumax=beta^2/(3*K^3+1);   %max value for timers
taumin=taumax/2;   %min value for timers
converged=zeros(numAgents,1);

% Initial Conditions
% Initial conditions must begin within boundaries otherwise flow will never
% occur.
z1o = 2*ones(numAgents,1);
z2o = z1o;
tauo = taumin + rand(1)*(taumax-taumin);
xi0 = [z1o;z2o;tauo];

%% Call HyEQsolver.m for 100 Agents
[t_100 j_100 x_100] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

convdiff_100 = ones(size(t_100))*0.5;
for i = 2:size(t_100)
    convdiff_100(i) = norm(x_100(i,numAgents+1:2*numAgents) - x_100(i-1,numAgents+1:2*numAgents),2);
    if convdiff_100(i) == 0
        convdiff_100(i) = convdiff_100(i-1);
    end
end

%% Setup for 500 Agents
% Problem Parameters
numAgents = 500;  %number of agents
[Q,R] = qr(rand(numAgents));
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 2;
v(2) = 4;
A = transpose(Q)*diag(v)*Q;
b=rand(numAgents,1)*5;
beta = min(eig(A));
K = max(eig(A));

% Global Parameters
taumax=beta^2/(3*K^3+1);   %max value for timers
taumin=taumax/2;   %min value for timers
converged=zeros(numAgents,1);

% Initial Conditions
% Initial conditions must begin within boundaries otherwise flow will never
% occur.
z1o = 2*ones(numAgents,1);
z2o = z1o;
tauo = taumin + rand(1)*(taumax-taumin);
xi0 = [z1o;z2o;tauo];

%% Call HyEQsolver.m for 500 Agents
[t_500 j_500 x_500] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

convdiff_500 = ones(size(t_500))*0.5;
for i = 2:size(t_500)
    convdiff_500(i) = norm(x_500(i,numAgents+1:2*numAgents) - x_500(i-1,numAgents+1:2*numAgents),2);
    if convdiff_500(i) == 0
        convdiff_500(i) = convdiff_500(i-1);
    end
end

%% Setup for 1000 Agents
% Problem Parameters
numAgents = 1000;  %number of agents
[Q,R] = qr(rand(numAgents));
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 2;
v(2) = 4;
A = transpose(Q)*diag(v)*Q;
b=rand(numAgents,1)*5;
beta = min(eig(A));
K = max(eig(A));

% Global Parameters
taumax=beta^2/(3*K^3+1);   %max value for timers
taumin=taumax/2;   %min value for timers
converged=zeros(numAgents,1);

% Initial Conditions
% Initial conditions must begin within boundaries otherwise flow will never
% occur.
z1o = 2*ones(numAgents,1);
z2o = z1o;
tauo = taumin + rand(1)*(taumax-taumin);
xi0 = [z1o;z2o;tauo];

%% Call HyEQsolver.m for 1000 Agents
[t_1000 j_1000 x_1000] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

convdiff_1000 = ones(size(t_1000))*0.5;
for i = 2:size(t_1000)
    convdiff_1000(i) = norm(x_1000(i,numAgents+1:2*numAgents) - x_1000(i-1,numAgents+1:2*numAgents),2);
    if convdiff_1000(i) == 0
        convdiff_1000(i) = convdiff_1000(i-1);
    end
end

%% Setup for 5000 Agents
% Problem Parameters
numAgents = 5000;  %number of agents
[Q,R] = qr(rand(numAgents));
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 2;
v(2) = 4;
A = transpose(Q)*diag(v)*Q;
b=rand(numAgents,1)*5;
beta = min(eig(A));
K = max(eig(A));

% Global Parameters
taumax=beta^2/(3*K^3+1);   %max value for timers
taumin=taumax/2;   %min value for timers
converged=zeros(numAgents,1);

% Initial Conditions
% Initial conditions must begin within boundaries otherwise flow will never
% occur.
z1o = 2*ones(numAgents,1);
z2o = z1o;
tauo = taumin + rand(1)*(taumax-taumin);
xi0 = [z1o;z2o;tauo];

%% Call HyEQsolver.m for 5000 Agents
[t_5000 j_5000 x_5000] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

convdiff_5000 = ones(size(t_5000))*0.5;
for i = 2:size(t_5000)
    convdiff_5000(i) = norm(x_5000(i,numAgents+1:2*numAgents) - x_5000(i-1,numAgents+1:2*numAgents),2);
    if convdiff_5000(i) == 0
        convdiff_5000(i) = convdiff_5000(i-1);
    end
end

%% Plot
figure(1); 
semilogy(t_5, convdiff_5, t_100, convdiff_100, t_500, convdiff_500, t_1000, convdiff_1000, t_5000, convdiff_5000)
legend('5 Agents','100 Agents', '500 Agents', '1000 Agents', '5000 Agents')
title('Effect of Network Size on Convergence Rate');
xlabel('Time'); ylabel('Distance Between Iterations');
xlim([0,10]);
