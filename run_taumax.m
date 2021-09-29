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
% Description: Run file to test different values of taumin and taumax on
%              the same problem.
%--------------------------------------------------------------------------

clear all                                                               

%% Global Parameters
global taumax
global taumin
global numAgents
global b
global A
global converged

% Problem Parameters
numAgents = 5;  %number of agents
%[Q,R] = qr(rand(numAgents));
Q = [   -0.547881311107979   0.365148780607783   0.076187197753691   0.745477022664536  -0.070370139526008;
        -0.375970848899789  -0.240036200120591  -0.635134609033948  -0.151610031301445  -0.612084120219850;
        -0.074203446272804  -0.051260302547968  -0.667299088026472   0.107812497819640   0.731405921754028;
        -0.061593929281199   0.871767383394488  -0.160481956244890  -0.458143735816345  -0.024034894992204;
        -0.741061418610573  -0.215506383396257   0.346059255250580  -0.446944280346138   0.291322458254830];
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 2;
v(2) = 4;
A = transpose(Q)*diag(v)*Q;
b=[1,2,3,-1,2];
beta = min(eig(A));
K = max(eig(A));

% Global Parameters
taumax=beta^2/(3*K^3+1);   %max value for timers
taumin=taumax/2;   %min value for timers
converged=zeros(numAgents,1);

%% Initial Conditions
% Initial conditions must begin within boundaries otherwise flow will never
% occur.
z1o = [2;2;2;2;2];
z2o = z1o;
tauo = taumax;
xi0 = [z1o;z2o;tauo];


%% HyEqsolver Inputs
% simulation horizon
TSPAN = [0 10];                                                                 
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

%% Call HyEQsolver.m for 3 Condition Numbers
% tau_max default
[t_5 j_5 x_5] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

%divided by 2
taumax = taumax/2;
taumin = taumax/2;

converged=zeros(numAgents,1);


[t_3 j_3 x_3] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

% divided by 4
taumax = taumax/2;
taumin = taumax/2;
converged=zeros(numAgents,1);

[t_1 j_1 x_1] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

%% Create distance between iteration arrays. 
convdiff_1 = ones(size(t_1))*0.5;
for i = 2:size(t_1)
    convdiff_1(i) = norm(x_1(i,numAgents+1:2*numAgents) - x_1(i-1,numAgents+1:2*numAgents),2);
    if convdiff_1(i) == 0
        convdiff_1(i) = convdiff_1(i-1);
    end
end

convdiff_3 = ones(size(t_3))*0.5;
for i = 2:size(t_3)
    convdiff_3(i) = norm(x_3(i,numAgents+1:2*numAgents) - x_3(i-1,numAgents+1:2*numAgents),2);
    if convdiff_3(i) == 0
        convdiff_3(i) = convdiff_3(i-1);
    end
end

convdiff_5 = ones(size(t_5))*0.5;
for i = 2:size(t_5)
    convdiff_5(i) = norm(x_5(i,numAgents+1:2*numAgents) - x_5(i-1,numAgents+1:2*numAgents),2);
    if convdiff_5(i) == 0
        convdiff_5(i) = convdiff_5(i-1);
    end
end

%% Plot
figure(1); 
semilogy(t_1, convdiff_1, t_3, convdiff_3, t_5, convdiff_5)
legend('smallest','medium','orig' )
title('Effect of Condition Number on Convergence');
% 
% A = beta^2*(1-2*taumax*K^3);
% B = 1-2*taumax*K;
% eterm = -(beta*A*B)/(8*K^2)