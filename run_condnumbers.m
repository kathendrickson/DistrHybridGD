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
% Description: Run file that adjusts the parameters beta and K and compares
%              the results. The same taumax and taumin are used for all
%              scenarios by using the most conservative requirement.
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
condn=5;
v = (condn-1).*rand(numAgents,1) + 1;
v(1) = 1;
v(2) = condn;
A = transpose(Q)*diag(v)*Q;
b=[1,2,3,-1,2];
beta = min(eig(A));
K = max(eig(A));

% Set Global Parameters Based on Most Conservative Estimate Above
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

%% Call HyEQsolver.m for 3 Condition Numbers
% Beta = 1, K=5
[t_15 j_15 x_15] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

% Beta=1, K=1
condn=1;
v = (condn-1).*rand(numAgents,1) + 1;
v(1) = 1;
v(2) = condn;
A = transpose(Q)*diag(v)*Q;
converged=zeros(numAgents,1);

[t_11 j_11 x_11] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

% Beta=1, K=2
condn=2;
v = (condn-1).*rand(numAgents,1) + 1;
v(1) = 1;
v(2) = condn;
A = transpose(Q)*diag(v)*Q;
converged=zeros(numAgents,1);

[t_12 j_12 x_12] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

% Beta=2, K=4
v = (4-2).*rand(numAgents,1) + 2;
v(1) = 4;
v(2) = 2;
A = transpose(Q)*diag(v)*Q;
converged=zeros(numAgents,1);

[t_24 j_24 x_24] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

% Beta=5, K=5
v = 5*ones(5,1);
A = transpose(Q)*diag(v)*Q;
converged=zeros(numAgents,1);

[t_55 j_55 x_55] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)


%% Create distance between iteration arrays. 
convdiff_11 = ones(size(t_11))*0.5;
for i = 2:size(t_11)
    convdiff_11(i) = norm(x_11(i,numAgents+1:2*numAgents) - x_11(i-1,numAgents+1:2*numAgents),2);
    if convdiff_11(i) == 0
        convdiff_11(i) = convdiff_11(i-1);
    end
end

convdiff_15 = ones(size(t_15))*0.5;
for i = 2:size(t_15)
    convdiff_15(i) = norm(x_15(i,numAgents+1:2*numAgents) - x_15(i-1,numAgents+1:2*numAgents),2);
    if convdiff_15(i) == 0
        convdiff_15(i) = convdiff_15(i-1);
    end
end

convdiff_12 = ones(size(t_12))*0.5;
for i = 2:size(t_12)
    convdiff_12(i) = norm(x_12(i,numAgents+1:2*numAgents) - x_12(i-1,numAgents+1:2*numAgents),2);
    if convdiff_12(i) == 0
        convdiff_12(i) = convdiff_12(i-1);
    end
end

convdiff_24 = ones(size(t_24))*0.5;
for i = 2:size(t_24)
    convdiff_24(i) = norm(x_24(i,numAgents+1:2*numAgents) - x_24(i-1,numAgents+1:2*numAgents),2);
    if convdiff_24(i) == 0
        convdiff_24(i) = convdiff_24(i-1);
    end
end

convdiff_55 = ones(size(t_55))*0.5;
for i = 2:size(t_55)
    convdiff_55(i) = norm(x_55(i,numAgents+1:2*numAgents) - x_55(i-1,numAgents+1:2*numAgents),2);
    if convdiff_55(i) == 0
        convdiff_55(i) = convdiff_55(i-1);
    end
end


%% Plot
figure(1); 
semilogy(t_11, convdiff_11, t_15, convdiff_15, t_12, convdiff_12,t_24, convdiff_24, t_55, convdiff_55)
legend('\beta=1, K=1','\beta=1, K=5','\beta=1, K=2','\beta=2, K=4','\beta=5, K=5' )
title('Effect of Problem Parameters on Convergence');
xlabel('Time'); ylabel('Distance Between Iterations')