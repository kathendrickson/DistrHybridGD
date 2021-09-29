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
% Description: Run file that tests the importance of agents knowing the
%              initial values of all other agents. This only impacts the
%              first jump. Inputs were designed so the first jump would
%              be expansive - agents thought that others were closer to
%              the optimum than they actually were.
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
v = 5*ones(5,1);
A = transpose(Q)*diag(v)*Q;
b=[    2.2279;
    3.2316;
    3.5468;
    3.7734;
    1.3801];
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
tauo = taumin + rand(1)*(taumax-taumin);


%% HyEqsolver Inputs
% simulation horizon
TSPAN = [0 16];                                                                 
JSPAN = [0 20];                                                                 
                                                                        
% rule for jumps                                                        
% rule = 1 -> priority for jumps                                        
% rule = 2 -> priority for flows                                        
% rule = 3 -> no priority, random selection when simultaneous conditions
rule = 1;
                                                                        
%solver tolerances
RelTol = 10e-8;
MaxStep = taumax/10;
options = odeset('RelTol', RelTol, 'MaxStep', MaxStep);

%% Call HyEQsolver.m for Two Tests
% Agents all know initialization values
z2o = z1o;
xi0 = [z1o;z2o;tauo];
[t_1,j_1,x_1] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

% Agents have different initial values
z2o = [-.2;-.4;-.599999999999999; .2; -.4];
xi0 = [z1o;z2o;tauo];

[t_2,j_2,x_2] = HyEQsolver(@f,@g,@C,@D,xi0,TSPAN,JSPAN,rule,options);
disp(converged)

zActual = [-0.445580000000000  -0.646320000000000  -0.709359999999999  -0.754679999999999  -0.276020000000001 -0.445580000000000  -0.646320000000000  -0.709359999999999  -0.754679999999999  -0.276020000000001];

error_1 = zeros(size(t_1));
for i=1:size(t_1)
    error_1(i) = norm(x_1(i,1:2*numAgents) - zActual,2);
end

error_2 = zeros(size(t_2));
for i=1:size(t_2)
    error_2(i) = norm(x_2(i,1:2*numAgents) - zActual,2);
end

%% Plot
figure(1); 
plot(t_1,error_1, t_2, error_2)
% x1 = xline(.0665, '--', 'Jump 1'); x1.LabelVerticalAlignment = 'bottom';
% x2 = xline(.0665*2, '--', 'Jump 2'); x2.LabelVerticalAlignment = 'bottom';
% x3 = xline(.0665*3, '--', 'Jump 3'); x3.LabelVerticalAlignment = 'bottom';
% x4 = xline(.0665*4, '--', 'Jump 4'); x4.LabelVerticalAlignment = 'bottom';
% x5 = xline(.0665*5, '--', 'Jump 5'); x5.LabelVerticalAlignment = 'bottom';

legend('Trial 1','Trial 2')
title('Comparing Results for Different Initial Values');
xlabel('Time'); ylabel('Distance from z^*')

%% Plot Using Toolbox
minarc = min([length(x_1),length(x_2)]);
t = [t_1(1:minarc),t_2(1:minarc)];
j = [j_1(1:minarc),j_2(1:minarc)];
x = [error_1(1:minarc),error_2(1:minarc)];

figure(2)
plotHarc(t,j,x);
title('Comparing Results for Different Initial Values');
xlabel('Time'); ylabel('Distance from z^*')