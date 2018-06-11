clc;
clear all
close all;
%%
% data with noise used as data to be fitted
load Mm.mat;
load Em.mat;
DataM=[y1 y2];

figure(); 
plot(tsamp,y1,'bo',tsamp,y2,'ro')
axis([0 72 0 7.5e5])
legend('M cells','E cells');
xlabel('tsamp')
ylabel('Cell number')
set(gca,'FontSize',15)

Ce_init(1) = 1;
Ce_init(2) = 1;

global init;
init=Ce_init;

%% with one initial guess
ParEst0=[2 2 2 2 5]; % initial guess
LB=[0 0 0 0 0]; % lower boundary
UB=[2 2 2 2 10]; % upper boundary
options = optimoptions('lsqcurvefit','MaxFunEvals',1200,'MaxIter',200,'display','off',...
                       'FinDiffRelStep',1e-9,'TolFun',1e-12,'TolX',1e-9); 
ParEst= lsqcurvefit(@fit_Cc,ParEst0,tsamp,DataM(:,1:2),LB,UB,options);

% run the forward model to check the fitting results
[t DataF]=ode45(@(t,y)forward_Cc(t,y,ParEst),tsamp,Ce_init);
% normalized/non-normalized residual sum of squares 
nRSS=sum(sum(((DataF(:,1:2)-DataM(:,1:2))./std(DataM(:,1:2))).^2));
RSS=sum(sum(((DataF(:,1:2)-DataM(:,1:2))).^2));
hold on;
plot(tsamp,DataF(:,1),'b');
plot(tsamp,DataF(:,2),'r');

%% with multiple start points
p0=[2 2 2 2 10];
LB=[0 0 0 0 0];
UB=[2 2 2 2 10];
% create the problem
problem = createOptimProblem('lsqcurvefit','x0',p0,'objective',@(p,xdata) fit_Cc(p,xdata)...
    ,'lb',LB,'ub',UB,'xdata',tsamp,'ydata',DataM(:,1:2));
problem.options.MaxFunctionEvaluations=2000;
ms = MultiStart('PlotFcns',@gsplotbestf,'Display','iter','StartPointsToRun','bounds')
% deciede number of initial guesses: 100
[xmulti,fval,exitflag,output] = run(ms,problem,100);
% exit with a positive flag is considered a success
success=output.localSolverSuccess;