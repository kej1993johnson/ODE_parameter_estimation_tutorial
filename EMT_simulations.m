% Example EMT Model for ODE/ Parameter Estimation Tutorial
close all; clear all; clc
% E and M subpopulation transition model
YAPnuc= [0.2 0.4 0.6 0.8];
%YAPnuc = 0.6;
ge = 0.5; % growth rate of epithelial cells
kem = 0.2; % baseline transition rate E to M
kme = 0.32; % transition rate from M to E
gm = 0.32; % growth rate of mesenchymal cells
carcap = 1e6;

Ce_init(1) = 1; % start with 1 epithelial cell
Ce_init(2) = 1; % start with 1 mesenchymal cell

params = vertcat(ge, kem, kme, gm, carcap, YAPnuc(4));
param_nms = {'g_{E}', 'k_{em}', 'k_{me}', 'g_{M}', '\theta', '\phi_{YAP}'} ;
tsamp = 0:1:72;

% first run it just once for YAPnuc = 0.8
f = @(t,Cc) [ge*(1-((Cc(1)+ Cc(2))./carcap))*Cc(1)- kem*YAPnuc(4).*Cc(1) + kme*Cc(2);  % dE/dt
             gm*(1-((Cc(1)+ Cc(2))./carcap))*Cc(2)+ kem*YAPnuc(4).*Cc(1) - kme*Cc(2)]; % dM/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[t,Cc]=ode45(f, tsamp,Ce_init, options);
namesc={'E', 'M'};
E= Cc(:,1);
M = Cc(:,2);

figure;
plot(tsamp, Cc(:,1), 'b', 'LineWidth', 2)
hold on
plot(tsamp, Cc(:,2), 'r', 'LineWidth',2)
xlabel ('time')
ylabel('cells')
legend ('E cells', 'M cells')
legend boxoff
title('E and M cell trajectories for \phi_{YAP}=0.8')

%% Now step through different YAP nuc fractions and observe the effect on 
% the fraction of M cells in time

figure;
for i = 1:length(YAPnuc)
f = @(t,Cc) [ge*(1-((Cc(1)+ Cc(2))./carcap))*Cc(1)- kem*YAPnuc(i).*Cc(1) + kme*Cc(2);  % dE/dt
             gm*(1-((Cc(1)+ Cc(2))./carcap))*Cc(2)+ kem*YAPnuc(i).*Cc(1) - kme*Cc(2)]; % dM/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[t,Cc]=ode45(f, tsamp,Ce_init, options);
namesc={'E', 'M'};
E= Cc(:,1);
M = Cc(:,2);




plot(tsamp, (Cc(:,2)./(Cc(:,1)+Cc(:,2))), 'LineWidth', 2)
text(tsamp(end), (Cc(end,2)/(Cc(end,1)+Cc(end,2))), ['\psi_{YAPnuc}=', num2str(YAPnuc(i))]);
hold on
xlabel('time')
ylabel('fraction of M cells')
title('Fraction of mesenchymal cells in time for varying \phi_{YAP}')
end

%% What is the effect of YAP nuclear localization on the fraction of M cells?
% Function to run the forward model
[cells]=forward_model_EMT(Ce_init, params, tsamp);
Eend= cells(end,1);
Mend = cells(end,2);
fracM = Mend./(Eend+Mend);

fracYAPvec = linspace(0, 1, 25); % make a vector of YAP fracs from 0 to 1;
p = params;
for i = 1:length(fracYAPvec)
    p(6) = fracYAPvec(i); % change YAP fraction each iteration
    [C]=forward_model_EMT(Ce_init, p, tsamp);
    fracMvec(i) = C(end,2)./(C(end,1)+C(end,2));
end

figure;
hold off
plot(fracYAPvec, fracMvec,'r', 'LineWidth',2)
xlabel('\phi_{YAP}')
ylabel('fraction of mesenchymal cells')
title('Fraction of mesenchymal cells as a function of YAP nuclear fraction')
%% Example map of parameter space to test effect of epithelial cell growth 
% rate and YAP nuclear fraction on the percent of cells in the M phenotype

% already have fracYAPvec from above, but lets update it to be higher
% resolution
fracYAPvec = linspace(0, 1, 50);
gevec = linspace(0.01*ge, ge, 50); % range from growth rate of E of 0 to growth rate of 5 x larger
p = params;

for i =1:length(gevec)
    p(6) = fracYAPvec(i);
    for j = 1:length(fracYAPvec)
        p(1) = gevec(j);
        [C]=forward_model_EMT(Ce_init, p, tsamp); % p is new each time
        fracmat(i,j) = C(end,2)./(C(end,1)+C(end,2));
    end
end
%%
xticklabels = linspace(0.01*ge, ge, 9);
xticks = linspace(1, length(gevec), numel(xticklabels));
yticklabels = linspace(0, 1, 9);
yticks = linspace(1, length(fracYAPvec), numel(yticklabels));


figure;
imagesc((fracmat))
colorbar
hold on
xlabel ('epithelial cell growth rate')
ylabel('\phi_{YAP}')
title('Effect of YAP nuclear fraction and E growth rate on fraction of M cells')
set(gca, 'XTick', xticks, 'XTickLabel', round(xticklabels,2))
set(gca, 'YTick', yticks, 'YTickLabel', round(yticklabels,2))
%% Sensitivity analysis
% Perform a sensitivity analysis on the 6 parameters

for j = 1:length(params)
    p = params;
  % first find for current p
[ Ci] = forward_model_EMT( Ce_init, p, tsamp );
% next find for incerased param
fraci = Ci(end,2)./(Ci(end,1)+Ci(end,2));
    p(j)=params(j)*1.01;
    deltap=0.01*params(j);
[ Cf] = forward_model_EMT( Ce_init, p, tsamp);
    fracf = Cf(end,2)./(Cf(end,1)+Cf(end,2));
sens_p(j) = norm(fracf-fraci)/deltap;
end
figure;
barh(1:1:length(p), (-log(sens_p)))
set(gca,'yticklabel',param_nms)
ylabel('Parameter')
xlabel('Sensitivity Score')
title('Sensitivity analysis on fraction of mesenchymal cells')

%% Maximum Likelihood Parameter Estimation
% Using simulated data from Jianchen
sigma = 1; % based on Jianchen's simulated data with standard deviation of 1
load Mm.mat; % this gets loaded in as y1
load Em.mat; % this gets loaded in as y2
ydata = horzcat(y1(2:end), y2(2:end));
t = tsamp(2:end); % don't fit on first time point, will assume we know the initial conditions

% Transform parameters and data into log space (and back)
pfxform = @(pval)[1 1 1 1 1 0].*log(pval)+[0 0 0 0 0 1].*log(pval./(1-pval)); %'forward' parameter transform into Reals
pbxform = @(phat)[1 1 1 1 1 0].*exp(phat)+[0 0 0 0 0 1].*(exp(phat)./(1+exp(phat)));  %'backward' parameter transform into model space
yfxform = @(y)log(y); % 'forward' transform for data and model output
ybxform = @(yhat)exp(yhat); % 'inverse' transform for data and model output

guess = [0.45 0.15 0.35 0.4 1.2e6 0.7]; % guess of parameters
modelfun = @(p)forward_model_EMT(Ce_init, p, t);

loglikelihood = @(phat)(sum(sum(log(normpdf(yfxform(ydata),yfxform(modelfun(pbxform(phat))), sigma)))));

% minimize objective function (-negative LL)
objfun = @(phat)-loglikelihood(phat);
phatbest = fminsearch(objfun, pfxform(guess)); % find best fitting parameters
params_fit = pbxform(phatbest);
%% Evaluate model fit

model_fit = modelfun(params_fit);

figure;
plot(t, model_fit(:,1), 'b-', 'LineWidth', 2)
hold on 
plot(tsamp, y1, 'bo')
plot(t, model_fit(:,2), 'r-', 'LineWidth',2)
plot(tsamp, y2, 'ro')
xlabel('time')
ylabel('cells')
legend('E cells model', 'E cells data', 'M cells model', 'M cells data', 'Location','NorthWest')
title('Model fit vs. simulated data using MLE')
xlim([0 72])

pct_error_params = (abs((params)-params_fit')./params)*100;

chi_squared = sum(sum(((ydata-model_fit).^2)./model_fit))

%% Write MCMC Search Algorithm (Simulated Annealing)
% For eacample where we only search for kme, kem, and phi YAP (params 2,3,
% and 6)
% Set the known values of guess
guess(1) = params(1);
guess(4) = params(4);
guess(5) = params(5);
J_init = objfun(pfxform(guess));
J_curr = J_init; % set current negLL value to inital value to initialize
count_accepted = 0; % start with no changes accepted

% Initialize a vector to store b and d params in
nruns = 3000;
store_params = zeros(nruns, 4); % store J + 3 params kem, kme, phiYAP
store_acc_params = zeros(nruns, 4 ); % store the converged parameters and J, + 3 params kem, kme, phiYAP

% Outside the loop: initialize temperature and iterations
T0= 2;
k = 1:1:nruns;
T= T0*exp(-5*k./(nruns)); % sets up gradual cooling, vary this function to accept or reject more stringently
% set the current value of parameters
guess1 = guess;

% check that temperature annealing is working
figure;
hold off
set(gca,'LineWidth',1.5,'FontSize',12);
plot(k, T, 'LineWidth', 3)
xlabel('Markov step')
ylabel('Temperature')
title('Temperature Annealing')

step = 0.005; % step size for searching for rate and growth parameters

% Set the first values of the transition rates and YAP nuclear localization
% as the initial guess you gave it before
kem1=guess1(2);
kme1=guess1(3);
YAPnuc1=guess1(6);
%%
for k = 1:nruns
    kem2 = kem1 + step*(2*rand-1);
    kme2 = kme1 + step*(2*rand-1);
    YAPnuc2 = YAPnuc1+ step*(2*rand-1);

    % Constrain search region
    if kem2<0 
        kem2=0;
    end
    if kem2>1
        kem2 = 1;
    end
    if kme2<0
        kme2=0;
    end
    if kme2>1
        kme2 = 1;
    end
    if YAPnuc2<0 
        YAPnuc2=0;
    end
    if YAPnuc2>1
        YAPnuc2 = 1;
    end

    % find the neg LL of the new params
    guess2 = guess;
    guess2(2)=kem2;
    guess2(3)=kme2;
    guess2(6)=YAPnuc2;
    J_new = objfun(pfxform(guess2));
    % store the  negLL and searched parameters
    store_params(k,1) = J_new;
    store_params(k,2) =kem2;
    store_params(k,3) = kme2;
    store_params(k,4) = YAPnuc2;

    prob(k) = exp((J_curr-J_new)./T(k));
    % if Jcurr > Jnew, Jnew is better--the numerator in the exponent is positive &
    % prob>1--> change will always be accepted
    % if Jcurr < Jnew, numerator is negative, & prob <1 --> chance change
    % will be accepted
    
    if rand < prob(k) % if true, accept the change!
        kem1 = kem2;
        kme1 = kme2;
        YAPnuc1=YAPnuc2;
        guess1=guess2;
        J_curr = J_new;
        count_accepted = count_accepted +1;
        % decrease search step size
        step = 0.999*step;
    else
        % increase search step size
        step = 1.001*step;
    end
    store_acc_params(k,1) = J_curr;
    store_acc_params(k,2)= kem1;
    store_acc_params(k,3)= kme1;
    store_acc_params(k,4)= YAPnuc1;

       
end
%% Check parameter search
p= [params(2) params(3) params(6)]
params_best_MC = store_acc_params(end,2:end)

figure;
plot(1:1:nruns, store_acc_params(:,1))
xlabel('iterations')
ylabel('negLL')
title('negLL')

 figure;
 plot(1:1:nruns, store_acc_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_acc_params(:,3), 'b.')
 plot(1:1:nruns, store_acc_params(:,4),'r.')
%   plot(1:1:nruns, store_acc_params(:,5),'m.')
 legend( 'kem', 'kme', '\phi_{YAP}')
 xlabel('iterations')
 ylabel('parameter')
 title('Accepted birth and death parameters')
 
 figure;
 plot(1:1:nruns, store_params(:,2), 'g.')
 hold on
 plot(1:1:nruns, store_params(:,3), 'b.')
 plot(1:1:nruns, store_params(:,4), 'r.')
 legend( 'kem', 'kme', '\phi_{YAP}')
 xlabel('iterations')
 ylabel('parameter')
 title('Explored birth and death parameters')
 
figure;
subplot(1,3,1)
hist(store_acc_params(:,2))
ylabel('frequency')
xlabel('kem')
subplot(1,3,2)
hist(store_acc_params(:,3))
ylabel('frequency')
xlabel('kme')
subplot(1,3,3)
hist(store_acc_params(:,4))
ylabel('frequency')
xlabel('YAPnuc')


% Plots of two parameters colored by likelihood
likelihoodval = exp(-store_acc_params(:,1));
likelihoodvalall = exp(-store_params(:,1));
pointsize = 30;
figure;
scatter(store_acc_params(:,2), store_acc_params(:,3), pointsize, likelihoodval,'filled');
colorbar
xlabel('kem')
ylabel('kme')
title(' kem versus kme colored by likelihood for accepted parameters')

