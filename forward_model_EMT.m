function[cells] = forward_model_EMT(Ce_init, params, t)

ge = params(1); % growth rate of epithelial cells
gm = params(2); % growth rate of mesenchymal cells
kem = params(3); % baseline transition rate E to M
kme = params(4); % transition rate from M to E
carcap = params(5);
YAPnuc = params(6);

f = @(t,Cc) [ge*(1-((Cc(1)+ Cc(2))./carcap))*Cc(1)- kem*YAPnuc.*Cc(1) + kme*Cc(2);  % dE/dt
             gm*(1-((Cc(1)+ Cc(2))./carcap))*Cc(2)+ kem*YAPnuc.*Cc(1) - kme*Cc(2)]; % dM/dt
options1 = odeset('Refine',1);  
options = odeset(options1,'NonNegative',1:2);
[tout,Cc]=ode45(f, t,Ce_init, options);

cells = Cc;
end