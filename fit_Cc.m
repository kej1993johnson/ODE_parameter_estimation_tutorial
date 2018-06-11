function y=fit_Cc(k,t)
global init;
% y0=iniCondition;
y0=init;
% define t here
[t Y]=ode45(@inv_Model,t,y0);

function dy=inv_Model(t,y)
dy=zeros(2,1);

ge = k(1); % growth rate of epithelial cells
gm = k(2); % growth rate of mesenchymal cells
kem = k(3); % baseline transition rate E to M
kme = k(4); % transition rate from M to E
carcap = k(5)*1e6; % carrying capacity
phi = 0.8;

dy(1)=ge*y(1)*(1-(y(1)+y(2))/carcap)-kem*phi*y(1)+kme*y(2);
dy(2)=gm*y(2)*(1-(y(1)+y(2))/carcap)+kem*phi*y(1)-kme*y(2);

end

y=Y(:,[1:2]);

end