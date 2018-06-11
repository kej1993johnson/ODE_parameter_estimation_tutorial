function dy=forward_Cc(t,y,para)
dy=zeros(2,1);

ge = para(1); % growth rate of epithelial cells
gm = para(2); % growth rate of mesenchymal cells
kem = para(3); % baseline transition rate E to M
kme = para(4); % transition rate from M to E
carcap = para(5)*1e6; % carrying capacity
phi = 0.8;

dy(1)=ge*y(1)*(1-(y(1)+y(2))/carcap)-kem*phi*y(1)+kme*y(2);
dy(2)=gm*y(2)*(1-(y(1)+y(2))/carcap)+kem*phi*y(1)-kme*y(2);
end