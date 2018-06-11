    RandNormD=normrnd(0,1,73,2);
    % add nosie to the generated data
    Ccm(:,1:2)=(RandNormD(:,1:2).*Cc(:,1:2)*0.05/1.96+Cc(:,1:2));
    Ccm=round(Ccm);
    
    y1=Ccm(:,1);
    y2=Ccm(:,2);
    save('Em','y1');
    save('Mm','y2');
    