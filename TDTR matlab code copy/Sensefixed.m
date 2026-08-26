for ii=1:length(lambda)
    %-------------Specific heat-----------------
    Ctemp=C;
    Ctemp(ii)=C(ii)*1.01;
    [deltaR_data_C(:,ii),ratio_data_C(:,ii)]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,Ctemp,t,eta,r_pump,r_probe,A_pump,nnodes);
    delta_C(:,ii)=ratio_data_C(:,ii)-ratio_data;
    Num=log(ratio_data_C(:,ii))-log(ratio_data);
    Denom=log(Ctemp(ii))-log(C(ii));
    S_C(:,ii)=Num/Denom;
    %-------------Thermal Conductivity (ky)----------
    lambdatemp=lambda;
    lambdatemp(ii)=lambda(ii)*1.01;
    etatemp = eta.*lambda./lambdatemp;
    [deltaR_data_L(:,ii),ratio_data_L(:,ii)]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambdatemp,C,t,etatemp,r_pump,r_probe,A_pump,nnodes);
    delta_L(:,ii)=ratio_data_L(:,ii)-ratio_data;
    Num=log(ratio_data_L(:,ii))-log(ratio_data);
    Denom=log(lambdatemp(ii))-log(lambda(ii));
    S_L(:,ii)=Num/Denom;
    %-------------Layer Thickess---------------
    ttemp=t;
    ttemp(ii)=t(ii)*1.01;
    [deltaR_data_t(:,ii),ratio_data_t(:,ii)]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,ttemp,eta,r_pump,r_probe,A_pump,nnodes);
    delta_t(:,ii)=ratio_data_t(:,ii)-ratio_data;
    Num=log(ratio_data_t(:,ii))-log(ratio_data);
    Denom=log(ttemp(ii))-log(t(ii));
    S_t(:,ii)=Num/Denom;
    
    %--------------------------------------------
    %-------------Anisotropy---------------
 
    etatemp(ii) = eta(ii)*1.01;
    %etatemp=eta;
    %etatemp(ii)=etatemp(ii)*1.01;
    [deltaR_data_eta(:,ii),ratio_data_eta(:,ii)]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,etatemp,r_pump,r_probe,A_pump,nnodes);
    delta_eta(:,ii)=ratio_data_eta(:,ii)-ratio_data;
    Num=log(ratio_data_eta(:,ii))-log(ratio_data);
    Denom=log(etatemp(ii))-log(eta(ii));
    S_eta(:,ii)=Num/Denom;
    %--------------------------------------------
    etatemp(ii) = eta(ii)*1.01;
    %etatemp=eta;
    %etatemp(ii)=etatemp(ii)*1.01;
    [deltaR_data_eta(:,ii),ratio_data_eta(:,ii)]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,etatemp,r_pump,r_probe,A_pump,nnodes);
    delta_eta(:,ii)=ratio_data_eta(:,ii)-ratio_data;
    Num=log(ratio_data_eta(:,ii))-log(ratio_data);
    Denom=log(etatemp(ii))-log(eta(ii));
    S_eta(:,ii)=Num/Denom;
    
end
% ftemp=f*1.01;
% [deltaR_data_f,ratio_data_f]=TDTR_REFL_VV3(tdelay,TCR,tau_rep,ftemp,lambda,C,t,eta,r_pump,r_probe,A_pump);
% delta_f(:,1)=ratio_data_f-ratio_data;
% Num=log(ratio_data_f)-log(ratio_data);
% Denom=log(ftemp)-log(f);
% S_f=Num/Denom;

r_pumptemp=r_pump*1.01;
[deltaR_data_r_pump,ratio_data_r_pump]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pumptemp,r_probe,A_pump,nnodes);
delta_r_pump(:,1)=ratio_data_r_pump-ratio_data;
Num=log(ratio_data_r_pump)-log(ratio_data);
Denom=log(r_pumptemp)-log(r_pump);
S_r_pump=Num/Denom;

r_probetemp=r_probe*1.01;
[deltaR_data_r_probe,ratio_data_r_probe]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probetemp,A_pump,nnodes);
delta_r_probe(:,1)=ratio_data_r_probe-ratio_data;
Num=log(ratio_data_r_probe)-log(ratio_data);
Denom=log(r_probetemp)-log(r_probe);
S_r_probe=Num/Denom;



close all