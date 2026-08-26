clear
clc
%------------------------------------------------------------------%
%%%-------------------------------------------------------------%%%
% Sample 24-4 ( 600 C deposited - O2 ratio 0.5% clean as-deposited%
%-----Al peak around 26 ps which is 80 nm Al---------------%
%-----Al deposited in sputter-II (5 mT, 50W, 140s)-Glass slide available%
%------------------------------------------------------------------%
%%%--------------------------------------------------------------%%%
T = [30 40 50 60 70 80 90 100 120];
t=T+273.15;
k_sap = 10.^(-1.047*log10(t)+2.1185)*100*1.126;
c_al = -0.536+1.1967*log10(t);
c_sap = -9.0786+4.922*log10(t);
c_si= 0.002*t+1.0;
k_si= 277.65-0.461*t;
c_mgo= -6.975+0.036*t++3.129e-4*t.^2-2.13873e-6*t.^3+4.60296e-9*t.^4-3.39944e-12*t.^5;
k_mgo= 100*(1.01-0.00165*t);
c_mgal2o4=15.225-0.73866*t+0.00747*t.^2-3.06852e-5*t.^3+5.72753e-8*t.^4-4.0378e-11*t.^5;
k_mgal2o4=100*(0.91346-0.00134*t);

GG1=0.2;

for n1=1:length(t)
    tt = t(n1);
    x = linspace(1e-3,750/tt,100);
    c_vo2_debye(n1) = 9*9.825e28*1.38e-23*(tt/750)^3*sum(x.^4.*exp(x)./(exp(x)-1).^2)*(x(3)-x(2))/1e6;
end

for n1=1:length(t)
    tt = t(n1);
    x = linspace(1e-3,425/tt,100);
    c_al_debye(n1) = 9*6.35e28*1.38e-23*(tt/425)^3*sum(x.^4.*exp(x)./(exp(x)-1).^2)*(x(3)-x(2))/1e6;
end
c_vo2_ryder(1:4)=[3.0883
3.14
3.29242
4.95524];

c_vo2_ryder(5:13)=[3.63487;3.56;
3.48973
3.45491
3.46657
3.47824
3.4899
3.50157
3.51323];


% linear

c_vo2=3+(T-25)*0.6/95;

for n1 = 1:length(T)
    [Xsol,tdelay_data,ratio_data,ratio_model] = TDTR_read(T(n1),c_sap(n1),c_al(n1),c_vo2_debye(n1),k_sap(n1));
     p1(n1) = Xsol(1); 
     p2(n1)=Xsol(2); 
     r_data(:,n1) = ratio_data; 
     r_model(:,n1) = ratio_model;
     k_eff1(n1)= 37/(1/GG1+37/p1(n1)+1/p2(n1));
     k_eff2(n1)= 37/(37/p1(n1)+1/p2(n1));
     clc
     perc = n1/length(T)*100
end
figure(120);plot(T,p1,'ko')
figure(130);plot(T,p2,'ko')
figure(160);plot(T,k_eff2,'ko')
% figure(140);plot(T,c_sap,'ko')
% figure(150);plot(T,k_sap,'ko')


CorrectedExpDataFileName = sprintf('%s.txt','output');    
fid1 = fopen(CorrectedExpDataFileName,'w'); 
for j=1:length(T)
fprintf(fid1,'\n  %6f %6f %6f %6f %6f %6f',t(j)+3,GG1,p1(j),p2(j),k_eff1(j),k_eff2(j)); 
end
fclose(fid1);



% T = [70 75 80 90 100 110 120 130 140];
% for n1 = 1:length(T)
%     [Xsol,tdelay_data,ratio_data,ratio_model] = TDTR_read(T(n1),c_vo2(n1),c_sap(n1),c_al(n1),k_sap(n1));
%      k_vo2(n1) = Xsol(1); r_data(:,n1) = ratio_data; r_model(:,n1) = ratio_model;
%      clc
%      perc = n1/length(T)*100
% endap