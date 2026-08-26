
function dT_SS = SS_Heating(lambda,C,h,eta,r,absorbance,A_tot_powermeter)

%Calculats the steady state heating of an anisotropic multilayer stack
%lambda:  thermal conductivity
%C: heat capacity
%h: film thicknesses
%r: laser 1/e2 radius
%abosorbance: fraction of laser power absorbed by the top layer
%A_tot_powermeter:  %total (pump + probe) power meter reading (assumes chopper is OFF!  If its not actually off, then you need to multiply the probe power by 2!)

%If you are running the function standalone, you can type the input
%parameters inside the "if" statement.  

%If you call the function from another program and
%provide inputs, it will skip over this is statement and use those inputs
%instead.

n_pt_783=2.8042;
k_pt_783= 4.8904;
Refl_Pt=((n_pt_783-1)^2+k_pt_783^2)/((n_pt_783+1)^2+k_pt_783^2)

n_VN_783=1.49;
k_VN_783= 3.44;
Refl_VN=((n_VN_783-1)^2+k_VN_783^2)/((n_VN_783+1)^2+k_VN_783^2)

if nargin==0
    absorbance =1-0.86; 

%-----------------------------
T_VN = 250;

t1=T_VN+273.15;

t11=t1/1000;
c_pt = 2.53603+1.17864*t11-0.82777*t11.^2+0.41152*t11.^3-0.05168*t11.^4
c_al = -0.536+1.1967*log10(t1);
c_VN = (43.93869 + 11.21546*t11 + -0.655382*t11.^2 + 0.044268*t11.^3 + -0.815717/t11.^2)*6.13/64.9482
y=0;%VN_1-y
c_VN2=exp(-0.428*y)*(10.367+26.46e-4*t1-1847e2/t1^2)*4.1833*6.13/64.9482;%"Thermodynamic properties of vanadium mononitride in its homogeneity range"


c_vo2=3+(T_VN-25)*0.6/95

A	=0.60877;
B	=13.99686;
C	=-38.55308;
D	=49.06315;
E	=-0.00341;
F	=-22.56093;
C_Al_gw = A+B*(t1/1000)^1+C*(t1/1000)^2+D*(t1/1000)^3+F*(t1/1000)^4+E/(t1/1000)^2



c_al = -0.536+1.1967*log10(t1);
c_sap= -9.0786+4.922*log10(t1); %very close to data in group web page 20-120C

k_sap_old = 10.^(-1.047*log10(t1)+2.1185)*100*1.126;

y0=9.01477;
A1=535.53105;
t111=47.54886;
A2=1657.24176;
t2=47.54922;
A3=69.32237;
t3=262.54622;
k_sap_gw= A1*exp(-t1/t111) + A2*exp(-t1/t2) + A3*exp(-t1/t3) + y0	

c_MgO=1.09829+10.13606*t11-8.72295*t11.^2-0.85468*t11.^3+2.97451*t11.^4
k_MgO=210.42151-1011.05592*t11+2193.98132*t11.^2-2235.47822*t11.^3+867.12576*t11.^4
%-----------------------------------------------------------
% SiC
lambda=[1800 180 0.15 380 0.1 220]; % W/m-K k_Al=237 Au estimated 153 K Al/G 300K from paper 45 MW/m2K   Al~61nm k~175 !!!!!L_Al matter when using small spot size++
C=[C_Al_gw*10 C_Al_gw 0.1 2.12 0.1 1.63]*1e6; % J/m^3-K C_V=2.93 C_SiO2=1.64 C_Al=2.42 C_Si=1.63
h=[1 83 1 5e6 1 1e6]*1e-9; % m ("A" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)  
 
%Pt-SiO2-MgO
% lambda=[170 17 0.2 1.3 0.2 146]; % W/m-K k_Al=237 Au estimated 153 K Al/G 300K from paper 45 MW/m2K   Al~61nm k~175
% C=[22.7 2.27 0.1 1.64 0.1 1.64]*1e6; % J/m^3-K C_V=2.93 C_SiO2=1.64 C_Al=2.42 C_Si=1.63
% h=[1 70 1 515 1 1e6]*1e-9; % m ("A" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)  

% VN
% lambda=[20.3 0.8 k_MgO]; % W/m-K k_Al=237 Au estimated 153 K Al/G 300K from paper 45 MW/m2K   Al~61nm k~175
% C=[c_VN2 0.1 c_MgO]*1e6; % J/m^3-K C_V=2.93 C_SiO2=1.64 C_Al=2.42 C_Si=1.63
% h=[32.6 1 5e6]*1e-9; % m ("A" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)  


% lambda=[200 20 1 k_MgO]; % W/m-K k_Al=237 Au estimated 153 K Al/G 300K from paper 45 MW/m2K   Al~61nm k~175
% C=[c_VN2*10 c_VN2 0.1 c_MgO]*1e6; % J/m^3-K C_V=2.93 C_SiO2=1.66 C_Al=2.42 C_Si=1.63
% h=[19 200 1 0.5e6]*1e-9; % m ("A" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
    %eta=ones(1,numel(lambda)); %isotropic layers, eta=kx/ky;
    eta=[1 1 1 1 1 1]
    r=10*1e-6;
    r_pump=r; %pump 1/e^2 radius, m
    r_probe=r; %probe 1/e^2 radius, m
    A_tot_powermeter = (20+2*10)*1e-3; %total (pump + probe) power meter reading (assumes chopper is OFF!  If its not actually off, then you need to multiply the probe power by 2!)
end

%-----------DO NOT MODIFY BELOW HERE------------------------------
f=0; %laser Modulation frequency, Hz
A_abs = absorbance*A_tot_powermeter %absorbed power
r_pump=r; %pump 1/e^2 radius, m
r_probe=r; %probe 1/e^2 radius, m
lambda_sub=lambda(length(lambda))*eta(length(eta));
kmin=1/(10000*r); %smallest wavevector (can't use exactly zero, because of a numerical issue -> 0/0 = NaN 
kmax=1/sqrt(r_pump^2+r_probe^2)*1.5; %maximum wavevector
dT_SS  =rombint_VV3(@(kvect) TDTR_TEMP_V4(kvect,f,lambda,C,h,eta,r_pump,r_probe,A_abs),kmin,kmax,1)

%-----------------------------------------------------------------

