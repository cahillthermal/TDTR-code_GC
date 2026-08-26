clear all;
clc;

tic

psc = 0; % TRUE if pump spot changes
frac = 0.3 ; % fractional change in pump spot size over delay time

% lambda=[160 0.15 1.32 140]; %W/m-K
% C=[2.42 0.1 1.62 1.6]*1e6; %J/m^3-K
% t=[70.41 1 311 1e9]*1e-9; %m

%-------------TYPE THERMAL SYSTEM PARAMTERS HERE--------------
%Anticipated system properties (initial guess for fitting, if you do fitting/errorbar estimation)


%quartz Al grpahene/SEI
% lambda=[1.4 0.1 182 0.1 6 0.1 5]
% C=[1.64 0.1 2.42 0.1 1.54 0.1 2]*1e6
% t=[1e6 1 80 1 200 1 1e3]*1e-9
% %test Al/HOPG
% lambda=[1500 150 0.06 5.4]; % W/m-K k_Al=237 Au estimated 153 K Al/G 300K from paper 45 MW/m2K
% C=[24.2 2.42 0.1 1.55]*1e6; % J/m^3-K C_V=2.93 C_SiO2=1.66 C_Al=2.42 C_Si=1.63
% t=[1 54 1 0.5*1e3]*1e-9; % m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)

T_VN = 27;
t1=T_VN+273.15;

%----Al heat capacity--------------------------
A	=0.59262;
B	=14.13493;
C	=-38.9046;
D	=49.38745;
E	=-0.00336;
F	=-22.6436;	
C_Al_gw = A+B*(t1/1000)^1+C*(t1/1000)^2+D*(t1/1000)^3+F*(t1/1000)^4+E/(t1/1000)^2;	

% thermal boundary conductance of Cr2O3/sapphire G = 100 MW/m^2-K
%4layer%
% lambda=[178 17.8  0.15 14]; %W/m-K
% C=     [26.5 2.65 0.1 3.92]*1e6; %J/m^3-K
% t=     [1.5 95 1 1e6]*1e-9

%3layer
% lambda=[123.5  0.2 140]; %W/m-K
% C=     [2.42 0.1 1.6]*1e6; %J/m^3-Kf
% t=     [51.4 1 1e6]*1e-9

lambda=[153 153 0.5 153 0.2 35]; %W/m-K
C=     [24.2 2.42 3.1 2.42 0.1 3.1]*1e6; %J/m^3-K
t=     [1 19.1 2 20.9 1 4.3e5]*1e-9



%pointprobe = input('Thermal conductivity measured by 4 point probe? \n(1=yes, 0=no)\n');
%if pointprobe == 1
%    lambda(1) = input('input measured Al thermal conductivity (W/m-K) \n');
%end

 %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
%lambda(1)=10*lambda(2); %top 1nm 10X the rest of layer (simulates absorbtion)
%C(1)=10*C(2); %top 1nm 10X rest of layer
eta=ones(1,numel(lambda)); %isotropic layers, eta=kx/ky;

r=9.7e-6; %1/e^2 radius spot, in m. (Option below if not the same size) 5X 9.7 / 10X 4.9
f=9.3e6; %laser Modulation frequency, Hz

r_pump=r; %pump 1/e^2 radius, m
r_probe=r; %probe 1/e^2 radius, m
tau_rep=1/74.88e6; %laser repetition period, s 74.8(tdtr2) or 80(tdtr1)
A_pump=12e-3; %laser power (Watts) . . . only used for amplitude est.
TCR=1e-4; %coefficient of thermal reflectance . . . only used for amplitude est.

nnodes = 35; %this is the number of nodes used for numerical integration. (DEFAULT=35)  

%choose time delays for Sensitivity plots
tdelay=logspace(log10(100e-12),log10(3600e-12),60)'; %vector of time delays (used to generate sensitivity plots) unit: second
% if this tdelay is used for calculating the errors, make sure the points
% are a lot, otherwise it will end up with large mistakes
%Choose range of time delays to fit, sec
tdelay_min=100e-12;
% tdelay_min=1050e-12;
% tdelay_min=1e-12;
tdelay_max=3600e-12;

%------------- CALCULATE STEADY STATE HEATING--------------
absorbance = 0.12*1.08*0.8; %fraction of incident light absorbed at surface    
% Optical power calibration is 1.08 times the reading of the model 835 power meter.
% Transmission coefficients of the objectives at 785 nm are 0.87, 0.90, 0.80, 0.70 for 2? 5? 10? and 20? respectively.
A_pump1=8e-3
A_probe=4e-3
A_tot_powermeter = A_pump1+2*A_probe; %incident light intensity, Watts
dT_SS = SS_Heating(lambda,C,t,eta,r_pump,absorbance,A_tot_powermeter) %max steady state temperature rise

PP_Heating=absorbance*A_pump1*2*tau_rep/(pi*r_pump^2)/(C(1)*t(1))
%----------------------------------PROGRAM OPTIONS BEGIN--------
%Generate Sensitivity Plots?
senseplot=1;
%Import Data? 0 for no, 1 for yes
importdata=1;
%If so, which variable(s) are you fitting for?
Xguess=[lambda(3),lambda(5)];
% Xguess=[lambda(3)];%, lambda(3)]; %initial guess for solution, could be for simulated data (if nothing imported) or real data
%Xguess=[lambda(3) lambda(2)];%, lambda(3)]; %initial guess for solution, could be for simulated data (if nothing imported) or real data
% Xguess=[lambda(4) lambda(3)];
%Calculate Errorbars? 0 for no, 1 for yes (takes longer)
ebar=0;
% Want to save results?
saveint=0;

%----------------------------------PROGRAM OPTIONS END--------

%---------------------ERRORBAR OPTIONS----------------------
if ebar==1
   % Init_ErrorBar;
   
     C_consider=ones(length(C),1);
    L_consider=ones(length(C),1);
    t_consider=ones(length(C),1);
    r_probe_consider=1;
    r_pump_consider=1;
    phase_consider=1;
    CErr=zeros(length(C),length(Xguess));
    lambdaErr=zeros(length(lambda),length(Xguess));
    tErr=zeros(length(lambda),length(Xguess));
    r_probeErr=zeros(1,length(Xguess));
    r_pumpErr=zeros(1,length(Xguess));
    
    %define parameters NOT to consider in error analysis (saves time)
    t_consider(length(t))=0; %last layer is semi-inf
    t_consider(3)=0; %thermal interface conductance at fixed t/vary lambda
    %t_consider(5)=0; %thermal interface conductance at fixed t/vary lambda
    
        
    C_consider(3)=0; %thermal interface layer has no capacitance
    %C_consider(5)=0; %thermal interface layer has no capacitance
    
       
    L_consider(4)=0; %solving for this %<-------------COMMENT-if solving for only lambda(4)!  
   % L_consider(3)=0;
  
    
        
    %define percent uncertainty in each layer/parameter
    Cperc_err=[0.02 0.02 0 0.02 0]; %percent uncertainty in specific heat %i.e. 0.02 -> 2% uncertainty
    lambdaperc_err=[0.05 0.05 0.1 0 0];% percent uncertainty in thermal conducitivyt
    tperc_err=[0.05 0.05 0 0.02 0];  % percent uncertainty in layer thickness
    r_err=0.05;  % percent uncertainty in beam radius
    radphase=0.006;  %phase error in degrees
    
    
end

%-----------------Make sensitivity plots--------------
[deltaR_data,ratio_data]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes);

if senseplot==1
    TDTR_SensPlots_V4
end

%--------------Import Data---------------------------
% Select Data File
if importdata==1
    [filename,pathname]=uigetfile('*.*',' ','/Users/gaeunchoi/Documents/TDTR');
    FileNames=strcat(pathname,filename);

% % Fix Data File
% if importdata==1
%     pathname = '/Users/gaeunchoi/Documents/TDTR/260121';
%     filename = 'pump12mW_pb6mW_92mVsig10X_1_1.txt';
%     FileNames = fullfile(pathname, filename);


    % [tdelay_raw,Vout_raw,Vin_raw,ratio_raw,Del_Vout,Del_Phase] = GetExpData(FileNames,filename,pathname,0);
    % Set the last number to 0 to hide the raw data plot / 1 to show it
    [tdelay_raw,Vout_raw,Vin_raw,ratio_raw,Del_Vout,Del_Phase] = GetExpData(FileNames,filename,pathname,0,0);


    tdelay_raw = tdelay_raw*1e-12;
    [DD,II] = min(abs(tdelay_raw-100e-12));
    Results_ratio = ratio_raw(II)
    
    [tdelay_data,Vin_data] = extract_interior_V4(tdelay_raw,Vin_raw,tdelay_min,tdelay_max);
    [tdelay_data,Vout_data] = extract_interior_V4(tdelay_raw,Vout_raw,tdelay_min,tdelay_max);
    [tdelay_data,ratio_data] = extract_interior_V4(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
    
if psc, r_pump = r_pump*(1 + frac*tdelay_data/tdelay_data(end)); end
    
%--------------Perform Fit (skips if no data imported)--------------------------
lambda_input = lambda;
C_input = C;
t_input = t;

    Xsol=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay_data,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes,0),Xguess,optimset('TolX',1e-4))
    Xguess=Xsol;
    tdelay=tdelay_data;
    fprintf('Data fit completed\n')
    TDTR_FIT_V4(Xsol, ratio_data, tdelay_data, ...
    TCR, tau_rep, f, lambda, C, t, eta, r_pump, r_probe, A_pump, nnodes, 1);

    %  Display the original input values using annotation
    % lambda_input, C_input, t_input: n-layer arrays
    nLayer = length(lambda_input);
    
    % Split the text into separate lines in a cell array
    txt_lines = cell(3,1);
    
    % lambda
    lambda_str = sprintf('%g ', lambda_input);
    txt_lines{1} = ['lambda = [', lambda_str, ']'];
    % C
    C_str = sprintf('%g ', C_input/1e6);
    txt_lines{2} = ['C = [', C_str, '] x10^6'];
    % t
    t_str = sprintf('%g ', t_input*1e9);
    txt_lines{3} = ['t = [', t_str, '] nm'];
    
    annotation('textbox', [0.67 0.75 0.3 0.15], 'String', txt_lines, ...
               'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k', ...
               'Interpreter', 'none');  % Use 'none' to disable special text interpretation


else
    [tdelay,Vin_data] = extract_interior_V4(tdelay,real(deltaR_data),tdelay_min,tdelay_max);
    [tdelay,Vout_data] = extract_interior_V4(tdelay,imag(deltaR_data),tdelay_min,tdelay_max);
    [tdelay,ratio_data] = extract_interior_V4(tdelay,ratio_data,tdelay_min,tdelay_max);
    Xsol = Xguess;
end

%--------------Compute Errorbars---------------------
if ebar==1
    TDTR_Errorbars_V4
end
fprintf('Fitting Solution:\n')

% 
% Results_K = Xsol(1);
% Results_G = Xsol(2);
% lambda(4) = Xsol(1);
% lambda(3) = Xsol(2);



if importdata==1

[Z,ratio_model]=TDTR_FIT_V4(Xsol,ratio_data,tdelay_data,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes);

else

[deltaR_data,ratio_model]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes);

end


CalcRatioFileName = sprintf('%s%s.txt','ModelRatio_',filename);
fullCalcRatioName=sprintf('%s%s',pathname,CalcRatioFileName);
tdelay_save=tdelay*1e12;
CalcRatio_data=[tdelay_save Vin_data Vout_data ratio_data ratio_model];
% save(fullCalcRatioName,'CalcRatio_data','-ascii'); 
XsolName = sprintf('%s%s.txt','Xsol',filename);
fullXsolName=sprintf('%s%s',pathname,XsolName);
% save(fullXsolName,'Xsol','lambda','C','t','-ascii');

%--------------resistance calculation -----------
r=t(2)*1e9/Xsol(1)


% CalcRatioFileName = sprintf('%s%s.txt','ModelRatio_',FileNames);
% dlmwrite(CalcRatioFileName,[tdelay*1e12,ratio_model],'delimiter','\t')

%----------------------------------------------------
if saveint == 1

    figure; 
    semilogx(tdelay_data, ratio_model, 'r-', 'LineWidth', 2);  % Model
    hold on;
    plot(tdelay_data, ratio_data, 'ok');                      % Experiment Data

    xlabel('time delay (s)', 'FontSize', 18);
    ylabel('Vratio', 'FontSize', 18);
    legend('Model', 'Experiment');
    legend boxoff;
    set(gca, 'FontSize', 18);

    % Autoset y-axis range
    axis([min(tdelay_data), max(tdelay_data), 0, 1.2*max([ratio_data; ratio_model])]);

    % Remove the print command here if you prefer to save figures manually
    % print(fignum,'-dpng',strcat(pathname,filename,'_FIT.png'))

end
%----------------------------------------------------
fprintf('Program Completed\n')
beep;


%% linear scann
% [filename,pathname]=uigetfile('*.*',' ','C:\Users\qiyez\Google Drive\data transfer\Exp Data\TDTR data\2017');
% FileNames=strcat(pathname,filename);
% [tdelay_raw,Vout_raw,Vin_raw,ratio_raw,Del_Vout,Del_Phase] = GetExpData(FileNames,filename,pathname,0);
% tdelay_raw = tdelay_raw*1e-12;
% [DD,II] = min(abs(tdelay_raw-100e-12));
% Results_ratio = ratio_raw(II)
% tdelay_min = 0;
% tdelay_max = 80e-12;
% [tdelay_data,ratio_data] = extract_interior_V4(tdelay_raw,ratio_raw,tdelay_min,tdelay_max);
% 
% figure(9)
% subplot(1,2,1)
% plot(1e12*tdelay_data,ratio_data,'ob');
% xlabel('t (ps)')
% ylabel('ratio')
% grid on
% title ('80 nm Al/ PMMA')
% subplot(1,2,2)
% plot(1e12*tdelay_data,ratio_data,'ob');
% xlabel('t (ps)')
% ylabel('ratio')
% grid on
% title ('93 nm Al/ PS')


