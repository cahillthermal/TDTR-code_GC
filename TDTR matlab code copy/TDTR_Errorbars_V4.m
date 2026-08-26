fprintf('Calculating Errobar\n')    
fprintf('YErr(n,:) = uncertainty (absolute) in X due to uncertainty in parameter Y(n)\n')
Xsol=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes),Xguess)
%Initial Guess, Variable to fit experiment to ... must change "TDTR_FIT_V4"program also
    for ii=1:length(lambda)
        %-------Specific Heat--------------
        if C_consider(ii)==1
            Ctemp=C;
            Ctemp(ii)=C(ii)*(1+Cperc_err(ii));
            Xsoltemp=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,Ctemp,t,eta,r_pump,r_probe,A_pump,nnodes),Xguess);
            for n=1:length(Xguess)
                CErr(ii,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable C(ii)
            end
            CErr
        end
        
        %-------Thermal Conductivity--------------
        if L_consider(ii)==1
            lambdatemp=lambda;
            lambdatemp(ii)=lambda(ii)*(1+lambdaperc_err(ii));;
            Xsoltemp=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambdatemp,C,t,eta,r_pump,r_probe,A_pump,nnodes),Xguess);
            for n=1:length(Xguess)
                lambdaErr(ii,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable lambda(ii)
            end
            lambdaErr
        end
        
        %-------Layer Thickness--------------
        if t_consider(ii)==1
            ttemp=t;
            ttemp(ii)=t(ii)*(1+tperc_err(ii));
            Xsoltemp=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,ttemp,eta,r_pump,r_probe,A_pump,nnodes),Xguess);
            for n=1:length(Xguess)
                tErr(ii,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable t(ii)
            end
            tErr
        end
    end
%-------Probe Radius--------------
if r_probe_consider==1
r_probetemp=r_probe*(1+r_err);
Xsoltemp=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probetemp,A_pump,nnodes),Xguess);
for n=1:length(Xguess)
    r_probeErr(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_probe(ii)
end
r_probeErr
end
%-------Pump Radius--------------
if r_pump_consider==1
r_pumptemp=r_pump*(1+r_err);
Xsoltemp=fminsearch(@(X) TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pumptemp,r_probe,A_pump,nnodes),Xguess);
for n=1:length(Xguess)
    r_pumpErr(1,n)=abs(Xsoltemp(n)-Xguess(n)); %Error in X(n) due to variable r_pump(ii)
end
r_pumpErr
end
%--------Phase Error-------------
if phase_consider==1
    radphase=pi/180*degphase;
    Vtemp=(Vin_data+sqrt(-1)*Vout_data)*exp(sqrt(-1)*radphase);
    Vin_phaseshifted=real(Vtemp);
    Vout_phaseshifted=imag(Vtemp);
    ratio_phaseshifted=-Vin_phaseshifted./Vout_phaseshifted;
    Xsoltemp=fminsearch(@(X) TDTR_FIT_V4(X,ratio_phaseshifted,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes),Xguess);
for n=1:length(Xguess)
    phaseErr(1,n)=abs(Xsoltemp(n)-Xguess(n));
end
phaseErr
end

ErrSummary_perc=[CErr;lambdaErr;tErr;r_probeErr;r_pumpErr;phaseErr]./(ones(3*length(lambda)+3,1)*Xsol); %percent error broken by variable

fprintf('Errorbar calculated\n')
fprintf('Errorbar breakdown:\n')
fprintf('Percent Err from C:\n')
CErr
fprintf('Percent Err from lambda:\n')
lambdaErr
fprintf('Percent Err from h:\n')
tErr
fprintf('Percent Err from spot size:\n')
r_probeErr
r_pumpErr
fprintf('Percent Err from phase:\n')
phaseErr
%----------------------------------------------------
fprintf('Total Percent Error:\n')
kErr_perc=sqrt(sum(ErrSummary_perc.^2,1)) %total percent error in each fitted parameter
fprintf('Total Absolute Error:\n')
kErr_abs=kErr_perc.*Xsol %total absolute error in each fitted parameter