%The main program tries to minimize "Z" by optimizing the variable(s) X
%This program Lets you:
%   1) Define the vector X: example, if lambda(3) is what you want to solve for,
%   then set X=lambda(3)....if you whish to simulatenous solve for more than one
%   variable, you can just define multiple variables (eg. X(1)=lambda(3), X(2)=lambda(4))
%
%   2) Define the fit.  Typically, this is the sum of the squares of the
%   residuals, but you might want to weight these by the sensitivity,
%   particularly if you don't intend to calculate the errorbars!

% function [Z,ratio_model]=TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes)
function [Z,ratio_model]=TDTR_FIT_V4(X,ratio_data,tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes,plotFlag)
if nargin < 15
    plotFlag = 0;
end



%Define the variables to be fit



% %fit absorption layer thickness by using the 

% VN
% h=X(1);
% C(1)=C(2)*h;
% t(2)=(194-h)*1e-9;
% lambda(1)=h*lambda(2);
% lambda(2)=X(2);
% % % 

% % t(2)=X(1);
%  lambda(2)=X(1);
%  lambda(1)=X(1)*17;
%  lambda(3)=X(2);
% fitting the FIB damaged Si
%  t(4)=X(1)/1e9;


lambda(2)=X(1);
lambda(4)=X(2);
% t(3)=X(2);



% 
% 
%fit thin FeRh C
% C(1)=X(1)*1e6;
% lambda(2)=X(2);


%--for graphite inplane fitting
% eta(3)=900/X(1);
%-------------------------------




  









[deltaR_model,ratio_model]=TDTR_REFL_V4(tdelay,TCR,tau_rep,f,lambda,C,t,eta,r_pump,r_probe,A_pump,nnodes);
%Uncomment the next three lines to see the non-linear optimization in
%action!

% figure (11);
% semilogx(1e12*tdelay, ratio_data, 'ob', 1e12*tdelay, ratio_model, 'g', 'LineWidth', 1.5);
% xlabel('t (ps)');
% ylabel('ratio');
% grid on;
% title('TDTR Fit Result');

if plotFlag
    figure; % MATLAB이 자동으로 번호 붙임
    semilogx(1e12*tdelay,ratio_data,'ob',1e12*tdelay,ratio_model,'r-');
    xlabel('t (ps)');
    ylabel('Vratio');
    grid on;
    title('TDTR Fit Result');

end



pause(0.001)

X
% X(1)/X(2)
res=(ratio_model-ratio_data).^2;
Z=sum(res)

format short g