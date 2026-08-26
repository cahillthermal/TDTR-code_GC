clear; clc;

%% ==========================================
%% Load averaged TDTR file
%% ==========================================

[file,path] = uigetfile('*.txt','Select TDTR_Average.txt');

if isequal(file,0)
    return
end

data = readmatrix(fullfile(path,file));

%% ==========================================
%% Columns
%% ==========================================

t    = data(:,2);   % ps
Vin  = data(:,3);
Vout = data(:,4);

%% ==========================================
%% Select signal
%% ==========================================

y = Vin;

% y = Vout;
% y = -Vin./Vout;

%% ==========================================
%% Fit region
%% ==========================================

fit_idx = (t>=35) & (t<=55);

t_fit = t(fit_idx);
y_fit = y(fit_idx);

%% ==========================================
%% Linear fit
%% ==========================================

p1 = polyfit(t_fit,y_fit,1);
baseline1 = polyval(p1,t);

%% ==========================================
%% Quadratic fit
%% ==========================================

p2 = polyfit(t_fit,y_fit,2);
baseline2 = polyval(p2,t);

%% ==========================================
%% Residuals
%% ==========================================

residual1 = y - baseline1;
residual2 = y - baseline2;

%% ==========================================
%% Figure 1
%% Original + Fits
%% ==========================================

figure(1); clf

plot(t_fit,y_fit,...
    'k','LineWidth',2)

hold on

plot(t_fit,baseline1(fit_idx),...
    '--','LineWidth',2)

plot(t_fit,baseline2(fit_idx),...
    '--','LineWidth',2)

xlim([35 55])

xlabel('Delay Time (ps)')
ylabel('Vratio')

title('Original Signal + Baseline Fits')

legend('Data',...
       'Linear Fit',...
       'Quadratic Fit',...
       'Location','best')

grid on

%% ==========================================
%% Figure 2
%% Linear residual
%% ==========================================

figure(2); clf

plot(t_fit,residual1(fit_idx),...
    'LineWidth',2)

hold on

xlim([35 55])

xlabel('Delay Time (ps)')
ylabel('Residual')

title('Linear Baseline Removed')

grid on

%% ==========================================
%% Figure 3
%% Quadratic residual
%% ==========================================

figure(3); clf

plot(t_fit,residual2(fit_idx),...
    'LineWidth',2)

hold on

xlim([35 55])

xlabel('Delay Time (ps)')
ylabel('Residual')

title('Quadratic Baseline Removed')

grid on

%% ==========================================
%% Figure 4
%% Residual comparison
%% ==========================================

figure(4); clf

plot(t_fit,residual1(fit_idx),...
    'LineWidth',2)

hold on

plot(t_fit,residual2(fit_idx),...
    'LineWidth',2)

xlim([35 55])

xlabel('Delay Time (ps)')
ylabel('Residual')

title('Residual Comparison')

legend('Linear Removed',...
       'Quadratic Removed',...
       'Location','best')

grid on