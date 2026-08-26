%% ============================================================
% TDTR_Independent_Analysis.m
%
% Independent TDTR analysis script
%
% TDTR_MAIN_V4.m and TDTR_FIT_V4.m are NOT modified.
%
% ------------------------------------------------------------
% USER SETTINGS ARE ONLY IN SECTION 1.
% ------------------------------------------------------------
%
% Example 1:
%   fit_params  = {'k2','k4'};
%   sens_params = {'k2','k4'};
%   scan_param  = 'k2';
%   refit_param = 'k4';
%
% Example 2:
%   fit_params  = {'k3','k5'};
%   sens_params = {'k3','k5'};
%   scan_param  = 'k3';
%   refit_param = 'k5';
%
% Example 3:
%   fit_params  = {'k2','t1'};
%   sens_params = {'k2','t1'};
%   scan_param  = 'k2';
%   refit_param = 't1';
%
% Sensitivity and fitting parameters are completely independent.
%
% ============================================================

clear;
clc;
close all;

tic


%% ============================================================
% 1. USER SETTINGS
% ============================================================

% ------------------------------------------------------------
% SAMPLE PARAMETERS
% ------------------------------------------------------------

lambda = [153 0.5 153 0.2 35];
% W/m-K

C = [2.42 3.1 2.42 0.1 3.1]*1e6;
% J/m^3-K

t = [29.1 2 20.9 1 4.3e5]*1e-9;
% m


% ------------------------------------------------------------
% FITTING PARAMETERS
%
% These are the parameters actually fitted to experimental data.
%
% Examples:
%   {'k2','k4'}
%   {'k3','k5'}
%   {'k2','t1'}
%   {'k3','C4'}
% ------------------------------------------------------------

fit_params = {'k2','k4'};


% ------------------------------------------------------------
% SENSITIVITY PARAMETERS
%
% Completely independent from fit_params.
%
% Examples:
%   {'k2','k4'}
%   {'k3','k5'}
%   {'k1','k2','k3','k4','k5'}
%   {'t1','t2','t3','t4','t5'}
%   {'k2','t1'}
% ------------------------------------------------------------

sens_params = {'k2','k4'};


% ------------------------------------------------------------
% SCAN PARAMETERS
%
% scan_param:
%   parameter that is fixed and scanned
%
% refit_param:
%   parameter that is refitted at each scan point
%
% Example:
%   scan_param  = 'k3';
%   refit_param = 'k5';
%
% IMPORTANT:
% Both parameters must be included in fit_params.
% ------------------------------------------------------------

scan_param  = 'k2';
refit_param = 'k4';


% ------------------------------------------------------------
% SCAN RANGE
%
% +/- 30% around the original fitted value
% ------------------------------------------------------------

scan_fraction = 0.30;

Nscan = 15;


%% ============================================================
% 2. TDTR PARAMETERS
% ============================================================

eta = ones(1,numel(lambda));

r = 9.7e-6;                    % m

f = 9.3e6;                     % Hz

r_pump = r;
r_probe = r;

tau_rep = 1/74.88e6;           % s

A_pump = 12e-3;                % W

TCR = 1e-4;

nnodes = 35;


%% ============================================================
% 3. TIME DELAY RANGE
% ============================================================

tdelay_min = 100e-12;
tdelay_max = 3600e-12;


% IMPORTANT:
% This is the same logarithmic time-delay spacing used
% in the original TDTR_MAIN_V4 sensitivity calculation.

tdelay_sens = logspace( ...
    log10(tdelay_min), ...
    log10(tdelay_max), ...
    60)';


%% ============================================================
% 4. CHECK USER SETTINGS
% ============================================================

if length(fit_params) ~= 2
    error('fit_params must contain exactly TWO parameters.');
end


if ~ismember(scan_param,fit_params)

    error(['scan_param must be one of fit_params. ', ...
           'Current scan_param = ',scan_param]);

end


if ~ismember(refit_param,fit_params)

    error(['refit_param must be one of fit_params. ', ...
           'Current refit_param = ',refit_param]);

end


if strcmp(scan_param,refit_param)

    error('scan_param and refit_param must be different.');

end


%% ============================================================
% 5. TDTR BASELINE FOR SENSITIVITY
% ============================================================

[~,ratio_sens_baseline] = TDTR_REFL_V4( ...
    tdelay_sens, ...
    TCR, ...
    tau_rep, ...
    f, ...
    lambda, ...
    C, ...
    t, ...
    eta, ...
    r_pump, ...
    r_probe, ...
    A_pump, ...
    nnodes);


%% ============================================================
% 6. CALCULATE SENSITIVITY
% ============================================================

S = zeros( ...
    length(tdelay_sens), ...
    length(sens_params));


for n = 1:length(sens_params)

    param = sens_params{n};

    type = param(1);
    layer = str2double(param(2:end));


    lambda_temp = lambda;
    C_temp = C;
    t_temp = t;
    eta_temp = eta;


    switch type

        case 'k'

            lambda_temp(layer) = ...
                lambda(layer)*1.01;

            % Keep anisotropy relation consistent
            eta_temp = ...
                eta .* lambda ./ lambda_temp;


            denominator = ...
                log(lambda_temp(layer)) ...
                - log(lambda(layer));


        case 't'

            t_temp(layer) = ...
                t(layer)*1.01;

            denominator = ...
                log(t_temp(layer)) ...
                - log(t(layer));


        case 'C'

            C_temp(layer) = ...
                C(layer)*1.01;

            denominator = ...
                log(C_temp(layer)) ...
                - log(C(layer));


        otherwise

            error('Unknown sensitivity parameter: %s',param);

    end


    [~,ratio_temp] = TDTR_REFL_V4( ...
        tdelay_sens, ...
        TCR, ...
        tau_rep, ...
        f, ...
        lambda_temp, ...
        C_temp, ...
        t_temp, ...
        eta_temp, ...
        r_pump, ...
        r_probe, ...
        A_pump, ...
        nnodes);


    S(:,n) = ...
        (log(ratio_temp) ...
        - log(ratio_sens_baseline)) ...
        ./ denominator;

end


%% ============================================================
% 7. SENSITIVITY PLOT
% ============================================================

figure;

for n = 1:length(sens_params)

    semilogx( ...
        tdelay_sens*1e12, ...
        S(:,n), ...
        'o', ...
        'MarkerSize',7, ...
        'LineWidth',1.0);

    hold on;

end


xlabel('td (ps)','FontSize',16);

ylabel('Ratio Sensitivity','FontSize',16);

legend( ...
    sens_params, ...
    'Location','best');

set(gca,'FontSize',16);

grid on;
box on;

xlim([100 3600]);

title('TDTR Sensitivity','FontSize',16);


%% ============================================================
% 8. IMPORT EXPERIMENTAL DATA
% ============================================================

fprintf('\nSelect TDTR experimental data file.\n');

[filename,pathname] = uigetfile( ...
    '*.*', ...
    'Select TDTR data file', ...
    '/Users/gaeunchoi/Documents/TDTR');


if isequal(filename,0)

    fprintf('File selection cancelled.\n');

    return;

end


FileNames = fullfile(pathname,filename);


fprintf('\nSelected file:\n%s\n\n',FileNames);


%% ============================================================
% 9. READ EXPERIMENTAL DATA
% ============================================================

[tdelay_raw,Vout_raw,Vin_raw,ratio_raw,...
    Del_Vout,Del_Phase] = ...
    GetExpData( ...
    FileNames, ...
    filename, ...
    pathname, ...
    0,0);


tdelay_raw = tdelay_raw*1e-12;


%% ============================================================
% 10. EXTRACT FITTING RANGE
% ============================================================

[tdelay_data,Vin_data] = ...
    extract_interior_V4( ...
    tdelay_raw, ...
    Vin_raw, ...
    tdelay_min, ...
    tdelay_max);


[tdelay_data,Vout_data] = ...
    extract_interior_V4( ...
    tdelay_raw, ...
    Vout_raw, ...
    tdelay_min, ...
    tdelay_max);


[tdelay_data,ratio_data] = ...
    extract_interior_V4( ...
    tdelay_raw, ...
    ratio_raw, ...
    tdelay_min, ...
    tdelay_max);


fprintf('Data imported successfully.\n');

fprintf('Number of data points = %d\n', ...
    length(tdelay_data));


%% ============================================================
% 11. ORIGINAL FIT
% ============================================================

Xguess = zeros(1,length(fit_params));


for n = 1:length(fit_params)

    Xguess(n) = get_parameter( ...
        fit_params{n}, ...
        lambda, ...
        C, ...
        t);

end


fprintf('\nInitial guesses:\n');

for n = 1:length(fit_params)

    fprintf( ...
        '%s = %.6g\n', ...
        fit_params{n}, ...
        Xguess(n));

end


%% ------------------------------------------------------------
% Perform two-parameter fit
% ------------------------------------------------------------

[Xsol,Z_original] = fminsearch( ...
    @(X) fit_function( ...
        X, ...
        fit_params, ...
        ratio_data, ...
        tdelay_data, ...
        TCR, ...
        tau_rep, ...
        f, ...
        lambda, ...
        C, ...
        t, ...
        eta, ...
        r_pump, ...
        r_probe, ...
        A_pump, ...
        nnodes), ...
    Xguess, ...
    optimset('TolX',1e-4));


%% ============================================================
% 12. ORIGINAL FIT RESULT
% ============================================================

fprintf('\n');

fprintf('============================================\n');

fprintf(' ORIGINAL TWO-PARAMETER FIT\n');

fprintf('============================================\n');


for n = 1:length(fit_params)

    fprintf( ...
        '%s = %.6g\n', ...
        fit_params{n}, ...
        Xsol(n));

end


fprintf('Z = %.6g\n',Z_original);

fprintf('============================================\n');


%% ============================================================
% 13. CALCULATE ORIGINAL FIT CURVE
% ============================================================

[~,ratio_fit] = fit_function( ...
    Xsol, ...
    fit_params, ...
    ratio_data, ...
    tdelay_data, ...
    TCR, ...
    tau_rep, ...
    f, ...
    lambda, ...
    C, ...
    t, ...
    eta, ...
    r_pump, ...
    r_probe, ...
    A_pump, ...
    nnodes);


%% ============================================================
% 14. ORIGINAL TDTR FIT PLOT
% ============================================================

figure;

semilogx( ...
    tdelay_data*1e12, ...
    ratio_data, ...
    'ob', ...
    'MarkerSize',5, ...
    'LineWidth',0.8);

hold on;

semilogx( ...
    tdelay_data*1e12, ...
    ratio_fit, ...
    'r-', ...
    'LineWidth',2);


xlabel('Time delay (ps)','FontSize',16);

ylabel('Vratio','FontSize',16);

legend( ...
    'Experiment', ...
    'Original fit', ...
    'Location','best');

set(gca,'FontSize',16);

grid on;
box on;

xlim([tdelay_min tdelay_max]*1e12);

title('Original TDTR Fit','FontSize',16);


%% ============================================================
% 15. GET ORIGINAL FIT VALUES FOR SCAN
% ============================================================

scan_index = find( ...
    strcmp(fit_params,scan_param));


refit_index = find( ...
    strcmp(fit_params,refit_param));


scan_fit = Xsol(scan_index);

refit_fit = Xsol(refit_index);


%% ============================================================
% 16. DEFINE SCAN RANGE
% ============================================================

scan_min = ...
    (1-scan_fraction)*scan_fit;

scan_max = ...
    (1+scan_fraction)*scan_fit;


scan_values = linspace( ...
    scan_min, ...
    scan_max, ...
    Nscan);


refit_values = zeros(size(scan_values));

Z_scan = zeros(size(scan_values));


%% ============================================================
% 17. SCAN ONE PARAMETER
%     AND REFIT THE OTHER
% ============================================================

fprintf('\n');

fprintf('============================================\n');

fprintf(' PARAMETER CORRELATION SCAN\n');

fprintf('============================================\n');

fprintf('Scan parameter : %s\n',scan_param);

fprintf('Refit parameter: %s\n',refit_param);

fprintf('Scan range     : %.6g to %.6g\n', ...
    scan_min,scan_max);

fprintf('============================================\n\n');


for ii = 1:length(scan_values)

    scan_value = scan_values(ii);

    % Start from original fitted value
    Xstart = refit_fit;

    % Objective function
    objective = @(X) fit_function( ...
        make_X(X,scan_value,scan_index,refit_index), ...
        fit_params, ...
        ratio_data, ...
        tdelay_data, ...
        TCR,tau_rep,f, ...
        lambda,C,t,eta, ...
        r_pump,r_probe,A_pump,nnodes);

    % Refit remaining parameter
    [refit_best,Z_best] = fminsearch( ...
        objective, ...
        Xstart, ...
        optimset('TolX',1e-4));

    refit_values(ii) = refit_best;
    Z_scan(ii) = Z_best;

    fprintf( ...
        '%2d/%2d: %s = %.6g --> %s = %.6g, Z = %.6g\n', ...
        ii, ...
        length(scan_values), ...
        scan_param, ...
        scan_value, ...
        refit_param, ...
        refit_best, ...
        Z_best);

end


%% ============================================================
% 18. PARAMETER CORRELATION PLOT
% ============================================================

figure;

plot( ...
    scan_values, ...
    refit_values, ...
    'o-', ...
    'LineWidth',1.5, ...
    'MarkerSize',6);

hold on;


% Original fitted point

plot( ...
    scan_fit, ...
    refit_fit, ...
    's', ...
    'MarkerSize',10, ...
    'LineWidth',1.5);


% Get axis limits AFTER plotting

xl = xlim;
yl = ylim;


% Horizontal dashed line from y-axis to original point

plot( ...
    [xl(1),scan_fit], ...
    [refit_fit,refit_fit], ...
    'r--', ...
    'LineWidth',1.2);


% Vertical dashed line from x-axis to original point

plot( ...
    [scan_fit,scan_fit], ...
    [yl(1),refit_fit], ...
    'r--', ...
    'LineWidth',1.2);


xlabel( ...
    [scan_param,' (fit value / input unit)'], ...
    'FontSize',16);

ylabel( ...
    [refit_param,' (fit value / input unit)'], ...
    'FontSize',16);


title( ...
    [refit_param,' vs. ',scan_param], ...
    'FontSize',16);


legend( ...
    'Refit scan', ...
    'Original fit', ...
    'Location','best');


grid on;
box on;

set(gca,'FontSize',14);


%% ============================================================
% 19. FIT ERROR vs SCAN PARAMETER
% ============================================================

figure;

plot( ...
    scan_values, ...
    Z_scan, ...
    'o-', ...
    'LineWidth',1.5, ...
    'MarkerSize',6);

hold on;


% Original fitted point

plot( ...
    scan_fit, ...
    Z_original, ...
    's', ...
    'MarkerSize',10, ...
    'LineWidth',1.5);


xl = xlim;
yl = ylim;


% Horizontal dashed line

plot( ...
    [xl(1),scan_fit], ...
    [Z_original,Z_original], ...
    'r--', ...
    'LineWidth',1.2);


% Vertical dashed line

plot( ...
    [scan_fit,scan_fit], ...
    [yl(1),Z_original], ...
    'r--', ...
    'LineWidth',1.2);


xlabel( ...
    [scan_param,' (fit value / input unit)'], ...
    'FontSize',16);

ylabel('Fit Error, Z','FontSize',16);


title( ...
    ['Z vs. ',scan_param], ...
    'FontSize',16);


legend( ...
    'Scan', ...
    'Original fit', ...
    'Location','best');


grid on;
box on;

set(gca,'FontSize',14);


%% ============================================================
% 20. SAVE RESULTS
% ============================================================

Results = table( ...
    scan_values(:), ...
    refit_values(:), ...
    Z_scan(:), ...
    'VariableNames', ...
    {'ScanParameter','RefitParameter','Z'});


save( ...
    'TDTR_Independent_Analysis_results.mat', ...
    'Results', ...
    'fit_params', ...
    'sens_params', ...
    'scan_param', ...
    'refit_param', ...
    'Xsol', ...
    'Z_original', ...
    'S', ...
    'tdelay_sens');


%% ============================================================
% 21. FINISH
% ============================================================

fprintf('\n');

fprintf('============================================\n');

fprintf(' ANALYSIS COMPLETED\n');

fprintf('============================================\n');

fprintf('\nOriginal fit:\n');

for n = 1:length(fit_params)

    fprintf( ...
        '%s = %.6g\n', ...
        fit_params{n}, ...
        Xsol(n));

end

fprintf('Z = %.6g\n',Z_original);

fprintf('\nSensitivity parameters:\n');

fprintf('%s\n',strjoin(sens_params,', '));

fprintf('\nScan parameter : %s\n',scan_param);

fprintf('Refit parameter: %s\n',refit_param);

fprintf('\nResults saved to:\n');

fprintf('TDTR_Independent_Analysis_results.mat\n');

fprintf('============================================\n');


toc


%% ============================================================
% LOCAL FUNCTION 1
% GET PARAMETER VALUE
% ============================================================

function value = get_parameter( ...
    param,lambda,C,t)


type = param(1);

layer = str2double(param(2:end));


switch type

    case 'k'

        value = lambda(layer);


    case 't'

        value = t(layer);


    case 'C'

        value = C(layer);


    otherwise

        error( ...
            'Unknown parameter: %s',param);

end

end


%% ============================================================
% LOCAL FUNCTION 2
% TDTR FIT FUNCTION
% ============================================================

function [Z,ratio_model] = fit_function( ...
    X, ...
    fit_params, ...
    ratio_data, ...
    tdelay, ...
    TCR,tau_rep,f, ...
    lambda,C,t,eta, ...
    r_pump,r_probe, ...
    A_pump,nnodes)


L = lambda;

CC = C;

tt = t;

EE = eta;


for n = 1:length(fit_params)

    param = fit_params{n};

    type = param(1);

    layer = str2double(param(2:end));


    switch type

        case 'k'

            L(layer) = X(n);


        case 't'

            tt(layer) = X(n);


        case 'C'

            CC(layer) = X(n);


        otherwise

            error( ...
                'Unknown fitting parameter: %s',param);

    end

end


% Keep eta consistent with conductivity changes

EE = eta .* lambda ./ L;


[~,ratio_model] = ...
    TDTR_REFL_V4( ...
        tdelay, ...
        TCR, ...
        tau_rep, ...
        f, ...
        L, ...
        CC, ...
        tt, ...
        EE, ...
        r_pump, ...
        r_probe, ...
        A_pump, ...
        nnodes);


Z = sum( ...
    (ratio_model-ratio_data).^2);

end


%% ============================================================
% LOCAL FUNCTION 3
% BUILD TWO-PARAMETER VECTOR DURING SCAN
% ============================================================

function X = make_X( ...
    refit_value, ...
    scan_value, ...
    scan_index, ...
    refit_index)


X = zeros(1,2);

X(scan_index) = scan_value;

X(refit_index) = refit_value;

end