%% =========================================================
%% TDTR AVERAGING FOR FITTING
%% - keeps original time axis
%% - fitting compatible
%% - one output txt only
%% =========================================================

clear; clc;

%% ===== FILE SELECT =====
start_path = fullfile(getenv('HOME'),'Documents','TDTR');

[filenames, pathname] = uigetfile('*.*',...
    'Select TDTR files',...
    start_path,...
    'MultiSelect','on');

if isequal(filenames,0)
    return
end

if ischar(filenames)
    filenames = {filenames};
end

nfiles = length(filenames);

%% =========================================================
%% FIRST FILE = REFERENCE AXIS
%% =========================================================

ref_data = readmatrix(...
    fullfile(pathname,filenames{1}),...
    'FileType','text');

time_ref = ref_data(:,2);   % ps axis 그대로 사용

N = length(time_ref);

%% ===== storage =====
Vin_all  = zeros(N,nfiles);
Vout_all = zeros(N,nfiles);

%% =========================================================
%% LOOP
%% =========================================================

for i = 1:nfiles

    fprintf('Processing %d / %d\n',i,nfiles);

    %% ===== load =====
    data = readmatrix(...
        fullfile(pathname,filenames{i}),...
        'FileType','text');

    time_exp = data(:,2)*1e-12;

    Vin  = data(:,3);
    Vout = data(:,4);

    %% =====================================================
    %% PHASE CORRECTION
    %% =====================================================

    Theta = fminsearch(@(X) phase_jump(Vout,Vin,time_exp,X),0);

    Vin_n  = Vin*cos(Theta) - Vout*sin(Theta);
    Vout_n = Vout*cos(Theta) + Vin*sin(Theta);

    Vin  = Vin_n;
    Vout = Vout_n;

    %% =====================================================
    %% TIME ZERO CORRECTION
    %% =====================================================

    Vin_max = max(Vin);

    half_val = Vin_max/2;

    [~,idx_peak] = max(Vin);

    idx_half = find(Vin(1:idx_peak) >= half_val,1);

    t_zero = time_exp(idx_half);

    time_exp = time_exp - t_zero;

    %% ===== convert to ps =====
    td_ps = time_exp*1e12;

    %% =====================================================
    %% OFFSET CORRECTION (BASELINE REMOVAL)
    %% =====================================================
    % 
    % % negative delay region (baseline)
    % neg_idx = time_exp < min(time_exp)*0.2;
    % 
    % % baseline
    % Vin_offset  = mean(Vin(neg_idx));
    % 
    % 
    % % correction
    % Vin  = Vin - Vin_offset;


    %% =====================================================
    %% INTERPOLATE TO REFERENCE AXIS
    %% =====================================================

    Vin_interp = interp1(...
        td_ps,...
        Vin,...
        time_ref,...
        'linear',...
        'extrap');

    Vout_interp = interp1(...
        td_ps,...
        Vout,...
        time_ref,...
        'linear',...
        'extrap');

    %% ===== store =====
    Vin_all(:,i)  = Vin_interp;
    Vout_all(:,i) = Vout_interp;

end

%% =========================================================
%% AVERAGE
%% =========================================================

Vin_avg  = mean(Vin_all,2);
Vout_avg = mean(Vout_all,2);

Vratio_avg = -Vin_avg ./ Vout_avg;

%% =========================================================
%% PLOT
%% =========================================================

figure(1); clf;

plot(time_ref,Vin_avg,...
    'k','LineWidth',3);

xlabel('Delay time (ps)');
ylabel('Vin (uV)');

grid on;

figure(2); clf;

plot(time_ref,Vout_avg,...
    'k','LineWidth',3);

xlabel('Delay time (ps)');
ylabel('Vout (uV)');

grid on;

figure(3); clf;

idx_positive = time_ref > 0;

semilogx(time_ref(idx_positive),...
         Vratio_avg(idx_positive),...
         'k','LineWidth',3);

xlim([100 max(time_ref)])

xlabel('Delay time (ps)');
ylabel('-Vin/Vout');

grid on;

%% =========================================================
%% EXPORT
%% SAME FORMAT AS RAW FILE
%% =========================================================

dummy = zeros(size(time_ref));

output = [dummy,...
          time_ref,...
          Vin_avg,...
          Vout_avg];

writematrix(output,...
    fullfile(pathname,'TDTR_Average.txt'),...
    'Delimiter','tab');

fprintf('\n');
fprintf('TDTR_Average.txt exported.\n');