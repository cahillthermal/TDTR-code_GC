%% ===== 설정 =====
tdelay_min = 100e-12;
lw = 2.5;

start_path = fullfile(getenv('HOME'),'Documents','TDTR');

[filenames, pathname] = uigetfile('*.*',...
    'Select TDTR data files',...
    start_path,...
    'MultiSelect','on');

if isequal(filenames,0)
    disp('File selection cancelled');
    return
end

if ischar(filenames)
    filenames = {filenames};
end

%% ===== figure 준비 =====
figure(1); clf; hold on; grid on;
xlabel('Delay time (ps)'); ylabel('Vin (mV)');

figure(2); clf; hold on; grid on;
xlabel('Delay time (ps)'); ylabel('Vout (mV)');

figure(3); clf; hold on; grid on;
set(gca,'XScale','log');
xlabel('Delay time (ps)');
ylabel('Vratio');

colors = lines(length(filenames));

%% ===== 파일 루프 =====
for i = 1:length(filenames)

    fullpath = fullfile(pathname, filenames{i});
    data = readmatrix(fullpath);

    fprintf('Processing file %d of %d: %s\n',i,length(filenames),filenames{i});

    %% ===== raw data =====
    time_exp = data(:,2)*1e-12;   % seconds

    Vin  = data(:,3);
    Vout = data(:,4);

    %% ===== zero index (phase_jump용) =====
    [~,II] = min(abs(time_exp));

    %% ===============================
    %% ===== VOUT AUTOCORRECTION =====
    %% ===============================

    Theta = fminsearch(@(X) phase_jump(Vout,Vin,time_exp,X),0);

    Vin_n  = Vin*cos(Theta) - Vout*sin(Theta);
    Vout_n = Vout*cos(Theta) + Vin*sin(Theta);

    Vin  = Vin_n;
    Vout = Vout_n;

    % outlier 제거
    if II>10
        STD  = std(Vout(1:(II-4)));
        MEAN = mean(Vout(1:min(2*II,length(Vout))));
    else
        STD  = std(Vout);
        MEAN = mean(Vout);
    end

    for k = 1:min(2*II,length(Vout))
        if abs(Vout(k)-MEAN) > 2.5*STD
            Vout(k) = MEAN;
        end
    end

    % 다시 phase correction
    Theta = fminsearch(@(X) phase_jump(Vout,Vin,time_exp,X),0);

    Vin_n  = Vin*cos(Theta) - Vout*sin(Theta);
    Vout_n = Vout*cos(Theta) + Vin*sin(Theta);

    Vin  = Vin_n;
    Vout = Vout_n;

  %% ===== TIME ZERO CORRECTION =====

    Vin_max = max(Vin);
    half_val = Vin_max/2;
    
    % peak 위치
    [~, idx_max] = max(Vin);
    
    % peak 이전 구간에서 half max 찾기 (rising edge)
    [~, idx_half] = min(abs(Vin(1:idx_max) - half_val));
    
    t_zero = time_exp(idx_half);
    
    % time shift
    time_exp = time_exp - t_zero;
    td_ps = time_exp*1e12;

    %% ===== ratio =====
    Vratio = -Vin ./ Vout;

    %% ===== legend 이름 =====
    [~, name, ~] = fileparts(filenames{i});
    token = regexp(name,'_(\d+)$','tokens');

    if ~isempty(token)
        legend_name = token{1}{1};
    else
        legend_name = name;
    end

    %% ===============================
    %% ===== plotting =====
    %% ===============================

    figure(1)
    plot(td_ps,Vin,...
        'LineWidth',lw,...
        'Color',colors(i,:),...
        'DisplayName',legend_name);

    figure(2)
    plot(td_ps,Vout,...
        'LineWidth',lw,...
        'Color',colors(i,:),...
        'DisplayName',legend_name);

    figure(3)

    idx_ratio = time_exp >= tdelay_min;

    plot(td_ps(idx_ratio),Vratio(idx_ratio),...
    'MarkerSize',5,...
    'LineWidth',2.5,...
    'Color',colors(i,:),...
    'DisplayName',legend_name);
end

%% ===== 마무리 =====
figure(1); legend show
figure(2); legend show
figure(3); legend show

xlim([100 max(td_ps)])