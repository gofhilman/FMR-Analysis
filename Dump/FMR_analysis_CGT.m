fileNames = ["CGT_CP2_allT.dat"];
data = [];
for it_var = 1:length(fileNames)
    data = [data; load(fileNames(it_var))];
end
% data = importdata("Ni_LiNbO_0deg.dat");
% data = table2array(TriNbRe20Py20NbRe20CP2allT2);
H_Col = 1;
freq_Col = 2;
S12_Col = 3;
T_Col = 4;

global next_freq;
global mode_manual;

next_freq = 0;
mode_manual = 1;

% Fem un vector amb el nombre de dades de cada experiment
it_var = 1;
clear("n_Mes");
n_Mes = [];
n_freq = [];
while sum(n_Mes) < size(data,1)
    data_left = data(sum(n_Mes)+1:end,:);

    [~, freq_exp] = max(data_left(2:end,freq_Col)<data_left(1:end-1,freq_Col));
    
    [~, index_var1] = max(data_left(2:end, H_Col) < data_left(1:end-1, H_Col));
    [~, index_var2] = max(data_left(2:end, H_Col) > data_left(1:end-1, H_Col));
    index_var = max(index_var1, index_var2);

    if index_var == freq_exp
        n_freq(it_var) = freq_exp;
        n_Mes(it_var) = size(data_left,1);
        break;
    end
    n_freq(it_var) = freq_exp;
    n_Mes(it_var) = index_var;
    it_var = it_var + 1;
end
% T_list(2:end) = [];
% n_Mes = size(data,1);
% n_freq(2:end) = [];
exp_start_index = [1];
for it_var = 1:length(n_Mes)
    exp_start_index(it_var + 1) = exp_start_index(it_var) + n_Mes(it_var);
end
freq_data = data(1:n_freq, freq_Col) * 1e-9; % en GHz

% Prenem un vector amb les diferents temperatures de mesura
T_list = [];
for it_var = 1:length(n_Mes)
    T_list(it_var) = mean(data(exp_start_index(it_var):exp_start_index(it_var + 1) - 1,T_Col));
end

fit_data = [];
%%
not_pH = 1;
% Iterem sobre les diferents temperatures
for exp_it_var = 4:length(n_Mes)
    % Prenem les dades per l'experiment número exp_it_var
    exp_data = data(exp_start_index(exp_it_var):exp_start_index(exp_it_var+1)-1,:);

    w = 0; c = 0;
    % ANÀLISI CAMPS NEGATIUS
    % Iterem sobre les freqüències de l'experiment
    for freq_it_var = 110:n_freq(exp_it_var)
%         try
%             if ~valid_freq(freq_it_var)
%                 continue;
%             end
%         catch
%             continue;
%         end
% 
%    -0.9288   24.2000
%    -0.9370   24.4000
%    -0.9434   24.6000
%    -0.9491   24.8000
%    -0.9555   25.0000
%    -0.9788   25.6000
%    -0.9856   25.8000
%    -0.9987   26.2000

        % Creem el vector de les freqüències
        freq_data = data(exp_start_index(exp_it_var):exp_start_index(exp_it_var)+n_freq(exp_it_var)-1, freq_Col) * 1e-9; % en GHz
        % Creem el vector de camps H i S12
        H_data = exp_data(freq_it_var:n_freq(exp_it_var):end, H_Col)./10;
        S12_data = exp_data(freq_it_var:n_freq(exp_it_var):end, S12_Col);

        % Divisió en metitat negativa de camp
        [~, zero_index] = min(abs(H_data));
        neg_H_data = H_data(1:zero_index);
        neg_S12_data = S12_data(1:zero_index);

%         % Intentem treure el background
        neg_S12_data_smooth = smoothdata(neg_S12_data,"gaussian");
        S12_data = neg_S12_data - neg_S12_data_smooth;

        % Anàlisi dades negatives
        [fit_Lorentz, fit_param] = fit_data_fcn(freq_data, neg_H_data, neg_S12_data, T_list(exp_it_var), freq_it_var, exp_it_var, c, w);

        if ~isempty(fit_Lorentz)
            c = fit_Lorentz.center; w = fit_Lorentz.w;
            intervals_inc= confint(fit_Lorentz);
            fit_data(freq_it_var,:,2*exp_it_var-1) = [freq_data(freq_it_var), c, w, abs(intervals_inc(1,3) - c), abs(intervals_inc(2,3) - c), abs(intervals_inc(1,6) - w), abs(intervals_inc(2,6) - w)];
        end
    end

    if not_pH
        continue;
    end
    w = 0; c = 0;
    % ANÀLISI CAMPS POSITIUS
    % Iterem sobre les freqüències de l'experiment
    for freq_it_var = 1:n_freq(exp_it_var)
        % Creem el vector de les freqüències
        freq_data = data(exp_start_index(exp_it_var):exp_start_index(exp_it_var)+n_freq(exp_it_var)-1, freq_Col) * 1e-9; % en GHz
        % Creem el vector de camps H i S12
        H_data = exp_data(freq_it_var:n_freq(exp_it_var):end, H_Col)./10;
        S12_data = exp_data(freq_it_var:n_freq(exp_it_var):end, S12_Col);

        % Divisió en metitat negativa de camp
        [~, zero_index] = min(abs(H_data));
        pos_H_data = H_data(zero_index:end);
        pos_S12_data = S12_data(zero_index:end);

%         % Intentem treure el background
        pos_S12_data_smooth = smoothdata(pos_S12_data,"gaussian");
        S12_data = pos_S12_data - pos_S12_data_smooth;

        % Anàlisi dades negatives
        [fit_Lorentz, fit_param] = fit_data_fcn(freq_data, pos_H_data, pos_S12_data, T_list(exp_it_var), freq_it_var, exp_it_var, c, w);

        if ~isempty(fit_Lorentz)
            c = fit_Lorentz.center; w = fit_Lorentz.w;
            intervals_inc= confint(fit_Lorentz);
            fit_data(freq_it_var,:,2*exp_it_var) = [freq_data(freq_it_var), c, w, abs(intervals_inc(1,3) - c), abs(intervals_inc(2,3) - c), abs(intervals_inc(1,6) - w), abs(intervals_inc(2,6) - w)];
        end
    end
end

%% KITTEL'S LAW PLOT (LEFT)

for data_set_it_var = 1:size(fit_data,3)

%     data_set_it_var = 11;
    center_data = nonzeros(abs(fit_data(:,2,data_set_it_var)));
    freq_data = nonzeros(abs(fit_data(:,1,data_set_it_var)));
    if length(freq_data) <= 3
        continue
    end
    ft = fittype( '13.99624354*g * sqrt(x * (x + 4*pi*Meff))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf 1.5]; % Meff g
    opts.Upper = [Inf 1e6];
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.StartPoint = [0 2];
    opts.TolFun = 1e-14;
    
    % Fit model to data.
    [kittel_fit, ~] = fit(center_data, freq_data, ft, opts );
    
    kittel_confint = confint(kittel_fit);
    index_var = ceil(data_set_it_var / 2);
    data_param_vals(mod(data_set_it_var+1,2) + 1, 1:7, index_var) = [T_list(index_var), kittel_fit.Meff, kittel_fit.g, abs(kittel_confint(1, 1) - kittel_fit.Meff), abs(kittel_confint(2, 1) - kittel_fit.Meff), abs(kittel_confint(1, 2) - kittel_fit.g), abs(kittel_confint(2, 2) - kittel_fit.g)]; 
    % Plot fit with data.
    figure( 'Name', ['Kittel plot. T = ', num2str(T_list(ceil(data_set_it_var/2))), ' K']);
    hold on
    errorbar(center_data, freq_data, nonzeros(fit_data(:, 4,data_set_it_var)), nonzeros(fit_data(:, 5,data_set_it_var)),'.', 'horizontal')
    h = plot( kittel_fit, center_data, freq_data );
    legend( h, 'Experimental data', "Kittel's Law fit", 'Location', 'NorthWest', 'Interpreter', 'none' );
    % Label axes
    title('Gyromagnetic fit');
    xlabel( '$H_r$ (T)', 'Interpreter', 'latex' );
    ylabel( '$f$ (GHz)', 'Interpreter', 'latex' );
    annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', [strcat("$M_{eff}$ = ", num2str(kittel_fit.Meff), " $\pm$ ", num2str(data_param_vals(mod(data_set_it_var+1,2) + 1, 4, index_var)), " T") strcat("$g$ = ", num2str(kittel_fit.g), " $\pm$ ", num2str(data_param_vals(mod(data_set_it_var+1,2) + 1, 5, index_var)))], 'Interpreter','latex')
    box on; grid on
    hold off
    pause(2)
    
    % close all;
end
%%
erase_freq = 7.8;
[~,freq_index] = min(abs(fit_data(:,1,data_set_it_var) - erase_freq));

if abs(fit_data(freq_index, 1, data_set_it_var) - erase_freq) > 0.1
    return
end

fit_data(freq_index, :, data_set_it_var) = [0 0 0 0 0 0 0];

%% KITTEL'S LAW PLOT (RIGHT)

for data_set_it_var = 1:size(fit_data,3)

    center_data = nonzeros(fit_data(:,2,data_set_it_var));
    freq_data = nonzeros(abs(fit_data(:,1,data_set_it_var)));
    if length(freq_data) <= 4
        continue
    end
    ft = fittype( 'gamma*((x + B0) - 4 * pi * Meff)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-1 0 27]; % B0 Meff gamma
    opts.Upper = [1 Inf 29];
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.StartPoint = [0 1.5 28];
    opts.TolFun = 1e-14;
    
    % Fit model to data.
    [kittel_fit, kittel_param] = fit(center_data, freq_data, ft, opts );
    
    kittel_confint = confint(kittel_fit);
    index_var = ceil(data_set_it_var / 2);
    data_param_vals(mod(data_set_it_var+1,2) + 1, 1:5, index_var) = [T_list(index_var), kittel_fit.B0, kittel_fit.Meff, abs(kittel_confint(1, 1) - kittel_fit.B0), abs(kittel_confint(1, 2) - kittel_fit.Meff)]; 

    % Plot fit with data.
    figure( 'Name', ['Kittel plot. T = ', num2str(T_list(mod(data_set_it_var+1,2) + 1)), ' K']);
    hold on
    % errorbar(center_data, freq_data, nonzeros(fit_data(:, 3,data_set_it_var)), nonzeros(fit_data(:, 3,data_set_it_var)), nonzeros(fit_data(:, 3,data_set_it_var)),'.')
    h = plot( kittel_fit, center_data, freq_data );
    legend( h, 'Experimental data', "Kittel's Law fit", 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    title('Kittel plot');
    xlabel( 'B (T)', 'Interpreter', 'none' );
    ylabel( 'f (GHz)', 'Interpreter', 'none' );
    annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', [strcat("$B_0$ = ", num2str(kittel_fit.B0), " $\pm$ ", num2str(data_param_vals(mod(data_set_it_var+1,2) + 1, 4, index_var)), " T") strcat("$M_{eff}$ = ", num2str(kittel_fit.Meff), " $\pm$ ", num2str(data_param_vals(mod(data_set_it_var+1,2) + 1, 5, index_var)), " T")], 'Interpreter','latex')
    box on; grid on
    hold off
    pause(2)
    
    % close all;
end
%%
erase_freq = 22.4;
[~,freq_index] = min(abs(fit_data(:,1,data_set_it_var) - erase_freq));

if abs(fit_data(freq_index, 1, data_set_it_var) - erase_freq) > 0.1
    return
end

fit_data(freq_index, :, data_set_it_var) = [0 0 0 0 0 0 0];
%% GILBERT'S DAMPING FACTOR
% clear("data_param_vals")
data_param_vals(2,1,1) = 0;
for data_set_it_var = 1:size(fit_data,3)

%     data_set_it_var = 11;
    index_var = ceil(data_set_it_var / 2);

    w_data = nonzeros(abs(fit_data(:,3,data_set_it_var))');
        freq_data = nonzeros(abs(fit_data(:,1,data_set_it_var))');
    if length(freq_data) <= 3
        continue
    end
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Lower = [0 0];
    
    % Fit model to data.
    [dampingfit, dampingparam] = fit( freq_data, w_data, ft, opts );
    
    alpha = 14 * dampingfit.p1;
    
    % Save results
    damping_confint = confint(dampingfit);

    data_param_vals(mod(data_set_it_var+1,2) + 1, 1, index_var) = [T_list(index_var)];
    data_param_vals(mod(data_set_it_var+1,2) + 1, 8:13, index_var) = [alpha, dampingfit.p2, 14 * damping_confint(1,1), 14 * damping_confint(2,1), damping_confint(1,2), damping_confint(2,2)]; 


    % Plot fit with data.
    figure( 'Name', ['Damping fit. T = ', num2str(T_list(index_var)), ' K'] );
    hold on;
    errorbar(freq_data, w_data, nonzeros(fit_data(:,6,data_set_it_var)), nonzeros(fit_data(:,7,data_set_it_var)), 'vertical' ,'.')
    h = plot( dampingfit, freq_data, w_data );
    legend( h, 'Experimental data', 'Linear fit', 'Location', 'NorthWest', 'Interpreter', 'none' );
    % Label axes
    if sign(sum(fit_data(:,2,data_set_it_var))) == 1
        title(['Damping factor. T = ', num2str(T_list(index_var)),' K. H > 0']);
    else
        title(['Damping factor. T = ', num2str(T_list(index_var)),' K. H < 0']);
    end
    
    xlabel( '$f$ (GHz)', 'Interpreter', 'latex' );
    ylabel( '$\mu_0\Delta H$ (T)', 'Interpreter', 'latex');
    annotation('textbox', [0.5, 0.15, 0.1, 0.1], 'String', [strcat("$\alpha$ = ", num2str(alpha), " $\pm$ ", num2str(abs(data_param_vals(mod(data_set_it_var+1,2) + 1, 10, index_var) - alpha))) strcat("$\alpha_{inh}$ = (", num2str(dampingfit.p2), " $\pm$ ", num2str(abs(data_param_vals(mod(data_set_it_var+1,2) + 1, 12, index_var) - dampingfit.p2)), ") T")], 'Interpreter','latex')
    box on; grid on
    hold off;

    pause(2)
    
%     close all;
end

%%
function [fit_Lorentz, fit_param] = fit_data_fcn(freq_data, H_data, S12_data, T_val, freq_it_var, exp_it_var, c, w)
    global mode_manual;
    global next_freq;

    if mode_manual == 1 || next_freq == 1 || freq_it_var == 1 || (c == 0 && w == 0)
        h = figure();
        plot(H_data, S12_data);
        title(['Freq = ', num2str(freq_data(freq_it_var)), ' GHz. T = ', num2str(T_val), ' K']);
        xlabel( 'H (T)', 'Interpreter', 'none' );
        ylabel( 'S12 (dB)', 'Interpreter', 'none' );
        try
            next_freq = 0;
            [x, ~] = ginput(2);
            close(h);
        catch
            next_freq = 1;
            fit_Lorentz = [];
            fit_param = [];
            return
        end
    else
        x = [c - 3 * w, c + 3 * w];
    end

    [~, h1_index] = min(abs(H_data - x(1)));
    [~, h2_index] = min(abs(H_data - x(2)));
    min_index = min(h1_index, h2_index); max_index = max(h1_index, h2_index);

    H_data_cropped = H_data(min_index:max_index);
    S12_data_cropped = S12_data(min_index:max_index);

    % Si no hi ha suficients punts, introduïm l'interval manualment
    if length(H_data_cropped) < 10
        [fit_Lorentz, fit_param] = fit_data_fcn(freq_data, H_data, S12_data, T_val, freq_it_var, exp_it_var, 0, 0);
        return
    end
    
    [S12_min, S12_min_index] = min(S12_data_cropped(ceil(0.15*length(S12_data_cropped)):floor(0.85*length(S12_data_cropped))));
    if c ~= 0 && (c < min(H_data) || c > max(H_data))
        c = H_data_cropped(S12_min_index);
    end
    if w ~= 0 && w < 1e-4
        w = 0.002;
    end
    if w == 0 && c == 0
        w = 0.002;
        c = H_data_cropped(S12_min_index);
    end

    a = 1/2 * w * abs(max(S12_data_cropped(end), S12_data_cropped(1)) - S12_min);
    b = 2 * a;
    k = (S12_data_cropped(end) - S12_data_cropped(1)) / (H_data_cropped(end) - H_data_cropped(1)); 
    off = (S12_data_cropped(1) + S12_data_cropped(end)) / 2;

    ft = fittype( 'A/(4*(x-center)^2+w^2)+offset-(w*B*(x-center))/(4*(x-center)^2+w^2)+k*x', 'independent', 'x', 'dependent', 'y' );
        
    startpoints = [a b c k off w];
    opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
                  'Algorithm','Trust-Region',...                
                  'Robust','off',...
                  'Startpoint', startpoints,...
                  'MaxFunEvals', 5000 ,...
                  'MaxIter', 5000 ,...
                  'TolFun', 10^-14);
    opts.Display = 'Off';
    
    % Fit model to data.
    [fit_Lorentz, fit_param] = fit( H_data_cropped, S12_data_cropped, ft, opts );

    figure();
    hold on;
        scatter( H_data_cropped, S12_data_cropped, 'Marker', '.' );
        plot(linspace(H_data_cropped(1),H_data_cropped(end),500), fit_Lorentz(linspace(H_data_cropped(1),H_data_cropped(end),500)));
        title(['Freq = ', num2str(freq_data(freq_it_var)), ' GHz. T = ', num2str(T_val), ' K']);
        xlabel( 'H (T)', 'Interpreter', 'none' );
        ylabel( 'S12 (dB)', 'Interpreter', 'none' );
    hold off;
    if fit_param.rsquare < 0.85
        [fit_Lorentz, fit_param] = fit_data_fcn(freq_data, H_data, S12_data, T_val, freq_it_var, exp_it_var, 0, 0);
    end
    pause(0.5);
    close(gcf);
end