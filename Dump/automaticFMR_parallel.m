%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Automatic FMR loop                         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for FMR fits
filen = 'devices_D.txt';
data = load(filen);
% Find the number of frequencies
[v, indf] = max(data(:,2), [],'omitnan');
nFreq = indf;
deltaFreq = 1; 
% 'y' = yes, 'n' = no
fullSpectrum = 'y'; % Do you have full spectrum (negative to positive fields)?
positive = 'y'; % Do you want to analyze positive or negative fields?
rangeX1 = 1; % Multiplier on Lorentizan fits
rangeX2 = 4; % Multiplier on Lorentizan fits
sampleText = filen(1:end-4); % Text in plot
degrees = 0;
aind = 1;
% Makes folder to save figures
mkdir(['Figures_' filen(1:end-4)])
fileDir = strcat(pwd,['\Figures_' filen(1:end-4)]);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

dataFit = [];
iniFreq = 22; 
finFreq = nFreq;
nFields = length(data)/nFreq;
ind = 1;
fieldOffset = 0;
% Find the 0 field index
for i = iniFreq:deltaFreq:finFreq
    freq = i  % index value of frequency
    val = smooth(data(freq:nFreq:end,4),1);
    h   = data(freq:nFreq:end,1);
    [d1, loc0] = min(abs(0-h));
    [~, iminh] = min(h);
    if fullSpectrum == 'y'
        if iminh < length(h)/2
            if positive == 'y'
                ini = loc0; % Index of first positive field
                rest = 0;
            else
                ini = 1;
                rest = loc0;
            end
        else
            if positive == 'y'
                 ini = 1;
                rest = loc0;
            else
                ini = loc0; % Index of first positive field
                rest = 0;
            end
        end
    end
%     rest = 0;
%     ini = 0;
    frequency = data(freq,2)*10^-9;
    if ini > 1 || rest > 1
        val = val(ini:end-rest);
        h   = h(ini:end-rest);
    end
    x2max = round(length(h)-rest);
    % Choose initial points
    if ind==1 
        plot(detrend(val),'.-')
        [x,y] = ginput(2);
        close;
        % Range of points in x
        x1 = round(x(1));
        x2 = round(x(2));
        dis = round((x2-x1)/2); % Approx. distance between center and each extreme
    end
    s21 = detrend(val(x1(1):x2(1),1));
    field = h(x1(1):x2(1),1);
    % Using kittels formula find b center:
    % Find minimum (detrend data, find minimum, find the index of it)
    s21d = detrend(s21);
    [v, indt] = min(s21d, [],'omitnan'); % find minimum in s21
    minField = field(indt);
    [d, ix] = min(abs(minField-h));
    [d1, ix1] = min(abs(minField-field));
    % Index of minimum in whole data
    if ind > 1
        [d, ix] = min(abs(fitLorent.center-h))
        fcenter = fitLorent.center
        w = fitLorent.w;
    else
        fcenter = field(ix1);
        w=0.005;
    end
    %---------------------------- Lorentizian Fitting ---------------------
    ft = fittype( ['A*w/(4*(x - center)^2+w^2) + offset - (w*B*(x-center)) ...' ...
        '(/4*(x - center)^2+w^2) + k*x'], 'independent', 'x', 'dependent', 'y' );
    c = fcenter;                                
    a = ((max(s21)-min(s21))*w)/2;            
    off = max(s21); b = a*2;k = 1; 
    startpoints = [a b c k off w];
    lower = [-2 -2 0   -5 -50 0.0001];
    upper = [ 2  2 0.5  5 10 2];                
    opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
                      'Algorithm','Trust-Region',...                
                      'Robust','off',...
                      'Startpoint', startpoints,...
                      'Lower', lower,...
                      'Upper',upper,...
                      'MaxFunEvals', 10000 ,...
                      'MaxIter', 10000 ,...
                      'TolFun', 10^-14);
    [fitLorent, gof] = fit(field, s21, ft, opts );
    fits{i} = fitLorent;
    plot( fitLorent, field, s21 );
%     plot( fitLorent)
%         hold on
    legend('S21 vs H', 'Lorentzian fit');
    xlabel( 'Field (T)'); ylabel( 'S21 (dB)');
    title(['FMR Frequency: ' num2str(frequency) 'GHz']);
    box on; grid on;
    pause
    %---------------------------- Lorentizian Fitting ---------------------
   % Save data only if its good
    [intervals] = confint(fitLorent); % Errors
    errorInt = abs(abs(intervals(1,6))-abs(intervals(2,6)));
    if gof.rsquare > 0.9 && errorInt < 0.25 && fitLorent.w > 1e-3
        dataFit(ind,1) = frequency*10^9;
        dataFit(ind,2) = abs(fitLorent.w);
        dataFit(ind,3) = fitLorent.center;
        dataFit(ind,4) = gof.rsquare;
        dataFit(ind,5) = errorInt;
    else
        ind = ind-1;
    end
    % Readjust parameters for next fit...
    if length(ix)>1
        ix = max(ix);
    end
    x1 = ix - dis*(rangeX1); x1=round(x1);
    x2 = ix + dis*(rangeX2); x2=round(x2)
    ind = ind + 1;
    % Safeguards!
    if x2 > x2max-1
        x2 = x2max-1;
    end
    if x1 < 1
        x1 = 1;
    end
end

%% Extract damping (DeltaH vs freq)
x = dataFit(:,1)*1e-9; % Frequency
% ffield = transpose(ff(round(x)));
y = dataFit(:,2); % Delta_H
err = dataFit(:,5);

des = 3;
men = 0; 


ft = fittype( {'x', '1'}, 'independent', 'x', 'dependent', 'y', 'coefficients', {'a', 'b'} );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
[fitL, gof] = fit( x(des:end-men), y(des:end-men), ft, opts );
figure1 = figure;
axes1 = axes('Parent',figure1);
hold on
errorbar(x,y,ones(size(y)).*err./2,'. blue','MarkerSize',18)
% xlim([0 max(x)])
% ylim([0 max(y)])
plot(fitL)   
hold off
grid on
xlabel( 'Frequency (GHz)', 'Interpreter', 'none','FontSize',20 );
ylabel( '\Delta B(T)', 'Interpreter', 'tex','FontSize',20 );

set(gca,'FontSize',20);box on; grid on
set(axes1,'FontSize',20,'LineWidth',1);
slope = fitL.a; %Slope... = 2*alhpa/(gyro (Hz/T))
gyro = 28;
alpha = (fitL.a*gyro)/2
alpha_0 = fitL.b
Legend{1} = [sampleText ', \alpha =' num2str(alpha)];
legend(Legend ,'Location','northwest','FontSize',15)
format long
[intervals] = confint(fitL);
error = abs(abs(intervals(1,1))-abs(intervals(2,1)));
savefig(figure1, fullfile(fileDir, ['alphaPlot_' num2str(degrees) '.fig']))
saveas(figure1, fullfile(fileDir, ['alphaPlot_' num2str(degrees) '.png']))
filename = ['DampingFit' '.mat'];
save(filename,'fitL','gof','intervals','dataFit')

%% Gyromagnetic fit (Kittel)
start = 1;
final = 0;
mu = 4*pi*1*10^-7; % [mu] = H/m = T*m/A
x = dataFit(start:end-final,3);  % [B] = T 
if positive == 'n'
    x = dataFit(start:end-final,3)*-1;
end
y = dataFit(start:end-final,1);  % frequency
y2 = y; y2 = y2(1:end,1)*10^-9;
figure2 = figure;
axes1 = axes('Parent',figure2);
% Kittel fit
ft = fittype( 'a*sqrt((1*x+b)*(1*x+b+c))', 'independent', 'x', 'dependent', 'y' );
lower = [-10 -1 -Inf];
startpoints = [0 0 0];
upper = [100 1 Inf];
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
                      'Algorithm','Trust-Region',...                
                      'Robust','off',...
                      'Startpoint', startpoints,...
                      'MaxFunEvals', 5000 ,...
                      'MaxIter', 5000 ,...
                      'TolFun', 10^-14,...
                      'Lower',lower,...
                      'Upper',upper);

% Fit model to data.
[fit1, gof] = fit( x, y2, ft, opts );
hold on
% xlim([0 0.20]); ylim([0 16])
plot(fit1); 
% errorbar(x,y2,ones(size(y2(:,1))).*(1-dataFit(start:end-final,4).^2),'^ blue','MarkerSize',5,'MarkerFaceColor','blue')
plot(x,y2,'. blue','MarkerSize',20,'MarkerFaceColor','blue')
hold off
box on; set(axes1,'FontSize',20,'LineWidth',1);
xlabel('$\mu_0 H$ (T)','Interpreter','latex'); ylabel('$f$ (GHz)','Interpreter','latex'); set(gca,'FontSize',20)
legend('Fitted curve','Data Points','Location','southeast','FontSize',20)
gyro = fit1.a; B0 = fit1.b; muMs = fit1.c; Gyro = fit1.a;
% Errors
[intervals] = confint(fit1);
deltab = abs(abs(intervals(1,2))-abs(intervals(2,2)));
deltac = abs(abs(intervals(1,3))-abs(intervals(2,3)));
filename = strcat(fileDir, ['\kittelResults' num2str(degrees) '.txt']);
fileID = fopen(filename,'w');
fprintf(fileID,'Gyromagnetic Factor = %f GHz/T\n', gyro);
fprintf(fileID,'B_0 = %f +- %f T \n', B0, deltab);
fprintf(fileID,'muM_s = %f +- %f T \n', muMs, deltac);
fclose(fileID);
savefig(figure2, fullfile(fileDir, ['gyroMag_' num2str(degrees) '.fig']))
saveas(figure2, fullfile(fileDir, ['gyroMag_' num2str(degrees) '.png']))


%% FMR imageSC
fs = data(:,2);
flds = data(:,1);
T = data(:,4);
nfields = round(length(T)/nFreq);
% x and y axis values
fields = flds(1:nFreq:end,1);
freq = fs(1:nFreq);
% Reshape to nfreq x nfields
matn = reshape(T,[],nfields);
[d, m] = size(matn);

% Normalize the data for each frequency
for i = 1:d
   ma = max(matn(i,:));
   matn(i,:) = matn(i,:)-ma;
end
figure3 = figure;
imagesc(fields,freq,matn)
xlabel('Magnetic field (T)')
ylabel('Frequency (Hz)')
%title('Hard Axis')
axis xy
savefig(figure3, fullfile(fileDir, ['imageSC_' num2str(degrees) '.fig']))
saveas(figure3, fullfile(fileDir, ['imageSC_' num2str(degrees) '.png']))
% 
% alphaList(aind,1) = degrees;
% alphaList(aind,2) = alpha;
% alphaList(aind,3) = error;

%%
p = 15;
j =1;
hold on
for i = 2:4:34
    plot(flds(i:34:end), detrend(T(i:34:end)),'LineWidth',1.3)
    Legend{j} = strcat(num2str(dataWG(i,2)*10^-9),'GHz');
    j = j+1;
end
legend(Legend)
xlabel('Magnetic field (mT)')
ylabel('S21 (arb. u.)')
box on; grid on;
hold off


%%
dataFitP211117 = dataFit;

%%
%% Alpha plot
alphaList = sortrows(alphaList,1);
x = alphaList(:,1);
y = alphaList(:,2);
errorbar(x,y,ones(size(y(:,1))).*alphaList(:,3)./2,'.-','MarkerSize',20)

% plot(alphaList(:,1),alphaList(:,2),'.-','MarkerSize',20)
xlabel('Degrees (\circ)','FontSize',20);
ylabel('\alpha','FontSize',20 ); 