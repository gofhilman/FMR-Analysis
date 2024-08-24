% Set the column number: Frequency, Field, S12, Index
fp = FileParameters(3, 2, 4, 0);

% Read experiment from file
FGT_exp = ExperimentCollection(["FGT_perp_RT01.txt", ...
    "FGTperp03.txt", "FGTperp04.txt", "FGTinplane05.txt", ...
    "FGT_25grade_06.txt", "FGT_25grade_07.txt", "FGT_45grade_08.txt", ...
    "FGT_60grade_09.txt", "FGT_80grade_10.txt"], fp);

% Set units of the data in file
FGT_exp.setfieldunits(MagneticUnit.oersted);
FGT_exp.setfrequencyunits(FrequencyUnit.hertz);

% Convert units
FGT_exp.convertfieldunits(MagneticUnit.tesla);
FGT_exp.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);

% Filter by positive or negative field
%   If you want the positive field, comment the 
%   first line and uncomment the second line
% FGT_exp.filternegativefield();
% FGT_exp.filterpositivefield();

%% Remove experiment
removeexperiment(FGT_exp, 12);

%% Add another one
new_experiment = Experiment("Fe2O3Pt_RT_01.txt", fp);
new_experiment.setfieldunits(MagneticUnit.oersted);
new_experiment.setfrequencyunits(FrequencyUnit.hertz);
new_experiment.convertfieldunits(MagneticUnit.tesla);
new_experiment.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);
addexperiment(FGT_exp, new_experiment);
FGT_exp.filternegativefield();

%% Make image plot
% When the plot pops up, press ENTER to be able to
% zoom in, and modify the plot options
for i = 1:1
    FGT_exp.experimentArray(i).getcrossingpoints(Filter.normalization);
    xlabel("$\mu_0 H$ (T)", 'Interpreter','latex');
    ylabel("Frequency (GHz)");
    exportgraphics(gca, "imageplot\figure_"+i+".png");
end

%% Test
x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298,298];
save = zeros(1,9);
for i = 1:12
    xData = FGT_exp.experimentArray(1, i).crossingPoints(:,1);
    yData = FGT_exp.experimentArray(1, i).crossingPoints(:,2);
    ft = fittype( 'gamma*sqrt(x*(x+dmF)+2*eaF)', 'independent', 'x', ...
                  'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = [-5 0 28];
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.TolFun = 1e-14;
    opts.Lower = [-Inf -Inf 28.024];
    opts.Upper = [Inf Inf 28.025];
    [fitresult, gof] = fit(xData, yData, ft, opts)
    confInt = confint(fitresult);
    dmConfInt = 0.5 * abs(confInt(1,1) - confInt(2,1));
    eaConfInt = 0.5 * abs(confInt(1,2) - confInt(2,2));
    gyroConfInt = 0.5 * abs(confInt(1,3) - confInt(2,3));
    figure( 'Name', 'Kittel fit' );
    plot(fitresult, xData, yData);
    xlabel(['$\mu_0 H$ (T)'], 'Interpreter', 'latex');
    ylabel(['Frequency (GHz)'], 'Interpreter', 'latex');
    title(['r^2 = ', num2str(gof.rsquare)]);

    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$H_{DM}$ = ", num2str(fitresult.dmF), ...
        " $\pm$ ", num2str(dmConfInt), " T"), ...
        strcat("$H_{E} H_{A}$ = ", ...
        num2str(fitresult.eaF), " $\pm$ ", ... 
        num2str(eaConfInt), " T^2")], 'Interpreter','latex');
    
    %{
annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        {sprintf("$H_{DM}$ = %.4f $/pm$ %.5f T", fitresult.dmF, dmConfInt), ...
        sprintf("$H_{E} H_{A}$ = %.5f $/pm$ %.6f T^2", fitresult.eaF, eaConfInt)}, 'Interpreter','latex');
    exportgraphics(gca, "kittelplot\kittel_alt_"+i+".png");
    %}
    save(i) = fitresult.dmF;
end

%% ABc
figure( 'Name', 'DM field vs temperature' );
plot(x,save,'.','MarkerSize',25);
title('H_{DM} vs temperature');
xlabel('Temperature (K)');
ylabel('H_{DM} (T)');
% xticks(0:20:max(extracted_data(:,1)));
% xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
% exportgraphics(gca, "DMFieldVsTemperature.png");

%% Plot S12 vs Field for each frequency and select peaks
% To go to the next frequency select a range for the peak,
% press ENTER or close the window
for i = 12:12
    FGT_exp.resonanceanalysis(i);
end

%% Make kittel plot
% Magnetization can be set to inPlane, outOfPlane, paramagnetic, aFM
for i = 1:12
    % if (i~=2) && (i~=7)
    FGT_exp.experimentArray(i).setmagnetization(Magnetization.aFM);
    FGT_exp.experimentArray(i).makekittel();
    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$H_{DM}$ = ", num2str(FGT_exp.experimentArray(1, i).kittelParameters.dmField.value), ...
        " $\pm$ ", num2str(FGT_exp.experimentArray(1, i).kittelParameters.dmField.error), " ", ...
        FGT_exp.experimentArray(1, i).kittelParameters.dmField.unit.tag) strcat("$H_{E} H_{A}$ = ", ...
        num2str(FGT_exp.experimentArray(1, i).kittelParameters.eaField.value), " $\pm$ ", ... 
        num2str(FGT_exp.experimentArray(1, i).kittelParameters.eaField.error), " ", ...
        FGT_exp.experimentArray(1, i).kittelParameters.eaField.unit.tag)], 'Interpreter','latex');
    exportgraphics(gca, "kittelplot\kittel_"+i+".png");
    % end
end

%% Make damping plot
for i = 1:12
    FGT_exp.experimentArray(i).makedampingfit();
    annotation('textbox', [0.5, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$\alpha$ = ", num2str(FGT_exp.experimentArray(1, i).dampingParameters.damping.value), ...
        " $\pm$ ", num2str(FGT_exp.experimentArray(1, i).dampingParameters.damping.error)) strcat("$\alpha_{inhomogenous}$ = ", ...
        num2str(FGT_exp.experimentArray(1, i).dampingParameters.inhomogeneousDamping.value), " $\pm$ ", ... 
        num2str(FGT_exp.experimentArray(1, i).dampingParameters.inhomogeneousDamping.error))], 'Interpreter','latex');
    exportgraphics(gca, "dampingplot\damping_"+i+".png");
end

%% vs temperature
x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298,298];
extracted_data = zeros(9,9);
j=1;
for i = 1:12
    if (i~=3) && (i~=7) && (i~=8)
        extracted_data(j,:) = [x(i), FGT_exp.experimentArray(1,i).kittelParameters.dmField.value, ...
            FGT_exp.experimentArray(1,i).kittelParameters.dmField.error, ...
            FGT_exp.experimentArray(1,i).dampingParameters.damping.value, ...
            FGT_exp.experimentArray(1,i).dampingParameters.damping.error, ...
            FGT_exp.experimentArray(1,i).kittelParameters.eaField.value, ...
            FGT_exp.experimentArray(1,i).kittelParameters.eaField.error, ...
            FGT_exp.experimentArray(1,i).dampingParameters.inhomogeneousDamping.value, ...
            FGT_exp.experimentArray(1,i).dampingParameters.inhomogeneousDamping.error];
        j = j+1;
    end
end

%% DM field vs temperature
figure( 'Name', 'DM field vs temperature' );
errorbar(extracted_data(:,1),extracted_data(:,2),extracted_data(:,3),'.','MarkerSize',25);
title('H_{DM} vs temperature');
xlabel('Temperature (K)');
ylabel('H_{DM} (T)');
xticks(0:20:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "DMFieldVsTemperature.png");

%% Damping vs temperature
figure( 'Name', 'Damping vs temperature' );
errorbar(extracted_data(:,1),extracted_data(:,4),extracted_data(:,5),'.','MarkerSize',15);
title('\alpha vs temperature');
xlabel('Temperature (K)');
ylabel('\alpha');
xticks(0:20:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "DampingVsTemperature.png");

%% EA field vs temperature
figure( 'Name', 'EA field vs temperature' );
errorbar(extracted_data(:,1),extracted_data(:,6),extracted_data(:,7),'.','MarkerSize',10);
title('H_{E} H_{A} vs temperature');
xlabel('Temperature (K)');
ylabel('H_{E} H_{A} (T^2)');
xticks(0:20:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "EAFieldVsTemperature.png");

%% Inhomogeneous Damping vs temperature
figure( 'Name', 'Inhomogeneous damping vs temperature' );
errorbar(extracted_data(:,1),extracted_data(:,8),extracted_data(:,9),'.','MarkerSize',15);
title('Inhomogeneous damping vs temperature');
xlabel('Temperature (K)');
ylabel('Inhomogeneous damping');
xticks(0:20:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "InhonogeneousDampingVsTemperature.png");

%% Make kittel plot
% Magnetization can be set to inPlane, outOfPlane, paramagnetic, aFM
for i = 1:12
    % if (i~=2) && (i~=7)
    FGT_exp.experimentArray(i).setmagnetization(Magnetization.inPlane);
    FGT_exp.experimentArray(i).makekittel();
    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$M_{eff}$ = ", num2str(FGT_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.value), ...
        " $\pm$ ", num2str(FGT_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.error), " ", ...
        FGT_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.unit.tag) strcat("$H_{A}$ = ", ...
        num2str(FGT_exp.experimentArray(1, i).kittelParameters.anisotropyField.value), " $\pm$ ", ... 
        num2str(FGT_exp.experimentArray(1, i).kittelParameters.anisotropyField.error), " ", ...
        FGT_exp.experimentArray(1, i).kittelParameters.anisotropyField.unit.tag)], 'Interpreter','latex');
    exportgraphics(gca, "kittelplot_alt\kittel_"+i+".png");
    % end
end

%% vs temperature (alt)
x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298,298];
extracted_data_alt = zeros(9,5);
j=1;
for i = 1:12
    if (i~=3) && (i~=7) && (i~=8)
        extracted_data_alt(j,:) = [x(i), FGT_exp.experimentArray(1,i).kittelParameters.effectiveMagnetization.value, ...
            FGT_exp.experimentArray(1,i).kittelParameters.effectiveMagnetization.error, ...
            FGT_exp.experimentArray(1,i).kittelParameters.anisotropyField.value, ...
            FGT_exp.experimentArray(1,i).kittelParameters.anisotropyField.error];
        j = j+1;
    end
end

%% Meff vs temperature
figure( 'Name', 'Meff vs temperature' );
errorbar(extracted_data_alt(:,1),extracted_data_alt(:,2),extracted_data_alt(:,3),'.','MarkerSize',25);
title('M_{eff} vs temperature');
xlabel('Temperature (K)');
ylabel('M_{eff} (T)');
xticks(0:20:max(extracted_data_alt(:,1)));
xlim([min(extracted_data_alt(:,1))-10 max(extracted_data_alt(:,1))+10]);
exportgraphics(gca, "MeffVsTemperature.png");

%% H ani vs temperature
figure( 'Name', 'Hani vs temperature' );
errorbar(extracted_data_alt(:,1),extracted_data_alt(:,4),extracted_data_alt(:,5),'.','MarkerSize',10);
title('H_{ani} vs temperature');
xlabel('Temperature (K)');
ylabel('H_{ani} (T)');
xticks(0:20:max(extracted_data_alt(:,1)));
xlim([min(extracted_data_alt(:,1))-10 max(extracted_data_alt(:,1))+10]);
exportgraphics(gca, "HaniVsTemperature.png");