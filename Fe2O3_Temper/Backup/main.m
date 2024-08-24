% Set the column number: Frequency, Field, S12, Index
fp = FileParameters(3, 2, 4, 0);

% Read experiment from file
Fe2O3_exp = ExperimentCollection(["Fe2O3Pt_30K_09.txt", "Fe2O3Pt_04_40K.txt", ...
    "Fe2O3Pt_04_57K.txt", "Fe2O3Pt_04_75K.txt", "Fe2O3Pt_04_92K.txt", ...
    "Fe2O3Pt_04_109K.txt", "Fe2O3Pt_04_126K.txt", "Fe2O3Pt_04_143K.txt", ...
    "Fe2O3Pt_185K_03.txt", "Fe2O3Pt_220K_06.txt", "FMR_Fe2O3Pt_RT_edited.txt"], fp);

% Set units of the data in file
Fe2O3_exp.setfieldunits(MagneticUnit.oersted);
Fe2O3_exp.setfrequencyunits(FrequencyUnit.hertz);

% Convert units
Fe2O3_exp.convertfieldunits(MagneticUnit.tesla);
Fe2O3_exp.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);

% Filter by positive or negative field
%   If you want the positive field, comment the 
%   first line and uncomment the second line
% Fe2O3_exp.filternegativefield();
% Fe2O3_exp.filterpositivefield();
% Will be done later

%% Remove experiment
removeexperiment(Fe2O3_exp, 11);

%% Add another one
new_experiment = Experiment("Fe2O3Pt_RT_01.txt", fp);
new_experiment.setfieldunits(MagneticUnit.oersted);
new_experiment.setfrequencyunits(FrequencyUnit.hertz);
new_experiment.convertfieldunits(MagneticUnit.tesla);
new_experiment.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);
addexperiment(Fe2O3_exp, new_experiment);
Fe2O3_exp.filternegativefield();

%% Add one more
new_experiment = Experiment("Fe2O3Pt_RT_detailed.txt", fp);
new_experiment.setfieldunits(MagneticUnit.oersted);
new_experiment.setfrequencyunits(FrequencyUnit.hertz);
new_experiment.convertfieldunits(MagneticUnit.tesla);
new_experiment.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);
addexperiment(Fe2O3_exp, new_experiment);
% Fe2O3_exp.filternegativefield();

%% Make image plot
% When the plot pops up, press ENTER to be able to
% zoom in, and modify the plot options
for i = 1:11
    Fe2O3_exp.experimentArray(i).getcrossingpoints(Filter.noiseCancelling);
    xlabel("$\mu_0 H$ (T)", 'Interpreter','latex');
    ylabel("Frequency (GHz)"); 
    exportgraphics(gca, "imageplot\figure_nc_"+i+".png");
end

%% Plot S12 vs Field for each frequency and select peaks
% To go to the next frequency select a range for the peak,
% press ENTER or close the window
for i = 11:11
    Fe2O3_exp.resonanceanalysis(i);
end

%% Make kittel plot
% Magnetization can be set to inPlane, outOfPlane, paramagnetic, aFM
for i = 1:11
    % if (i~=2) && (i~=7)
    Fe2O3_exp.experimentArray(i).setmagnetization(Magnetization.aFM);
    Fe2O3_exp.experimentArray(i).makekittel();
    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$H_{DM}$ = ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.dmField.value), ...
        " $\pm$ ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.dmField.error), " ", ...
        Fe2O3_exp.experimentArray(1, i).kittelParameters.dmField.unit.tag) strcat("$H_{E} H_{A}$ = ", ...
        num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.eaField.value), " $\pm$ ", ... 
        num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.eaField.error), " ", ...
        Fe2O3_exp.experimentArray(1, i).kittelParameters.eaField.unit.tag)], 'Interpreter','latex');
    exportgraphics(gca, "kittelplot\kittel_"+i+".png");
    % end
end

%% Make damping plot
for i = 1:11
    Fe2O3_exp.experimentArray(i).makedampingfit();
    annotation('textbox', [0.5, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$\alpha$ = ", num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.damping.value), ...
        " $\pm$ ", num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.damping.error)) strcat("$\alpha_{inhomogenous}$ = ", ...
        num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.inhomogeneousDamping.value), " $\pm$ ", ... 
        num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.inhomogeneousDamping.error))], 'Interpreter','latex');
    exportgraphics(gca, "dampingplot\damping_"+i+".png");
end

%% vs temperature
x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298];
extracted_data = zeros(8,9);
j=1;
for i = 1:11
    if (i~=3) && (i~=7) && (i~=8)
        extracted_data(j,:) = [x(i), Fe2O3_exp.experimentArray(1,i).kittelParameters.dmField.value, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.dmField.error, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.damping.value, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.damping.error, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.eaField.value, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.eaField.error, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.inhomogeneousDamping.value, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.inhomogeneousDamping.error];
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
for i = 1:11
    % if (i~=2) && (i~=7)
    Fe2O3_exp.experimentArray(i).setmagnetization(Magnetization.inPlane);
    Fe2O3_exp.experimentArray(i).makekittel();
    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$M_{eff}$ = ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.value), ...
        " $\pm$ ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.error), " ", ...
        Fe2O3_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.unit.tag) strcat("$H_{A}$ = ", ...
        num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.anisotropyField.value), " $\pm$ ", ... 
        num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.anisotropyField.error), " ", ...
        Fe2O3_exp.experimentArray(1, i).kittelParameters.anisotropyField.unit.tag)], 'Interpreter','latex');
    exportgraphics(gca, "kittelplot_alt\kittel_"+i+".png");
    % end
end

%% vs temperature (alt)
x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298];
extracted_data_alt = zeros(8,5);
j=1;
for i = 1:11
    if (i~=3) && (i~=7) && (i~=8)
        extracted_data_alt(j,:) = [x(i), Fe2O3_exp.experimentArray(1,i).kittelParameters.effectiveMagnetization.value, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.effectiveMagnetization.error, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.anisotropyField.value, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.anisotropyField.error];
        j = j+1;
    end
end

%% Meff vs temperature
figure( 'Name', 'Meff vs temperature' );
errorbar(extracted_data_alt(:,1),extracted_data_alt(:,2),extracted_data_alt(:,3),'.','MarkerSize',10);
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

%% Alternative method
x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298];
save = zeros(11,4);
for i = 1:11
    xData = Fe2O3_exp.experimentArray(1, i).crossingPoints(:,1);
    yData = Fe2O3_exp.experimentArray(1, i).crossingPoints(:,2);
    ft = fittype( 'gamma*sqrt(x*(x+dmF)+2*eaF)', 'independent', 'x', ...
                  'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = [0 0 28];
    opts.Display = 'Off';
    opts.MaxFunEvals = 5000;
    opts.MaxIter = 5000;
    opts.TolFun = 1e-14;
    opts.Lower = [-Inf -Inf 28];
    opts.Upper = [Inf Inf 28.0001];
    [fitresult, gof] = fit(xData, yData, ft, opts)
    confInt = confint(fitresult);
    dmConfInt = 0.5 * abs(confInt(1,1) - confInt(2,1));
    eaConfInt = 0.5 * abs(confInt(1,2) - confInt(2,2));
    gyroConfInt = 0.5 * abs(confInt(1,3) - confInt(2,3));
    figure( 'Name', 'Kittel fit' );
    plot(fitresult, xData, yData);
    xlabel('$\mu_0 H$ (T)', 'Interpreter', 'latex');
    ylabel('Frequency (GHz)', 'Interpreter', 'latex');
    title(['r^2 = ', num2str(gof.rsquare)]);
    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("H_{DM} = ", num2str(fitresult.dmF), ...
        " \pm ", num2str(dmConfInt), " T"), ...
        strcat("H_{E} H_{A} = ", ...
        num2str(fitresult.eaF), " \pm ", ... 
        num2str(eaConfInt), " T^2")], 'Interpreter','tex');
    exportgraphics(gca, "kittelplot\kittel_m2_"+i+".png");

    Fe2O3_exp.experimentArray(i).trialanderror(Filter.noiseCancelling, ...
        fitresult.dmF, fitresult.eaF, fitresult.gamma);
    clim([-0.1 0]); % If using noise cancelling
    exportgraphics(gca, "kittelplot\kittel_m2_2dimg_"+i+".png");
    save(i,1) = fitresult.dmF;
    save(i,2) = dmConfInt;
    save(i,3) = fitresult.eaF;
    save(i,4) = eaConfInt;
end

%% DM field vs temperature (alt)
figure( 'Name', 'DM field vs temperature' );
errorbar(x,save(:,1),save(:,2),'.','MarkerSize',10);
title('H_{DM} vs Temperature');
xlabel('Temperature (K)');
ylabel('H_{DM} (T)');
xticks(0:15:max(x));
xlim([min(x)-10 max(x)+10]);
exportgraphics(gca, "DMFieldVsTemperature_alt.png");

%% EA field vs temperature (alt)
figure( 'Name', 'EA field vs temperature' );
errorbar(x,save(:,3),save(:,4),'.','MarkerSize',10);
title('H_{E} H_{A} vs Temperature');
xlabel('Temperature (K)');
ylabel('H_{E} H_{A} (T^2)');
xticks(0:15:max(x));
xlim([min(x)-10 max(x)+10]);
exportgraphics(gca, "EAFieldVsTemperature_alt.png");

%% Testing
i = 1;
dmF = -2.26906719034998;
eaF = 0.112108347402597;
gamma = 28;
Fe2O3_exp.experimentArray(i).trialanderror(Filter.normalization, dmF, eaF, gamma);