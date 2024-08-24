% Set the column number: Frequency, Field, S12, Index
fp = FileParameters(2, 1, 4, 0);

% Read experiment from file
Fe2O3_exp = ExperimentCollection(["Fe2O3_0.dat", "Fe2O3_15.dat", ...
    "Fe2O3_30.dat", "Fe2O3_45.dat", "Fe2O3_60.dat", "Fe2O3_75.dat",...
    "Fe2O3_90.dat"], fp);

% Set units of the data in file
Fe2O3_exp.setfieldunits(MagneticUnit.tesla);
Fe2O3_exp.setfrequencyunits(FrequencyUnit.hertz);

% Convert units
Fe2O3_exp.convertfieldunits(MagneticUnit.tesla);
Fe2O3_exp.convertfrequencyunits(Scale.giga * FrequencyUnit.hertz);

% Filter by positive or negative field
%   If you want the positive field, comment the 
%   first line and uncomment the second line
Fe2O3_exp.filternegativefield();
% Fe2O3_exp.filterpositivefield();

%% Make image plot
% When the plot pops up, press ENTER to be able to
% zoom in, and modify the plot options
for i = 1:7
    % There are 3 methods of filtering to produce the image plots: (ranked
    % from the least smooth to smoothest method)
    % 1. Filter.none (without filter)
    % 2. Filter.normalization (using normalization / Marc method)
    % 3. Filter.noiseCancelling (Pablo method / ask Pablo to know more)
    Fe2O3_exp.experimentArray(i).getcrossingpoints(Filter.noiseCancelling);
    clim([-0.1 0]); % If using noise cancelling
    xlabel("$\mu_0 H$ (T)", 'Interpreter','latex');
    ylabel("Frequency (GHz)");
    exportgraphics(gca, "imageplot\figure_hef_"+i+".png");
    % exportgraphics(gca, "imageplot_lf\figure_lef_"+i+".png");
end

%% Plot S12 vs Field for each frequency and select peaks
% To go to the next frequency select a range for the peak,
% press ENTER or close the window
for i = 1:6
    Fe2O3_exp.resonanceanalysis(i);
end

%% Make kittel plot
% Magnetization can be set to inPlane, outOfPlane, paramagnetic, aFM
for i = 1:6
    % Change Magnetization.aFM to Magnetization.inPlane for obtaining M_eff instead of H_DM
    Fe2O3_exp.experimentArray(i).setmagnetization(Magnetization.aFM);
    Fe2O3_exp.experimentArray(i).makekittel();
    annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$H_{DM}$ = ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.dmField.value), ...
        " $\pm$ ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.dmField.error), " ", ...
        Fe2O3_exp.experimentArray(1, i).kittelParameters.dmField.unit.tag) strcat("$H_{E} H_{A}$ = ", ...
        num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.eaField.value), " $\pm$ ", ... 
        num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.eaField.error), " ", ...
        Fe2O3_exp.experimentArray(1, i).kittelParameters.eaField.unit.tag)], 'Interpreter','latex');
    % exportgraphics(gca, "kittelplot\kittel_hef_3_"+i+".png");
    exportgraphics(gca, "kittelplot_lf\kittel_lef_"+i+".png");
    % end
end

%% Make damping plot
for i = 1:6
    Fe2O3_exp.experimentArray(i).makedampingfit();
    annotation('textbox', [0.5, 0.15, 0.1, 0.1], 'String', ...
        [strcat("$\alpha$ = ", num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.damping.value), ...
        " $\pm$ ", num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.damping.error)) strcat("$\alpha_{inhomogenous}$ = ", ...
        num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.inhomogeneousDamping.value), " $\pm$ ", ... 
        num2str(Fe2O3_exp.experimentArray(1, i).dampingParameters.inhomogeneousDamping.error))], 'Interpreter','latex');
    % exportgraphics(gca, "dampingplot\damping_hef_2_"+i+".png");
    exportgraphics(gca, "dampingplot_lf\damping_lef_"+i+".png");
end

%% vs degree
x = [0, 15, 30, 45, 60, 75];
extracted_data = zeros(6,9);
j=1;
for i = 1:6
    % if (i~=3) && (i~=7) && (i~=8)
        extracted_data(j,:) = [x(i), Fe2O3_exp.experimentArray(1,i).kittelParameters.dmField.value, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.dmField.error, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.damping.value, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.damping.error, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.eaField.value, ...
            Fe2O3_exp.experimentArray(1,i).kittelParameters.eaField.error, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.inhomogeneousDamping.value, ...
            Fe2O3_exp.experimentArray(1,i).dampingParameters.inhomogeneousDamping.error];
        j = j+1;
    % end
end

%% DM field vs plane orientation
figure( 'Name', 'DM field vs plane orientation' );
errorbar(extracted_data(:,1),extracted_data(:,2),extracted_data(:,3),'.','MarkerSize',10);
title('H_{DM} vs Plane Orientation');
xlabel('Orientation (degree)');
ylabel('H_{DM} (T)');
xticks(0:15:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "DMFieldVsOrientation_hef_2.png");
% exportgraphics(gca, "DMFieldVsOrientation_lef.png");

%% Damping vs plane orientation
figure( 'Name', 'Damping vs plane orientation' );
errorbar(extracted_data(:,1),extracted_data(:,4),extracted_data(:,5),'.','MarkerSize',15);
title('\alpha vs Plane Orientation');
xlabel('Orientation (degree)');
ylabel('\alpha');
xticks(0:15:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "DampingVsOrientation_hef_2.png");
% exportgraphics(gca, "DampingVsOrientation_lef.png");

%% EA field vs plane orientation
figure( 'Name', 'EA field vs plane orientation' );
errorbar(extracted_data(:,1),extracted_data(:,6),extracted_data(:,7),'.','MarkerSize',10);
title('H_{E} H_{A} vs Plane Orientation');
xlabel('Orientation (degree)');
ylabel('H_{E} H_{A} (T^2)');
xticks(0:15:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "EAFieldVsOrientation_hef_2.png");
% exportgraphics(gca, "EAFieldVsOrientation_lef.png");

%% Inhomogeneous Damping vs plane orientation
figure( 'Name', 'Inhomogeneous damping vs plane orientation' );
errorbar(extracted_data(:,1),extracted_data(:,8),extracted_data(:,9),'.','MarkerSize',15);
title('Inhomogeneous damping vs Plane Orientation');
xlabel('Orientation (degree)');
ylabel('Inhomogeneous damping');
xticks(0:15:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
exportgraphics(gca, "InhonogeneousDampingVsOrientation_hef_2.png");
% exportgraphics(gca, "InhonogeneousDampingVsOrientation_lef.png");

%% Make kittel plot (alt to get Me_ff instead of H_DM)
% % Magnetization can be set to inPlane, outOfPlane, paramagnetic, aFM
% for i = 1:7
%     Fe2O3_exp.experimentArray(i).setmagnetization(Magnetization.inPlane);
%     Fe2O3_exp.experimentArray(i).makekittel();
%     annotation('textbox', [0.15, 0.15, 0.1, 0.1], 'String', ...
%         [strcat("$M_{eff}$ = ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.value), ...
%         " $\pm$ ", num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.error), " ", ...
%         Fe2O3_exp.experimentArray(1, i).kittelParameters.effectiveMagnetization.unit.tag) strcat("$H_{A}$ = ", ...
%         num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.anisotropyField.value), " $\pm$ ", ... 
%         num2str(Fe2O3_exp.experimentArray(1, i).kittelParameters.anisotropyField.error), " ", ...
%         Fe2O3_exp.experimentArray(1, i).kittelParameters.anisotropyField.unit.tag)], 'Interpreter','latex');
%     exportgraphics(gca, "kittelplot_alt\kittel_hef_alt_"+i+".png");
%     % exportgraphics(gca, "kittelplot_alt\kittel_lef_alt_"+i+".png");
%     % end
% end

% %% vs temperature (alt)
% x = [30, 40, 57, 75, 92, 109, 126, 143, 185, 220, 298];
% extracted_data_alt = zeros(9,5);
% j=1;
% for i = 1:11
%     if (i~=3) && (i~=7) && (i~=8)
%         extracted_data_alt(j,:) = [x(i), Fe2O3_exp.experimentArray(1,i).kittelParameters.effectiveMagnetization.value, ...
%             Fe2O3_exp.experimentArray(1,i).kittelParameters.effectiveMagnetization.error, ...
%             Fe2O3_exp.experimentArray(1,i).kittelParameters.anisotropyField.value, ...
%             Fe2O3_exp.experimentArray(1,i).kittelParameters.anisotropyField.error];
%         j = j+1;
%     end
% end
% 
% %% Meff vs temperature
% figure( 'Name', 'Meff vs temperature' );
% errorbar(extracted_data_alt(:,1),extracted_data_alt(:,2),extracted_data_alt(:,3),'.','MarkerSize',25);
% title('M_{eff} vs temperature');
% xlabel('Temperature (K)');
% ylabel('M_{eff} (T)');
% xticks(0:20:max(extracted_data_alt(:,1)));
% xlim([min(extracted_data_alt(:,1))-10 max(extracted_data_alt(:,1))+10]);
% exportgraphics(gca, "MeffVsTemperature.png");
% 
% %% H ani vs temperature
% figure( 'Name', 'Hani vs temperature' );
% errorbar(extracted_data_alt(:,1),extracted_data_alt(:,4),extracted_data_alt(:,5),'.','MarkerSize',10);
% title('H_{ani} vs temperature');
% xlabel('Temperature (K)');
% ylabel('H_{ani} (T)');
% xticks(0:20:max(extracted_data_alt(:,1)));
% xlim([min(extracted_data_alt(:,1))-10 max(extracted_data_alt(:,1))+10]);
% exportgraphics(gca, "HaniVsTemperature.png");

%% Marking points for the second method
for i = 1:7
    Fe2O3_exp.experimentArray(i).getcrossingpoints(Filter.noiseCancelling);
end

%% Alternative method
x = [0, 15, 30, 45, 60, 75, 90];
save = zeros(7,4);
for i = 1:7
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
    exportgraphics(gca, "kittelplot\kittel_m2_2_"+i+".png");
    % exportgraphics(gca, "kittelplot_lf\kittel_lef_m2_"+i+".png");

    Fe2O3_exp.experimentArray(i).trialanderror(Filter.noiseCancelling, ...
        fitresult.dmF, fitresult.eaF, fitresult.gamma);
    clim([-0.1 0]); % If using noise cancelling
    exportgraphics(gca, "kittelplot\kittel_m2_2dimg_2_"+i+".png");
    % exportgraphics(gca, "kittelplot_lf\kittel_lef_m2_2dimg_"+i+".png");
    save(i,1) = fitresult.dmF;
    save(i,2) = dmConfInt;
    save(i,3) = fitresult.eaF;
    save(i,4) = eaConfInt;
end

%% DM field vs plane orientation (alt)
figure( 'Name', 'DM field vs plane orientation' );
errorbar(x,save(:,1),save(:,2),'.','MarkerSize',10);
title('H_{DM} vs Plane Orientation');
xlabel('Orientation (degree)');
ylabel('H_{DM} (T)');
xticks(0:15:max(x));
xlim([min(x)-10 max(x)+10]);
exportgraphics(gca, "DMFieldVsOrientation_hef_alt_2.png");
% exportgraphics(gca, "DMFieldVsOrientation_lef_alt.png");

%% EA field vs plane orientation (alt)
figure( 'Name', 'EA field vs plane orientation' );
errorbar(x,save(:,3),save(:,4),'.','MarkerSize',10);
title('H_{E} H_{A} vs Plane Orientation');
xlabel('Orientation (degree)');
ylabel('H_{E} H_{A} (T^2)');
xticks(0:15:max(x));
xlim([min(x)-10 max(x)+10]);
% exportgraphics(gca, "EAFieldVsOrientation_hef_alt_2.png");
exportgraphics(gca, "EAFieldVsOrientation_lef_alt.png");

%% Testing
i = 1;
dmF = -2; %-1.8714;
eaF = 0.016263;
gamma = 28;
Fe2O3_exp.experimentArray(i).trialanderror(Filter.noiseCancelling, dmF, eaF, gamma);

%% Method comparison
figure( 'Name', 'Method comparison for DM field');
errorbar(extracted_data(:,1),extracted_data(:,2),extracted_data(:,3), ...
    'r.','MarkerSize',10,'DisplayName','Method 1');
title('Method comparison for H_{DM}');
xlabel('Orientation (degree)');
ylabel('H_{DM} (T)');
xticks(0:15:max(extracted_data(:,1)));
xlim([min(extracted_data(:,1))-10 max(extracted_data(:,1))+10]);
hold on;
errorbar(x,save(:,1),save(:,2),'b.','MarkerSize',10,'DisplayName','Method 2');
legend('show','Location','southeast');
hold off;
exportgraphics(gca, "MethodComp_dmF_hef_2.png");
% exportgraphics(gca, "MethodComp_dmF_lef.png");

%% Approximation for method 1
figure('Name','Approximation for method 1')
title('Approximation for method 1');
xlabel('$\mu_0 H$ (T)', 'Interpreter', 'latex');
ylabel('Frequency (GHz)', 'Interpreter', 'latex');
hold on
for i=1:7
plot(Fe2O3_exp.experimentArray(1, i).resonanceParameters(:,2)*cos(x(i)*pi/180), ...
    Fe2O3_exp.experimentArray(1, i).resonanceParameters(:,1),"x",'LineWidth',1.5, ...
    'DisplayName', sprintf('%i%c',x(i),char(176)));
end
legend('show','Location','southwest');
hold off
exportgraphics(gca, "Approx_m1_hef.png");
% exportgraphics(gca, "Approx_m1_lef.png");

%% Approximation for method 2
figure('Name','Approximation for method 2')
title('Approximation for method 2');
xlabel('$\mu_0 H$ (T)', 'Interpreter', 'latex');
ylabel('Frequency (GHz)', 'Interpreter', 'latex');
hold on
for i=1:7
plot(Fe2O3_exp.experimentArray(1, i).crossingPoints(:,1)*cos(x(i)*pi/180), ...
    Fe2O3_exp.experimentArray(1, i).crossingPoints(:,2),"x",'LineWidth',1.5, ...
    'DisplayName', sprintf('%i%c',x(i),char(176)));
end
legend('show','Location','southwest');
hold off
exportgraphics(gca, "Approx_m2_hef.png");
% exportgraphics(gca, "Approx_m2_lef.png");
