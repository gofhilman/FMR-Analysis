figure()
hold on
for i=1:7
plot(Fe2O3_exp.experimentArray(1, i).resonanceParameters(:,2)*cos(x(i)*pi/180), ...
    Fe2O3_exp.experimentArray(1, i).resonanceParameters(:,1),"x",'LineWidth',2, ...
    'DisplayName', sprintf('%i%c',x(i),char(176)));
end
legend('show','Location','southwest');
hold off