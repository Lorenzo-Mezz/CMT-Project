% read CSV file and clean data
co2_data = readtable('../data/co2_annmean_mlo.csv'); 
co2_data.x______________________________________________________________ = [];
co2_data(1:29,:) = [];
c = split(co2_data.Var1, ',');
c(:,3) = num2cell(zeros(size(c, 1), 1));
year = str2double(c(:,1));
co2_cc = str2double(c(:,2));

% create clean table
c_table = table(year, co2_cc, 'VariableNames', {'Annee', 'CO2'});

% fit a 2nd-degree polynomial (quadratic trend) to historical CO2 data
model = polyfit(year, co2_cc, 2);
disp(model)

% time vector and generating monthly intervals
t = []
future_years = (2025:2050)';
t_min = min(year);  
for j = 1: length(future_years)
    for i = 1:12
        t = [t; future_years(j) + (i-1)/12];
    end
end

% predict future CO2 concentrations using the polynomial model and store in table
co2_pred = polyval(model, t);
co2_future = table(t, co2_pred, 'VariableNames', {'Annee','CO2'});

% combine historical data with future predictions and plot
c = [c_table; co2_future];
figure; 
plot(c.Annee, c.CO2)

xlabel('Time');
ylabel('Atmospheric CO2 (ppm)');
title('Atmospheric CO2 Prediction');

saveas(gcf, "../data/co2_pred.png");

