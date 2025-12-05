co2_data = readtable('co2_annmean_mlo.csv');
co2_data.x______________________________________________________________ = [];
co2_data(1:29,:) = [];
c = split(co2_data.Var1, ',');
c(:,3) = num2cell(zeros(size(c, 1), 1));

year = str2double(c(:,1));
co2_cc = str2double(c(:,2));
c_table = table(year, co2_cc, 'VariableNames', {'Annee', 'CO2'});

model = polyfit(year, co2_cc, 2);
disp(model)

t = []
future_years = (2025:2050)';
t_min = min(year);  
for j = 1: length(future_years)
    for i = 1:12
        t = [t; future_years(j) + (i-1)/12];
    end
end

co2_pred = polyval(model, t)

mean_pred = polyval(model, t);
co2_future = table(t, co2_pred, 'VariableNames', {'Annee','CO2'});
disp(co2_future)
c = [c_table; co2_future];

figure; 
plot(c.Annee, c.CO2)
