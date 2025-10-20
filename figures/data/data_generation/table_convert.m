load('table.mat');
data = [Nvec', ovec', evec'];
writematrix(data, '../harmonic_table.csv');