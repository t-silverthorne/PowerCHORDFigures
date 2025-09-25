fig = openfig('test_Nmeas48_Amp10_MaxIter100_fmin1_fmax24_Nfreqch64_Nacroch64_Nsampch100_NfqTinf64_NfqT25000_Nsampmc500_Nfreqmc64_Nacromc64_Npermmc1000.fig', ...
    'invisible');
ax = findall(fig,'type','axes');
ax_sorted = flipud(ax);
top_ax = ax_sorted(3);        % top panel
dots = findobj(top_ax,'Type','Line','Marker','.');
x = get(dots,'XData');
y = get(dots,'YData');
x = x(:)
vec = cell2mat(x(:));
