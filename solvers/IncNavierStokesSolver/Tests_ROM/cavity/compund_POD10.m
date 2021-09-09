
% compute a compund POD from several trajectories which already underwent a POD (because they had so many time steps)

x_100 = load('100_triple/VCS_fields_TT_pod_x.txt');
x_200 = load('200_triple/VCS_fields_TT_pod_x.txt');
x_300 = load('300_triple/VCS_fields_TT_pod_x.txt');
x_400 = load('400_triple/VCS_fields_TT_pod_x.txt');
x_500 = load('500_triple/VCS_fields_TT_pod_x.txt');
x_600 = load('600_triple/VCS_fields_TT_pod_x.txt');
x_700 = load('700_triple/VCS_fields_TT_pod_x.txt');

y_100 = load('100_triple/VCS_fields_TT_pod_y.txt');
y_200 = load('200_triple/VCS_fields_TT_pod_y.txt');
y_300 = load('300_triple/VCS_fields_TT_pod_y.txt');
y_400 = load('400_triple/VCS_fields_TT_pod_y.txt');
y_500 = load('500_triple/VCS_fields_TT_pod_y.txt');
y_600 = load('600_triple/VCS_fields_TT_pod_y.txt');
y_700 = load('700_triple/VCS_fields_TT_pod_y.txt');


all_x = [ x_100' x_200' x_300' x_400' x_500' x_600' x_700'];
%all_x = [ x_120'  x_140' ];


[U,S,V] = svd(all_x);
dim_x = size(all_x, 2);

Ut = U';

sing_vals = diag(S(1:dim_x,1:dim_x))';
modes_taken_x = sum(cumsum(sing_vals) ./ sum(sing_vals) < 0.9999);

modes_taken_x
dim_x

writematrix(Ut(1:modes_taken_x,:), 'all_x.txt', 'Delimiter', ' ')
writematrix(modes_taken_x, 'dim_x.txt');


all_y = [ y_100' y_200' y_300' y_400' y_500' y_600' y_700'];
%all_y = [ y_120'  y_140' ];

[U,S,V] = svd(all_y);
dim_y = size(all_y, 2);
sing_vals = diag(S(1:dim_y,1:dim_y))';
modes_taken_y = sum(cumsum(sing_vals) ./ sum(sing_vals) < 0.9999);

modes_taken_y
dim_y
Ut = U';

writematrix(Ut(1:modes_taken_y,:), 'all_y.txt', 'Delimiter', ' ')
writematrix(modes_taken_y, 'dim_y.txt');


