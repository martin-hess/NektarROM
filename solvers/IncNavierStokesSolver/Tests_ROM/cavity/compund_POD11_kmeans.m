
% compute a compund POD from several trajectories which already underwent a POD (because they had so many time steps)

x_100 = load('100_triple/VCS_fields_TT_pod_x.txt');
x_200 = load('200_triple/VCS_fields_TT_pod_x.txt');
x_300 = load('300_triple/VCS_fields_TT_pod_x.txt');
x_400 = load('400_triple/VCS_fields_TT_pod_x.txt');
x_500 = load('500_triple/VCS_fields_TT_pod_x.txt');
x_600 = load('600_triple/VCS_fields_TT_pod_x.txt');
x_700 = load('700_triple/VCS_fields_TT_pod_x.txt');


all_x = [  x_100' x_200'  x_300'  x_400' x_500' x_600' x_700' ];

%all_x = [  x_100' mean(x_200)'  mean(x_300)'  mean(x_400)' mean(x_500)' mean(x_600)' mean(x_700)' ];


size(all_x)

idx = kmeans( all_x', 5 );


idx'

