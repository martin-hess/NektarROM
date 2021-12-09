

l2_r0 = load('ROM_cluster_L2reduc0.txt');
linf_r0 = load('ROM_cluster_Linfreduc0.txt');

l2_r1 = load('ROM_cluster_L2reduc1.txt');
linf_r1 = load('ROM_cluster_Linfreduc1.txt');

l2_r2 = load('ROM_cluster_L2reduc2.txt');
linf_r2 = load('ROM_cluster_Linfreduc2.txt');

l2_r3 = load('ROM_cluster_L2reduc3.txt');
linf_r3 = load('ROM_cluster_Linfreduc3.txt');

l2_r4 = load('ROM_cluster_L2reduc4.txt');
linf_r4 = load('ROM_cluster_Linfreduc4.txt');

l2_r5 = load('ROM_cluster_L2reduc5.txt');
linf_r5 = load('ROM_cluster_Linfreduc5.txt');

l2_r6 = load('ROM_cluster_L2reduc6.txt');
linf_r6 = load('ROM_cluster_Linfreduc6.txt');

l2_r7 = load('ROM_cluster_L2reduc7.txt');
linf_r7 = load('ROM_cluster_Linfreduc7.txt');

l2_r8 = load('ROM_cluster_L2reduc8.txt');
linf_r8 = load('ROM_cluster_Linfreduc8.txt');

l2_r9 = load('ROM_cluster_L2reduc9.txt');
linf_r9 = load('ROM_cluster_Linfreduc9.txt');

l2_r10 = load('ROM_cluster_L2reduc10.txt');
linf_r10 = load('ROM_cluster_Linfreduc10.txt');

l2_r11 = load('ROM_cluster_L2reduc11.txt');
linf_r11 = load('ROM_cluster_Linfreduc11.txt');

l2_r12 = load('ROM_cluster_L2reduc12.txt');
linf_r12 = load('ROM_cluster_Linfreduc12.txt');

l2_conv_mean = [   mean(l2_r1) mean(l2_r2) mean(l2_r3) mean(l2_r4) mean(l2_r5) mean(l2_r6) mean(l2_r7) mean(l2_r8) mean(l2_r9) mean(l2_r10) mean(l2_r11) mean(l2_r12)];
l2_conv_max = [  max(l2_r1) max(l2_r2) max(l2_r3) max(l2_r4) max(l2_r5) max(l2_r6) max(l2_r7) max(l2_r8) max(l2_r9) max(l2_r10) max(l2_r11) max(l2_r12)];

linf_conv_mean = [ mean(linf_r1) mean(linf_r2) mean(linf_r3) mean(linf_r4) mean(linf_r5) mean(linf_r6) mean(linf_r7) mean(linf_r8) mean(linf_r9) mean(linf_r10) mean(linf_r11) mean(linf_r12)];
linf_conv_max = [ max(linf_r1) max(linf_r2) max(linf_r3) max(linf_r4) max(linf_r5) max(linf_r6) max(linf_r7) max(linf_r8) max(linf_r9) max(linf_r10) max(linf_r11) max(linf_r12)];

l2_conv_mean = flip(l2_conv_mean);
l2_conv_max = flip(l2_conv_max);




figure
plot((l2_conv_mean))
hold on
plot((l2_conv_max))



%plot(l2_adv(1:1000))
%plot(linf_adv(1:1000))




