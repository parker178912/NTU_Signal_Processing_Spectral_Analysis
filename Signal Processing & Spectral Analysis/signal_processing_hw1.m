clc;clear;close all;
% (a)
subplot(3, 1, 1);
fs1 = 10;
N1 = 256;
Ndata_a = 256;
n1 = 0:Ndata_a-1;
dt1 = n1/fs1;
x1 = sin(2*pi*dt1) + sin(2*pi*1.02*dt1) + sin(2*pi*1.07*dt1);
y1 = fft(x1, N1);
mag1 = abs(y1);
ang1 = angle(y1);
f1 = (0:N1-1)*fs1/N1;
plot(f1(1:N1/2), mag1(1:N1/2)*2/N1);
xlabel('Hz'); ylabel('AMP');
title('Ndata = 256, NFFT = 256, fs = 10');
xlim([0 2]);

% (b)
subplot(3, 1, 2);
fs2 = 10;
N2 = 2048;
Ndata_b = 256;
n2 = 0:Ndata_b-1;
dt2 = n2/fs2;
x2 = sin(2*pi*dt2) + sin(2*pi*1.02*dt2) + sin(2*pi*1.07*dt2);
y2 = fft(x2, N2);
mag2 = abs(y2);
ang2 = angle(y2);
f2 = (0:N2-1)*fs2/N2;
plot(f2(1:N2/2), mag2(1:N2/2)*2/N2);
xlabel('Hz'); ylabel('AMP');
title('Ndata = 256, NFFT = 2048, fs = 10');
xlim([0 2]);

% (c)
subplot(3, 1, 3);
fs3 = 100;
N3 = 2048;
Ndata_c= 256;
n3 = 0:Ndata_c-1;
dt3 = n3/fs3;
x3 = sin(2*pi*dt3) + sin(2*pi*1.02*dt3) + sin(2*pi*1.07*dt3);
y3 = fft(x3, N2);
mag3 = abs(y3);
ang3 = angle(y3);
f3 = (0:N3-1)*fs3/N3;
plot(f3(1:N3/2), mag3(1:N3/2)*2/N3);
xlabel('Hz'); ylabel('AMP');
title('Ndata = 256, NFFT = 2048, fs = 100');
xlim([0 2]);
