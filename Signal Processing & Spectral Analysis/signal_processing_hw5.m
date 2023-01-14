clc; clear all; close all;

xlsFile = 'HW-5-Pre-event Data of No-110103 NCREE.xlsx';
data = xlsread(xlsFile);
SR = 200;
T = 1/SR;
L = length(data);

%% zero mean
for i= 1:14
    data(:,i) = data(:,i) - mean(data(:,i));
end

%% Data Pre-Processing (Low pass filter and then down sampling the data)
cff = 20;
Acc = zeros(L, 14);
[Para_B, Para_A] = butter(14, cff/(0.5*SR), 'low');
for i = 1:size(data, 2)
    Acc(:, i) = filtfilt(Para_B, Para_A, data(:, i));
end
clear Para_A Para_B

ratio = 4;
SR = SR/ratio;
data = Acc(1:ratio:end,:);

%% zero mean
for i= 1:14
    data(:,i) = data(:,i) - mean(data(:,i));
end

%% correlation
for i = 1:14
    for j = 1:14
        R_acc(i,j,:) = xcorr(data(:,i),data(:,j),'unbiased');
    end
end
pt = size(R_acc,3);
R_acc = R_acc(:,:,1:round(pt/2));
pt = size(R_acc,3); 
R_acc = R_acc(:,:,end:-1:1);
R_acc = R_acc(:,:,1:round(pt/2));
pt = size(R_acc,3);


%% power spectrum density
PSD = zeros(14,14,pt);
for i = 1:14
    for j = 1:14
        Lu(i,j,:) = (R_acc(i,j,:) + R_acc(j,i,:))/2;
        Qu(i,j,:) = (R_acc(i,j,:) - R_acc(j,i,:))/2;
        Lf(i,j,:) = real(fft(Lu(i,j,:)));
        Qf(i,j,:) = imag(fft(Qu(i,j,:)));
    end
end
PSD = Lf - Qf * 1i;

%% Eigen Value Spectrum(singular value distribution)
for i = 1:size(PSD,3)
    [U(:,:,i),SS,V(:,:,i)] = svd(PSD(:,:,i).');
    S(:,i) = diag(SS);
    U1(:,i) = U(:,1,i);
    PSD(:,:,i) = (U(:,:,i)*SS*U(:,:,i)').';
end
f = (0:pt-1) / (pt/SR);

%% plot singular value spectrum
figure
semilogy(f,S(1,:));
hold on
semilogy(f,S(2,:)/5)
semilogy(f,S(3,:)/10)
xlim([0,25])

%% Extract indexes for each principle frequency
for i = 1:length(f)
    if f(i)>1.464 && f(i)<1.465
        m1 = i;
    elseif f(i)>3.417 && f(i)<3.418
        m2 = i;
    elseif f(i)>3.808 && f(i)<3.809
        m3 = i;
    elseif f(i)>5.273 && f(i)<5.274
        m4 = i;
    elseif f(i)>5.566 && f(i)<5.567
        m5 = i;
    elseif f(i)>5.859 && f(i)<5.86
        m6 = i;
    elseif f(i)>7.421 && f(i)<7.422
        m7 = i;
    elseif f(i)>10.156 && f(i)<10.157
        m8 = i;
    elseif f(i)>10.351 && f(i)<10.352
        m9 = i;
    elseif f(i)>13.085 && f(i)<13.086
        m10 = i;
    elseif f(i)>13.378 && f(i)<13.379
        m11 = i;
    elseif f(i)>15.136 && f(i)<15.137
        m12 = i;
    elseif f(i)>17.480 && f(i)<17.481
        m13 = i;
    elseif f(i)>20.019 && f(i)<20.02
        m14 = i;
    end
end

%% mode shape
mode=zeros(14,14);
% fill mode shape array using the real values at the identified indices
for i=1:14
    mode(i,:) = [U1(i,m1),U1(i,m2),U1(i,m3),U1(i,m4),U1(i,m5),U1(i,m6),U1(i,m7),U1(i,m8),U1(i,m9),U1(i,m10),U1(i,m11),U1(i,m12),U1(i,m13),U1(i,m14)];
end
% floor number
floorplot=(0:1:13)';
% vertical line
vertical = zeros(14,1);

% plot mode shape
figure
plot(vertical,floorplot);
hold on
plot(mode(:,1),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 1.46484 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,2),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 3.41797 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,3),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 3.80859 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,4),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 5.27344 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,5),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 5.56641 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,6),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 5.85938 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,7),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 7.42188 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,8),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 10.1562 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,9),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 10.3516 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,10),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 13.0859 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,11),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 13.3789 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,12),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 15.1367 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,13),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 17.4805 Hz')
hold off

figure
plot(vertical,floorplot);
hold on
plot(mode(:,14),floorplot)
xlim([-1,1])
ylabel('Floor Number')
xlabel('Amplitude')
title('Mode shape for 20.0195 Hz')
hold off