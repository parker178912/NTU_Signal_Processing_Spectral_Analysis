clc; clear all; close all;
%% Load Data
% Case 1
a11 = importdata('./case- 1/WN_50gal _SPECIMEN1_65_A1a__.txt');
a21 = importdata('./case- 1/WN_50gal _SPECIMEN1_67_A2a__.txt');
a31 = importdata('./case- 1/WN_50gal _SPECIMEN1_69_A3a__.txt');
a41 = importdata('./case- 1/WN_50gal _SPECIMEN1_71_A4a__.txt');
a51 = importdata('./case- 1/WN_50gal _SPECIMEN1_73_A5a__.txt');
a61 = importdata('./case- 1/WN_50gal _SPECIMEN1_75_A6a__.txt');

b11 = importdata('./case- 1/WN_50gal _SPECIMEN1_66_A1b__.txt');
b21 = importdata('./case- 1/WN_50gal _SPECIMEN1_68_A2b__.txt');
b31 = importdata('./case- 1/WN_50gal _SPECIMEN1_70_A3b__.txt');
b41 = importdata('./case- 1/WN_50gal _SPECIMEN1_72_A4b__.txt');
b51 = importdata('./case- 1/WN_50gal _SPECIMEN1_74_A5b__.txt');
b61 = importdata('./case- 1/WN_50gal _SPECIMEN1_76_A6b__.txt');

data_a1 = [a11 a21 a31 a41 a51 a61];
data_b1 = [b11 b21 b31 b41 b51 b61];
average1 = (data_a1 + data_b1)./2;
% Case 8
a18 = importdata('./case- 8/WN_50gal _SPECIMEN1_65_A1a__.txt');
a28 = importdata('./case- 8/WN_50gal _SPECIMEN1_67_A2a__.txt');
a38 = importdata('./case- 8/WN_50gal _SPECIMEN1_69_A3a__.txt');
a48 = importdata('./case- 8/WN_50gal _SPECIMEN1_71_A4a__.txt');
a58 = importdata('./case- 8/WN_50gal _SPECIMEN1_73_A5a__.txt');
a68 = importdata('./case- 8/WN_50gal _SPECIMEN1_75_A6a__.txt');

b18 = importdata('./case- 8/WN_50gal _SPECIMEN1_66_A1b__.txt');
b28 = importdata('./case- 8/WN_50gal _SPECIMEN1_68_A2b__.txt');
b38 = importdata('./case- 8/WN_50gal _SPECIMEN1_70_A3b__.txt');
b48 = importdata('./case- 8/WN_50gal _SPECIMEN1_72_A4b__.txt');
b58 = importdata('./case- 8/WN_50gal _SPECIMEN1_74_A5b__.txt');
b68 = importdata('./case- 8/WN_50gal _SPECIMEN1_76_A6b__.txt');

data_a8 = [a18 a28 a38 a48 a58 a68];
data_b8 = [b18 b28 b38 b48 b58 b68];
average8 = (data_a8 + data_b8)./2;
%% System Parameter
L = length(a11);
Fs = 200; % (hz)
T = 1/Fs; % (s)
time = (0: L-1) *T;
%% Plot the data in time domain Case 1 & 8, Floor1~6
figure(1)
for i = 1:6
    subplot(6,1,i)
    plot(time(5001:6001), average1(5001:6001, i), 'g', time(5001:6001), average8(5001:6001, i), 'm')
    ylabel('Acceleration (g)')
    title( ['Floor ', num2str(i), ' in Time Domain, x axis (time in sec)' ] )
    ylim( [-0.15 0.15] )
    legend('Case 1', 'Case 8')
end
 %% Data Pre-Processing (Low pass filter and then down sampling the data)
% Case 1
cff = 20;
Acc1 = zeros(L, 6);
[Para_B, Para_A] = butter(6, cff/(0.5*Fs), 'low');
for I = 1:size(data_a1, 2)
    Acc1(:, I) = filtfilt(Para_B, Para_A, average1(:, I));
end
clear Para_A Para_B
% zero mean
for i = 1:6
    Acc1(:, i) = Acc1(:, i) - mean(average1(:, i));
end

% Case 8
cff = 20;
Acc8 = zeros(L, 6);
[Para_B, Para_A] = butter(6, cff/(0.5*Fs), 'low');
for I = 1:size(data_a8, 2)
    Acc8(:, I) = filtfilt(Para_B, Para_A, average8(:, I));
end
clear Para_A Para_B

% zero mean
for i = 1:6
    Acc8(:, i) = Acc8(:, i) - mean(average8(:, i));
end
%% Auto-spectral density function
n = L; % correlation data length
f = Fs*(0: (L-1)/ 2)' /L;

% create zero matrix after pre-processing
% Case 1
[~, ch1] = size(Acc1);
for i = 1: ch1
    for j = 1:ch1
        [Ryy1(i,j,:) ] = xcorr(Acc1(:,i), Acc1(:,j), 'unbiased');   
    end
end

pt1 = size(Ryy1, 3);
Ryy1 = Ryy1(:, :, 1:round(pt1/2));
pt1 = size(Ryy1, 3);
Ryy1 = Ryy1(: ,: , end: -1:1);
Ryy1 = Ryy1(:, :, 1:round(pt1/2));
pt1 = size(Ryy1, 3);
Freq1 = (0: pt1-1) / (pt1/Fs);
PSD1 = zeros(ch1, ch1, pt1);

for i = 1: ch1
    for j = 1:ch1
        Lu1(i, j, :) = (Ryy1(i, j, :) + Ryy1(j, i, :)) /2;
        Qu1(i, j, :) = (Ryy1(i, j, :) - Ryy1(j, i, :)) /2;   
        Lf1(i, j, :) = real(fft(Lu1(i, j, :)));
        Qf1(i, j, :) = imag(fft(Qu1(i, j, :)));
    end
end

PSD1 = Lf1 - (Qf1 * 1i);

% Case 8
[~, ch8] = size(Acc8);
for i = 1: ch8
    for j = 1:ch8
        [Ryy8(i,j,:) ] = xcorr(Acc8(:,i), Acc8(:,j), 'unbiased');   
    end
end

pt8 = size(Ryy8, 3);
Ryy8 = Ryy8(:, :, 1:round(pt8/2));
pt8 = size(Ryy8, 3);
Ryy8 = Ryy8(: ,: , end: -1:1);
Ryy8 = Ryy8(:, :, 1:round(pt8/2));
pt8 = size(Ryy8, 3);
Freq8 = (0: pt8-1) / (pt8/Fs);
PSD8 = zeros(ch8, ch8, pt8);

for i = 1: ch8
    for j = 1:ch8
        Lu8(i, j, :) = (Ryy8(i, j, :) + Ryy8(j, i, :)) /2;
        Qu8(i, j, :) = (Ryy8(i, j, :) - Ryy8(j, i, :)) /2;   
        Lf8(i, j, :) = real(fft(Lu8(i, j, :)));
        Qf8(i, j, :) = imag(fft(Qu8(i, j, :)));
    end
end

PSD8 = Lf8 - (Qf8 * 1i);
%% Single Value Decomposition
% Case1
L_PSD1 = size(PSD1, 3);

for i = 1: L_PSD1
    % perform singular value decomposition
    [ U1(:, :, i), SS1, V1(:, :, i) ] = svd(PSD1(:, :, i).' );
    % Extract Diagonal
    Sdiag1(:, i) = diag(SS1);
    % Extract Values
    U11(:, i) = U1(:, 1, i);
end

% plot SVD
figure()
semilogy(Freq1, Sdiag1(1, :), 'm', Freq1, Sdiag1(2, :)/5 , 'g',  Freq1, Sdiag1(3, :)/10, 'k')
xlim([0 20])
xlabel('Freqency (Hz)')
ylabel('Amplitude')
title( ['Singular Value Spectrum Case 1' ] )

% Case 8
L_PSD8 = size(PSD8, 3);

for i = 1: L_PSD8
    % perform singular value decomposition
    [ U8(:, :, i), SS8, V8(:, :, i) ] = svd(PSD8(:, :, i).' );
    % Extract Diagonal
    Sdiag8(:, i) = diag(SS8);
    % Extract Values
    U18(:, i) = U8(:, 1, i);
end

% plot SVD
figure()
semilogy(Freq8, Sdiag8(1, :), 'm', Freq1, Sdiag8(2, :)/5 , 'g',  Freq1, Sdiag8(3, :)/10, 'k')
xlim([0 20])
xlabel('Freqency (Hz)')
ylabel('Amplitude')
title( ['Singular Value Spectrum Case 8' ] )

%% Extract indexes for each principle frequency
% Case 1
for i = 1: size(PSD1,3)
    if Freq1(i) > 1.1267 && Freq1(i) < 1.1269
        M11 = i;
    elseif Freq1(i) > 3.6657 && Freq1(i) < 3.6659
        M21 = i;
    elseif Freq1(i) > 6.29506 && Freq1(i) < 6.29508
        M31 = i;
    elseif Freq1(i) > 9.1796 && Freq1(i) < 9.1797
        M41 = i;
    elseif Freq1(i) > 12.0943 && Freq1(i) < 12.0945
        M51 = i;
    elseif Freq1(i) > 14.2878 && Freq1(i) < 14.288
        M61 = i;
    end
end

% Case 8
for i = 1: size(PSD8,3)
    if Freq8(i) > 1.11177 && Freq8(i) < 1.11179
        M18 = i;
    elseif Freq8(i) > 3.54566 && Freq8(i) < 3.54568
        M28 = i;
    elseif Freq8(i) > 6.23497 && Freq8(i) < 6.23499
        M38 = i;
    elseif Freq8(i) > 9.11958 && Freq8(i) < 9.1196
        M48 = i;
    elseif Freq8(i) > 11.9290 && Freq8(i) < 11.9291
        M58 = i;
    elseif Freq8(i) > 14.2427 && Freq8(i) < 14.2429
        M68 = i;
    end
end
    
% Create empty mode shape array
fr1 = zeros(7, 6);
fr8 = zeros(7, 6);
    
% Fill mode shape array at the identified indices
for i = 1:6
    fr1(i+1, :) = [ U11(i, M11), U11(i, M21) ,U11(i, M31) ,U11(i, M41) ,U11(i, M51) ,U11(i, M61) ] ;
    fr8(i+1, :) = [ U18(i, M18), U18(i, M28) ,U18(i, M38) ,U18(i, M48) ,U18(i, M58) ,U18(i, M68) ] ;
end
    
 % create array of floor number
floorplot = (0:1:6)' ;
 % creat vertical line to plot against
vertical = [0,0,0,0,0,0,0];
 
% Radian and Angle
R1 = abs(fr1); R8 = abs(fr8);
theta1 = angle(fr1);theta8 = angle(fr8);
%%
% plot mode shape for the 1st frequency
figure()
subplot(2,1,1)
plot(vertical, floorplot, '-', fr1(:, 1), floorplot, '-mo', fr8(:, 1), floorplot, '-ko')
xlim([-1 1])
ylabel('Floor Number')
xlabel('Amplitude')
title( ['1st Mode Shape' ] )
legend('Floor plot', 'Case 1', 'Case 8')
% 1st mode location of poles
subplot(2,1,2)
polarplot(theta1(:, 1), R1(:, 1), 'mo', theta8(:, 1), R8(:, 1), 'k*');
rlim([0 1])
title('Location of poles')

 % plot mode shape for the 2nd frequency
figure()
subplot(2,1,1)
plot(vertical, floorplot, '-', fr1(:, 2), floorplot, '-mo', fr8(:, 2), floorplot, '-ko')
 xlim([-1 1])
ylabel('Floor Number')
xlabel('Amplitude')
title( ['2nd Mode Shape' ] )
legend('Floor plot', 'Case 1', 'Case 8')
% 2nd mode location of poles
subplot(2,1,2)
polarplot(theta1(:, 2), R1(:, 2), 'mo', theta8(:, 2), R8(:, 2), 'k*');
rlim([0 1])
title('Location of poles')
%
 % plot mode shape for the 3rd frequency
figure()
% subplot(2, 1, 1, 'position', [100,500,200,1000])
% pos1 = set(p1,'position',[100,0,500,1500])
subplot(2,1,1)
plot(vertical, floorplot, '-', fr1(:, 3), floorplot, '-mo', fr8(:, 3), floorplot, '-ko')
xlim([-1 1])
ylabel('Floor Number')
xlabel('Amplitude')
title( ['3rd Mode Shape' ] )
legend('Floor plot', 'Case 1', 'Case 8')
% 3nd mode location of poles
subplot(2, 1, 2)
% set(gcf,'position',[100,0,200,200])
polarplot(theta1(:, 3), R1(:, 3), 'mo', theta8(:, 3), R8(:, 3), 'k*');
rlim([0 1])
title('Location of poles')

% plot mode shape for the 4nd frequency
figure()
subplot(2,1,1)
plot(vertical, floorplot, '-', fr1(:, 4), floorplot, '-mo', fr8(:, 4), floorplot, '-ko')
xlim([-1 1])
ylabel('Floor Number')
xlabel('Amplitude')
title( ['4nd Mode Shape' ] )
legend('Floor plot', 'Case 1', 'Case 8')
% 4nd mode location of poles
subplot(2,1,2)
polarplot(theta1(:, 4), R1(:, 4), 'mo', theta8(:, 4), R8(:, 4), 'k*');
rlim([0 1])
title('Location of poles')

% plot mode shape for the 5nd frequency
figure()
subplot(2,1,1)
plot(vertical, floorplot, '-', fr1(:, 5), floorplot, '-mo', fr8(:, 5), floorplot, '-ko')
xlim([-1 1])
ylabel('Floor Number')
xlabel('Amplitude')
title( ['5nd Mode Shape' ] )
legend('Floor plot', 'Case 1', 'Case 8')
% 5nd mode location of poles
subplot(2,1,2)
polarplot(theta1(:, 5), R1(:, 5), 'mo', theta8(:, 5), R8(:, 5), 'k*');
rlim([0 1])
title('Location of poles')
    
% plot mode shape for the 6nd frequency
figure()
subplot(2,1,1)
plot(vertical, floorplot, '-', fr1(:, 6), floorplot, '-mo', fr8(:, 6), floorplot, '-ko')
xlim([-1 1])
ylabel('Floor Number')
xlabel('Amplitude')
title( ['6nd Mode Shape' ] )
legend('Floor plot', 'Case 1', 'Case 8')
% 6nd mode location of poles
subplot(2,1,2)
polarplot(theta1(:, 6), R1(:, 6), 'mo', theta8(:, 6), R8(:, 6), 'k*');
rlim([0 1])
title('Location of poles')