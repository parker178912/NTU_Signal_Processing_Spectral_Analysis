clc;clear;close all;
% property define
fs = 200; % sampling rate

%Only 10 seconds of data used for analysis (t=10 to t=20 seconds)
%Select 20 sec data for analysis
%the following 4 statement need to be modified
all_data = importdata('6-story ambient test structure with cut member-HW7.txt');
data = all_data.data;
pt_t15 = fs*15;
pt_t20 = fs*20;
t = data(pt_t15+1:pt_t20,1);
accel=data(pt_t15+1:pt_t20,2:end);

%Sample period
dt = 1/fs; %Unit: Sec

[pt,ch] = size(accel);

%Anticipated pairs of singular values
totalpair = 6;

%% Embedding
nrh = 100;
nch = pt-nrh+1;
X=zeros(ch*nrh,nch);
for i=1:ch
    X(i:ch:end,:) = hankel(accel(1:nrh,i),accel(nrh:pt,i));
end

%% Singular Value Decomposition
if true
    [U,S,V] = svd(X,'econ');
else
    Scov = X*X';
    [U,S,V] = svd(Scov,'econ');
    clear Scov
end
S = diag(S);

%% Grouping
figure
plot(1:length(S),S,'-o','LineWidth',1.5)
xlim([0,20]);
xlabel('Number of Singular Values')
ylabel('Singular Value')
title('Singular Value Plot')

figure
plot(1:length(S),100*cumsum(S)/sum(S),'-o','LineWidth',1.5);
xlim([0 20]);
xlabel('Number of Singular Value')
ylabel('Percentage of Accumulated Singular Values')
title('Singular Value Distribution, L=100')

%Based on the plot, select pairs of SV to reconstruct data
X_recon = zeros([size(X),totalpair]);
for I = 1:totalpair
    tmpidx = [2*I-1,2*I]; %Creates SV pairs (1,2; 3,4; 5,6; 7,8; 9,10; 11,12)
    if true
        tmp = U(:,tmpidx)*diag(S(tmpidx))*V(:,tmpidx)'; %Pulls USV values associated with pair
    else
        tempv1 = X'*U(:,tmpidx(1))/sqrt(S(tmpidx(1)));
        tempv2 = X'*U(:,tmpidx(2))/sqrt(S(tmpidx(2)));
        tmp = U(:,tmpidx)*diag(Sqrt(S(tmpidx)))*[tempv1,tempv2];
    end
    X_recon(:,:,I)=tmp; %Places matrix generated by U*S*V' for each pair into 3D matrix.
end

clear tmpidx tmp tmpv1 tmpv2

%% Step 4: Reconstruction

x_recon = zeros(pt,ch,totalpair);
for K = 1:totalpair  %Loop through for each pair of SVs
    for I = 1:ch  %Loop through for each channel of data
        %Pull individual Hankel matrix associated with each channel
        tmph = X_recon(I:ch:end,:,K);
        tmph = flipud(tmph);
        
        %Reconstruct measurements
        for J = 1:pt  %Loop through for each measurement point
            tmp=diag(tmph,J-nrh);
            x_recon(J,I,K)=mean(tmp);
        end
    end
end

clear tmph

%Fully reconstructed signals generated by summing the SV pair results
ReCon_1 = x_recon(:,1,1) + x_recon(:,1,2) + x_recon(:,1,3) + x_recon(:,1,4) + x_recon(:,1,5)+x_recon(:,1,6);
ReCon_2 = x_recon(:,2,1) + x_recon(:,2,2) + x_recon(:,2,3) + x_recon(:,2,4) + x_recon(:,2,5)+x_recon(:,2,6);
ReCon_3 = x_recon(:,3,1) + x_recon(:,3,2) + x_recon(:,3,3) + x_recon(:,3,4) + x_recon(:,3,5)+x_recon(:,3,6);
ReCon_4 = x_recon(:,4,1) + x_recon(:,4,2) + x_recon(:,4,3) + x_recon(:,4,4) + x_recon(:,4,5)+x_recon(:,4,6);
ReCon_5 = x_recon(:,5,1) + x_recon(:,5,2) + x_recon(:,5,3) + x_recon(:,5,4) + x_recon(:,5,5)+x_recon(:,5,6);
ReCon_6 = x_recon(:,6,1) + x_recon(:,6,2) + x_recon(:,6,3) + x_recon(:,6,4) + x_recon(:,6,5)+x_recon(:,6,6);



%% Plot Reconstructed Signals

%Show all SV pair contribution to reconstructing the original signal for
%each measurement channel

%Channel 1 (Floor 1) components
figure 
subplot(7,1,1)
plot(t,accel(:,1),'LineWidth',1.5)
hold on
plot(t,ReCon_1,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Resp.')
ylim([-0.1, 0.1])
legend('Original','Reconstructed')

for I = 1:totalpair
    subplot(7,1,I+1)
    plot(t,x_recon(:,1,I),'LineWidth',1.5)
    ylabel('Resp.')
    ylim([-0.1, 0.1])
    legend(['Group ',sprintf('%d',I)])
end

sgtitle('Channel 1: First Floor Response')

%Channel 2 (Floor 2) components
figure 
subplot(7,1,1)
plot(t,accel(:,2),'LineWidth',1.5)
hold on
plot(t,ReCon_2,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Resp.')
ylim([-0.1, 0.1])
legend('Original','Reconstructed')

for I = 1:totalpair
    subplot(7,1,I+1)
    plot(t,x_recon(:,2,I),'LineWidth',1.5)
    xlabel('Time (s)')
    ylabel('Resp.')
    ylim([-0.1, 0.1])
    legend(['Group ',sprintf('%d',I)])
end

sgtitle('Channel 2: Second Floor Response')

%Channel 3 (Floor 3) components
figure 
subplot(7,1,1)
plot(t,accel(:,3),'LineWidth',1.5)
hold on
plot(t,ReCon_3,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Resp.')
ylim([-0.1, 0.1])
legend('Original','Reconstructed')

for I = 1:totalpair
    subplot(7,1,I+1)
    plot(t,x_recon(:,3,I),'LineWidth',1.5)
    xlabel('Time (s)')
    ylabel('Resp.')
    ylim([-0.1, 0.1])
    legend(['Group ',sprintf('%d',I)])
end

sgtitle('Channel 3: Third Floor Response')

%Channel 4 (Floor 4) components
figure 
subplot(7,1,1)
plot(t,accel(:,4),'LineWidth',1.5)
hold on
plot(t,ReCon_4,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Resp.')
ylim([-0.1, 0.1])
legend('Original','Reconstructed')

for I = 1:totalpair
    subplot(7,1,I+1)
    plot(t,x_recon(:,4,I),'LineWidth',1.5)
    xlabel('Time (s)')
    ylabel('Resp.')
    ylim([-0.1, 0.1])
    legend(['Group ',sprintf('%d',I)])
end

sgtitle('Channel 4: Fourth Floor Response')

%Channel 5 (Floor 5) components
figure 
subplot(7,1,1)
plot(t,accel(:,5),'LineWidth',1.5)
hold on
plot(t,ReCon_5,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('SV0')
ylim([-0.1, 0.1])
legend('Original','Reconstructed')

for I = 1:totalpair
    subplot(7,1,I+1)
    plot(t,x_recon(:,5,I),'LineWidth',1.5)
    xlabel('Time (s)')
    ylabel('Resp.')
    ylim([-0.1, 0.1])
    legend(['Group ',sprintf('%d',I)])
end

sgtitle('Channel 5: Fifth Floor Response')

%Channel 6 (Floor 6) components
figure 
subplot(7,1,1)
plot(t,accel(:,6),'LineWidth',1.5)
hold on
plot(t,ReCon_6,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Resp.')
ylim([-0.1, 0.1])
legend('Original','Reconstructed')

for I = 1:totalpair
    subplot(7,1,I+1)
    plot(t,x_recon(:,6,I),'LineWidth',1.5)
    xlabel('Time (s)')
    ylabel('Resp.')
    ylim([-.1,.1])
    legend(['Group ',sprintf('%d',I)])
end

sgtitle('Channel 6: Sixth Floor Response')

%Plot comparison of original signals and reconstructed signals for better
%visualization
figure
subplot(3,1,1)
plot(t,accel(:,1),'LineWidth',1.5)
hold on
plot(t,ReCon_1,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp.')
legend('Original','Reconstructed')
title('Channel 1')

subplot(3,1,2)
plot(t,accel(:,2),'LineWidth',1.5)
hold on
plot(t,ReCon_2,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp.')
legend('Original','Reconstructed')
title('Channel 2')

subplot(3,1,3)
plot(t,accel(:,3),'LineWidth',1.5)
hold on
plot(t,ReCon_3,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp.')
legend('Original','Reconstructed')
title('Channel 3')

sgtitle('Channels 1-3: Original vs. Reconstructed')


figure
subplot(3,1,1)
plot(t,accel(:,4),'LineWidth',1.5)
hold on
plot(t,ReCon_4,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp.')
legend('Original','Reconstructed')
title('Channel 4')

subplot(3,1,2)
plot(t,accel(:,5),'LineWidth',1.5)
hold on
plot(t,ReCon_5,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp.')
legend('Original','Reconstructed')
title('Channel 5')

subplot(3,1,3)
plot(t,accel(:,6),'LineWidth',1.5)
hold on
plot(t,ReCon_6,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Amp.')
legend('Original','Reconstructed')
title('Channel 6')

sgtitle('Channels 4-6: Original vs. Reconstructed')


%Plot FRF comparison
f = t - t(1);
f = fs*f/f(end);

%Channel 1
figure
subplot(ch+1,1,1)
plot(f,2*abs(fft(accel(:,1))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_1)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Mag.')
xlim([0 18])
legend('Original', 'Reconstructed')

for I = 1:totalpair
    subplot(ch+1,1,I+1)
    plot(f,2*abs(fft(accel(:,1))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,1,I))),'LineWidth',1.5)
    xlabel('Frequency (Hz)')
    ylabel('Mag.')
    xlim([0 18])
    legend('Original', ['Group ',sprintf('%d',I)])
end

sgtitle('Channel 1: FFT of SV Components')

%Channel 2
figure
subplot(ch+1,1,1)
plot(f,2*abs(fft(accel(:,2))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_2)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Mag.')
xlim([0 18])
legend('Original', 'Reconstructed')

for I = 1:totalpair
    subplot(ch+1,1,I+1)
    plot(f,2*abs(fft(accel(:,2))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,2,I))),'LineWidth',1.5)
    ylabel('Mag.')
    xlim([0 18])
    legend('Original', ['Group ',sprintf('%d',I)])
end

sgtitle('Channel 2: FFT of SV Components')

%Channel 3
figure
subplot(ch+1,1,1)
plot(f,abs(fft(accel(:,3))),'--','LineWidth',1.5)
hold on
plot(f,abs(fft(ReCon_3)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Mag.')
xlim([0 18])
legend('Original', 'Reconstructed')

for I = 1:totalpair
    subplot(ch+1,1,I+1)
    plot(f,2*abs(fft(accel(:,3))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,3,I))),'LineWidth',1.5)
    xlabel('Frequency (Hz)')
    ylabel('Mag.')
    xlim([0 18])
    legend('Original', ['Group ',sprintf('%d',I)])
end

sgtitle('Channel 3: FFT of SV Components')

%Channel 4
figure
subplot(ch+1,1,1)
plot(f,2*abs(fft(accel(:,4))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_4)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Mag.')
xlim([0 18])
legend('Original', 'Reconstructed')

for I = 1:totalpair
    subplot(ch+1,1,I+1)
    plot(f,2*abs(fft(accel(:,4))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,4,I))),'LineWidth',1.5)
    xlabel('Frequency (Hz)')
    ylabel('Mag.')
    xlim([0 18])
    legend('Original', ['Group ',sprintf('%d',I)])
end

sgtitle('Channel 4: FFT of SV Components')

%Channel 5
figure
subplot(ch+1,1,1)
plot(f,2*abs(fft(accel(:,5))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_5)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Mag.')
xlim([0 18])
ylim([0,50])
legend('Original', 'Reconstructed')

for I = 1:totalpair
    subplot(ch+1,1,I+1)
    plot(f,2*abs(fft(accel(:,5))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,5,I))),'LineWidth',1.5)
    xlabel('Frequency (Hz)')
    ylabel('Mag.')
    xlim([0 18])
    ylim([0,50])
    legend('Original', ['Group ',sprintf('%d',I)])
end

sgtitle('Channel 5: FFT of SV Components')

%Channel 6
figure
subplot(ch+1,1,1)
plot(f,2*abs(fft(accel(:,6))),'--','LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_6)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel('Mag.')
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')

for I = 1:totalpair
    subplot(ch+1,1,I+1)
    plot(f,2*abs(fft(accel(:,6))),'--','LineWidth',1.5)
    hold on
    plot(f,2*abs(fft(x_recon(:,6,I))),'LineWidth',1.5)
    xlabel('Frequency (Hz)')
    ylabel('Mag.')
    xlim([0 18])
    ylim([0,50])
    legend('Original', ['Group ',sprintf('%d',I)])
end

sgtitle('Channel 6: FFT of SV Components')

%Original vs. Reconstructed FFT for each channel, larger for easier
%visualization

figure
subplot(3,1,1)
sgtitle('Ch 1-3: FFT of Original vs. Reconstructed Signal')
plot(f,2*abs(fft(accel(:,1))),'LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_1)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel(['Amp.'])
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')
title('Channel 1')

subplot(3,1,2)
plot(f,2*abs(fft(accel(:,2))),'LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_2)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel(['Amp.'])
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')
title('Channel 2')


subplot(3,1,3)
plot(f,2*abs(fft(accel(:,3))),'LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_3)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel(['Amp.'])
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')
title('Channel 3')


figure
subplot(3,1,1)
sgtitle('Ch 4-6: FFT of Original vs. Reconstructed Signal')
plot(f,2*abs(fft(accel(:,4))),'LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_4)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel(['Amp.'])
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')
title('Channel 4')


subplot(3,1,2)
plot(f,2*abs(fft(accel(:,5))),'LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_5)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel(['Amp.'])
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')
title('Channel 5')


subplot(3,1,3)
plot(f,2*abs(fft(accel(:,6))),'LineWidth',1.5)
hold on
plot(f,2*abs(fft(ReCon_6)),'LineWidth',1.5)
xlabel('Frequency (Hz)')
ylabel(['Amp.'])
xlim([0 18])
ylim([0 50])
legend('Original', 'Reconstructed')
title('Channel 6')