clc,clear,close all
opengl hardware; 
% >>>>> Data <<<<<
Data = importdata('Homework-2 data set-RCF-Four Specimen Test Data.xlsx');
Data = Data.data; % ！！！！！！！！
Time = Data(4001:1:18001,1);
t = Data(:,1);

SR = round(1/(t(2)-t(1)));
% Case1
Data6_04 = Data(4001:1:18001,3) ;  % 單位: g /手動改Case
Data6_10 = Data(4001:1:18001,2) ;  % 單位: g 
Acc_rel1 = Data6_04-Data6_10;
% Case2
%Data4_04 = Data(:,7) ;  % 單位: g 
%Data4_10 = Data(:,6) ;  % 單位: g 
%Acc_rel2 = Data4_04-Data4_10;

% >>>>> WPT <<<<<
wpt = wpdec(Acc_rel1,3,'db1');
figure
plot(wpt) % Wavelet Packet Tree
% saveas(gcf,fullfile(cd,'Wavelet Packet Tree.png'))
% saveas(gcf,fullfile(cd,'Wavelet Packet Tree.fig'))

cfs32 = wpcoef(wpt,[3,2]); % Wavelet packet coefficients
rcfs32 = wprcoef(wpt,[3,2]); % 重建wavelet packet

cfs33 = wpcoef(wpt,[3,3]);
rcfs33 = wprcoef(wpt,[3,3]);

cfs34 = wpcoef(wpt,[3,4]);
rcfs34 = wprcoef(wpt,[3,4]);


figure
subplot(3,1,1)
plot(Acc_rel1) ; title('Original signal')
axis tight

subplot(3,1,2)
plot(cfs32) ; title('Wavelet Packet (3,2) coefficients')
axis tight

subplot(3,1,3)
plot(rcfs32) ; title('Reconstruct Wavelet Packet (3,2) coefficients')
axis tight
% saveas(gcf,fullfile(cd,'Packet (3,2).png'))
% saveas(gcf,fullfile(cd,'Packet (3,2).fig'))

figure
subplot(3,1,1)
plot(Acc_rel1) ; title('Original signal')
axis tight

subplot(3,1,2)
plot(cfs33) ; title('Wavelet Packet (3,3) coefficients')
axis tight

subplot(3,1,3)
plot(rcfs33) ; title('Reconstruct Wavelet Packet (3,3) coefficients')
axis tight
% saveas(gcf,fullfile(cd,'Packet (3,3).png'))
% saveas(gcf,fullfile(cd,'Packet (3,3).fig'))

figure
subplot(3,1,1)
plot(Acc_rel1) ; title('Original signal')
axis tight

subplot(3,1,2)
plot(cfs34) ; title('Wavelet Packet (3,4) coefficients')
axis tight

subplot(3,1,3)
plot(rcfs34) ; title('Reconstruct Wavelet Packet (3,4) coefficients')
axis tight
% saveas(gcf,fullfile(cd,'Packet (3,4).png'))
% saveas(gcf,fullfile(cd,'Packet (3,4).fig'))


%%
[t,F,E,A,T,C] = wptspect(Acc_rel1,200,11); %要分號嗎
% SPectrogram using WPT
function [t,F,E,A,T,C] = wptspect(u,SR,Level)

if size(u,1)==1
    u=u';
end
wname = 'bior6.8';

dt=1/SR;
N=max(size(u));
t=[0:dt:(N-1)*dt];

% Index of Player order
index=0;
for I=1:1:Level
    index=[index fliplr(index)+2^(I-1)];
end
index=index+1;

% Wavelet packet transform
T=wpdec(u,Level,wname);

A=zeros(length(u),2^Level);
for I = 1:1:2^Level
    %Reconstruction of all components
    C(:,I)=wprcoef(T,[Level,index(I)-1]);

    %Hilbert amplitude of all components
    A(:,I)=abs(hilbert(wprcoef(T,[Level,index(I)-1])));
end
% Energy = square of Hilbert amplitude
E=A.*A;

%Frequency resolution
F=linspace(0,1/dt/2,2^Level)';

%Interest freq
temp=find(F<=50);

% plot T-F spectrogram
figure('Position',[200 50 500 300])
surf(t,F(temp),E(:,temp)','LineStyle','none','FaceColor','interp','FaceLighting','phong',...
    'EdgeColor','interp'),view(0,90),colormap('jet'),ylim([0 50]),xlabel('Time(sec)'),...
    ylabel('Frequency(Hz)'),xlim([0 t(end)])

%adjust color display of spectogram(optional)
%caxis([0 0.2])

end





