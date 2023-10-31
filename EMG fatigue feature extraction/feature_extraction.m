%% frequcy features
clc
clear
close all
%% load data
load('sample.mat')
data=datablock1.data;
Fs=sampfreq;
ch_name = ["quadro rectus femoris", "quadro vastus femoris M",  "quadro vastus femoris L", "biceps femoirs"];
figure(1)
t=(1:1:size(data,1))/Fs;
plot(t,data(:,1))
title('EMG signal')
xlabel("time(s)")
ylabel("amplitude")
legend(ch_name(1))
%% signal preoprocessing
%50 Hz notch filter 
WL=48*2/Fs;
WU=52*2/Fs;
Wn=[WL,WU];           % frequency band
n=4;                  % butterworth order
type='stop';
[b,a]= butter(n,Wn,type);
% 100 Hz  notch filter
WL=98*2/Fs;
WU=102*2/Fs;
Wn=[WL,WU];
n=4;                  % butterworth order
[b2,a2]= butter(n,Wn,type);

% 200 Hz notch filter
WL=198*2/Fs;
WU=202*2/Fs;
Wn=[WL,WU];
n=4;                   % butterworth order
[b3,a3]= butter(n,Wn,type);
%% frequncy feature
figure(2)
img=imread('FreqFeat.JPG');
imshow(img)
wt = 2;    % 2 second window
window = wt*Fs;
for i=1:channels
    filtered_signal_1=filter(b,a,data(:,i));
    filtered_signal_2=filter(b2,a2,filtered_signal_1);
    filtered_signal=filter(b3,a3,filtered_signal_2);
    for j=1:length(filtered_signal)/window
        S=filtered_signal((j-1)*window+1:j*window);
        L=length(S);
        NFFT=2^nextpow2(L);
        FFT=fft(S,NFFT);
        freq=Fs/2*linspace(0,1,NFFT/2+1);
        abs_FFT=abs(FFT(1:NFFT/2+1));
%         mean frequency
        MNF(i,j)=(freq*abs_FFT)/sum(abs_FFT);
%         median frequency
        c=0;
        for k=1:size(freq,2)
            if c< sum(abs_FFT)/2
                c=c+abs_FFT(k);
            elseif c> sum(abs_FFT)/2
                MDF(i,j)=freq(k);
                break
            end
        end
%        peak frequency
        [amp loc]=max(abs_FFT);
        PF(i,j)=freq(loc);
%         total power
        power(i,j)=sum(abs_FFT);
    end
    
end
%% ploting features
t = 1:wt:size(MNF,2)*wt;

figure(3)
for i=1:channels
    plot(t, MNF(i,:))
    hold on
end
xlabel("time(s)")
ylabel("frequency(Hz)")
legend(ch_name)
title ('mean freuency feature')

figure(4)
for i=1:channels
    plot(t, MDF(i,:))
    hold on
end
xlabel("time(s)")
ylabel("frequency(Hz")
legend(ch_name)
title ('Median frequency feature')


figure(5)

for i=1:4
    subplot(4,1,i)
    bar(t,PF(i,:))
    ylabel("frequency(Hz)")
    legend(ch_name(i))
end

xlabel("time(s)")
subplot(4,1,1)
title ('peak frequency feature')

figure(6)
for i=1:channels
    plot(t, power(i,:))
    hold on
end
xlabel("time(s)")
ylabel("magnitude")
legend(ch_name)
title ('power of EMG signal')
%% save features
feature=[];
feature(1:4,:)= MNF;
feature(5:8,:)= MDF;
feature(9:12,:)= PF;
feature(13:16,:)= power;

writematrix(feature,'extracted_features.xlsx','Sheet',1)