%% start and stop MDF
clc
clear
close all
%% load  legpress 30% dataset
num=1:24;
num([5,6,16])=[];
feature = [];
Fs = 1000;
wo = 50/(Fs/2);
bw = wo/30;
[b0,a0] = iirnotch(wo,bw);          %% 50 Hz notch filter

w1 = 100/(Fs/2);
bw1 = wo/30;
[b1,a1] = iirnotch(w1,bw1);          %% 100 Hz notch filter

n=50;
f=[0 0.03 0.04  1];
a2=[0 0 0.5 1];
b2=firls(n,f,a2);                   %% band pass filter 20 to 500 Hz
window=4*Fs;
%% feature extraction
% we extract the start and stop sample time for our dataset by hand
start = [2400,1001,6300,3500,2500,2900,11000,3500,3200,9500,4000,3500,2100,2600,5000,4200,16000,7800,4000,3500,2000];
stop = [138500 52000 56500 67000 100000 42000 114000 86000 89500 127000 77500 110000 65500 65000 85000 86500 72000 165000 154500 157000 170000];

for I=1:size(num,2)
    name=([num2str(num(I)),'.mat']);
    load(name)
    Fs=sampfreq;
    EMG=datablock1.data';
    ch_name = sources;
    %% fatigue feature
    count = 1;
    for J=1:channels
        normalize_EMG = (EMG(J,:)-mean(EMG(J,:)))/std(EMG(J,:));
        signal_notch_50 = filtfilt(b0,a0,normalize_EMG);
        signal_notch_100 = filtfilt(b1,a1,signal_notch_50);
        filterd_signal = filtfilt(b2,1,signal_notch_100);
        %% calculate MDF feature for the begining and the end of the excercise within 2 second window
        start_signal = filterd_signal(start(I)-1000:start(I)+1000);
        stop_signal = filterd_signal(stop(I)-1000:stop(I)+1000);
        L = length(start_signal);
        NFFT = 2^nextpow2(L);
        fft0 = fft(start_signal,NFFT);
        fft1 = fft(stop_signal,NFFT);
        freq = Fs/2*linspace(0,1,NFFT/2+1);
        abs_fft0 = abs(fft0/L);
        mdf(1,:) = abs_fft0(1:NFFT/2+1);
        mdf(1,2:end-1) = 2*mdf(1,2:end-1);
        abs_fft1 = abs(fft1/L);
        mdf(2,:) = abs_fft1(1:NFFT/2+1);
        mdf(2,2:end-1) = 2*mdf(2,2:end-1);
        for k=1:2
            c=0;
            for i=1:size(freq,2)
                if c< sum(mdf(k,:))/2
                    c=c+mdf(i);
                elseif c> sum(mdf(k,:))/2
                    MDF_ss(I,count)=freq(i);
                    count = count+1;
                    break
                end
            end
        end
    end
end
%% for 
result = round(MDF_ss,2);
colm_name = ["subject"; "quadro_rectus_femoris_beginning";"quadro_rectus_femoris_end";
    "quadro_vastus_m_beginning";"quadro_vastus_m_end";"quadro_vastus_l_beginning";
    "quadro_vastus_l_end";"bicep_femoris_beginning";"bicep_femoris_end"];

subject = num';
result_table = table(subject,result(:,1),result(:,2),result(:,3),result(:,4),result(:,5),result(:,6),result(:,7),result(:,8), 'VariableNames', colm_name);
disp('MEdian ferquency of quadroceps muscles in leg press')
disp(result_table)
writetable(result_table,"MDF_leg_press")