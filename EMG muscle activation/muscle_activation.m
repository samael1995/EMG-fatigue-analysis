clc
clear
close all
%% load data
EMG = load('sample.mat');
ch_name = EMG.sources;
data=EMG.datablock1.data;
Fs=EMG.sampfreq;
[L, chanel] = size(data);
t = (1:L)/Fs;

%% preprocessing
% zero mean and rectified signal
mu = mean(data(1:5*Fs,:));
MAV = abs(data-mu);

% signal smoothing with MA filter
window=500;        %% sample equivalant to window/Fs second
for i = 1:chanel
    EMG = MAV(:,i);
    for j = 1:L - window
        MA(j,i)= mean(EMG(j:j+window));
    end
end
%% identifying activation of EMG muscle based on three method of RMS, quartile and mean thresholding
% RMS thresholding
L2= length(MA);
for i=1:chanel
    window = 200;
    for j=1:floor(L2/window)-1
        mu(j) = sqrt(mean(MA(1+(j-1)*window: window*j,i).^2));
    end
    MU(i) = mean(mu);
end
threshold = MU*0.8;

for i=1:chanel
    for j= 1:L2
        if MA(j,i)>threshold(i)
            RMS_T(j,i)=1;
        else
            RMS_T(j,i)=0;
        end
    end
end
% figure(1)
subplot(2,1,1)
plot((1:length(MA(:,1)))/Fs, MA(:,1))
title("smoothed EMG signal")
xlabel("time")
ylabel("amplitude")
subplot(2,1,2)
plot((1:length(RMS_T(:,1)))/Fs, RMS_T(:,1))
ylim([-0.5,1.5])
xlabel("time")
ylabel("amplitude")
title("RMS thresholding")

        
% quartile threshold
Q_threshold = 0.50;
Quartile = quantile(MA,Q_threshold);

for i=1:chanel
    for j= 1:L2
        if MA(j,i)>Quartile(i)
            Q_T(j,i)=1;
        else
            Q_T(j,i)=0;
        end
    end
end

figure(2)
subplot(2,1,1)
plot((1:length(MA(:,1)))/Fs, MA(:,1))
title("smoothed EMG signal")
xlabel("time")
ylabel("amplitude")
subplot(2,1,2)
plot((1:length(Q_T(:,1)))/Fs, Q_T(:,1))
ylim([-0.5,1.5])
xlabel("time")
ylabel("amplitude")
title([num2str(Q_threshold)," quartile thresholding"])


% mean threshold
Mean = mean(MA);
for i=1:chanel
    for j= 1:L2
        if MA(j,i)>Mean(i)
            m_T(j,i)=1;
        else
            m_T(j,i)=0;
        end
    end
end

figure(3)
subplot(2,1,1)
plot((1:length(MA(:,1)))/Fs, MA(:,1))
title("smoothed EMG signal")
xlabel("time")
ylabel("amplitude")
subplot(2,1,2)
plot((1:length(m_T(:,1)))/Fs, m_T(:,1))
ylim([-0.5,1.5])
xlabel("time")
ylabel("amplitude")
title("mean thresholding")

%% feature extraction:
% The first muscle to be activated at the beginning and the last muscle to be inactive
% duration of active time of muscle
% mean of EMG activity

offset = [];
onset = [];
for i= 1:chanel
    c=1;
    c1=1;
    for j= 2:L2-1
        if Q_T(j-1,i)==0 && Q_T(j+1,i)==1
            onset(i,c)=j;
            c=c+1;
        end
         if Q_T(j-1,i)==1 && Q_T(j+1,i)==0
             offset(i,c1)=j;
             c1=c1+1;
         end
    end
end
end_sample = max(offset');
start_sample = onset(:,1)';
end_time = end_sample/Fs;
start_time = start_sample/Fs;
active_time = (end_time-start_time);
for i=1:chanel
    Mean_amplitude_EMG(i) = mean(MA(start_sample(i):end_sample(i),i));
end
% 
%% result table
electrode = [ch_name(1,:);ch_name(2,:);ch_name(3,:);ch_name(4,:);ch_name(5,:);ch_name(6,:);ch_name(7,:);ch_name(8,:)];
colm_name = ["elctrode","start time","end time","active time", "mean active EMG"];
result_table = table(electrode,start_time',end_time',active_time',Mean_amplitude_EMG', 'VariableNames', colm_name);
disp(result_table)

