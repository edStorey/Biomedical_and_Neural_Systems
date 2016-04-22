clear all
% Set segment size, sampling frequency and load ramp and hold data
seg_size = 1024;
fs = 1000;
load('jm201a.mat');

% Reshape the data so that the trials are segmented and consecutive, the
% ramp stages are not included
for resize = 1:length(trig_start_hold)
    if resize == 85
        trial_start = trig_start_hold(resize); 
        trial_stop = trig_stop_hold(resize); 
        trial_length = samp_hold(resize);
    else
    trial_start = trig_start_hold(resize); 
    trial_stop = trig_stop_hold(resize); 
    trial_length = samp_hold(resize);
    end
    seg_num = floor(trial_length/seg_size);
    for seg_sep = 1:seg_num
        
        if resize == 1
            s1 = eeg(trial_start+(seg_sep-1)*seg_size+1:trial_start+seg_sep*seg_size);
            s2 = emg(trial_start+(seg_sep-1)*seg_size+1:trial_start+seg_sep*seg_size);
            s3 = force(trial_start+(seg_sep-1)*seg_size+1:trial_start+seg_sep*seg_size)/max(force(trial_start:trial_stop));
        else
            s1(end+1:end+seg_size) = eeg(trial_start+(seg_sep-1)*seg_size+1:trial_start+seg_sep*seg_size);
            s2(end+1:end+seg_size) = emg(trial_start+(seg_sep-1)*seg_size+1:trial_start+seg_sep*seg_size);
            s3(end+1:end+seg_size) = force(trial_start+(seg_sep-1)*seg_size+1:trial_start+seg_sep*seg_size)/max(force(trial_start:trial_stop));
        end
    end
    
end

% Full-wave rectify the emg signal
s2 = abs(s2- mean(s2));

% Create Periodograms and cross-periodograms
[Px,seg_tot,Ixx] = PSD2(s1,seg_size);
[Px,seg_tot,Iyy] = PSD2(s2,seg_size);
[CP Ixy] = CPSD2(s1, s2, seg_size);
clear Px CP

% Calculate R and Q
R = psi(1,1);
Q = 0.1;

% Reshape the data and take the log
z = log([Ixx(1:seg_size/2,:);Iyy(1:seg_size/2,:);real(Ixy(1:seg_size/2,:));imag(Ixy(1:seg_size/2,:))]);

% If there are any 0s then they will be converted to -Inf here the lower
% limit has been set to -100
size_z = size(z);
for a = 1:size_z(2)
    for b = 1:size_z(1)
        if z(b,a) == -Inf || z(b,a) == NaN
            z(b,a) = -100;
        end
    end
end

% Initialise variables z, x and P
x = zeros(size(z));
x(:,1) = z(:,1);
P(1) = R;

% initialise x2 and Pp2
xp(:,2) = x(:,1);
Pp(2) = P(1)+Q; 

% Run the filter equations for each segment 
for seg = 2: seg_tot
K(seg) =  Pp(seg)/(Pp(seg) + R);
x(:,seg) = xp(:,seg) + K(seg)*(z(:,seg) - xp(:,seg));

P(seg) = (1 - K(seg))*Pp(seg);
xp(:,seg+1) = x(:,seg);
Pp(seg+1) = P(seg) + Q;
end




% Initialise the first collumns and values of the estimates for x and P
x_hat = zeros(size(z));
x_hat(:,seg_tot) = x(:,seg_tot);
P_hat(seg_tot) = P(seg_tot);

% run the filter in the reverse direction
for seg = seg_tot-1:-1:1
    A(seg) = P(seg)/(Pp(seg+1));
    x_hat(:,seg) = x(:,seg) + A(seg)*(x_hat(:,seg+1)-xp(:,seg+1));
    P_hat(seg) = P(seg) + A(seg)*(P_hat(seg+1)-Pp(seg+1))*A(seg)'; 
    
end


% Here we convert the data out of the log form by taking the exponential
% for the forward filtered data
cxx1 = exp(x(1:seg_size/2,:));
cyy1 = exp(x(seg_size/2+1:seg_size,:));
cyx1 = complex(real(exp(x(seg_size+1:seg_size+seg_size/2,:))) ,  real((exp(x(seg_size+seg_size/2+1:seg_size*2,:)))));
c_hat_size = size(cxx1);

% Calculate the coherence for the forward filtered data
c_hat1 = (abs(cyx1).^2)./(cyy1.*cxx1);
c_hat1 = c_hat1(2:c_hat_size(1),:);
c_hat1_10 = (c_hat1(10,:));
c_hat1_mean = mean(c_hat1,1);

% Here we convert the data out of the log form by taking the exponential
% for the reverse filtered data
cxx = exp(x_hat(1:seg_size/2,:));
cyy = exp(x_hat(seg_size/2+1:seg_size,:));
cyx = complex(real(exp(x_hat(seg_size+1:seg_size+seg_size/2,:))) ,  real((exp(x_hat(seg_size+seg_size/2+1:seg_size*2,:)))));

% Calculate the coherence for the reverse filtered data
c_hat = (abs(cyx).^2)./(cyy.*cxx);
c_hat = c_hat(2:c_hat_size(1),:);
c_hat1_mean = mean(c_hat1,1);
c_hat_10 = c_hat(10,:);




% Initialse axis variables
samp_period = 1/fs;
T = linspace(0, length(s1)/fs, length(s1)); % Data for the length
Ts = transpose(T);
deltaf=1/(seg_size*samp_period);
freq_max = 100;
freq=(1:freq_max)'*deltaf;


% From here images of the various data collected are plotted with a white
% background and labelled axis. All plot the data after it has been reverse filtered

% The EEG auto-spectra
figure(1);
[h] = imagesc(Ts,flipud(freq),log10(cxx(:,2:round(freq_max/(seg_size/fs))))');

h_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(h_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('Filtered EEG Spectrum', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)

% The EMG auto-spectra
figure(2);
[h] = imagesc(Ts,flipud(freq),log10(cyy(:,2:round(freq_max/(seg_size/fs))))');

h_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(h_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('Filtered EMG Spectrum', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)

% The EMG-EEG cross-spectra
figure(3);
[h] = imagesc(Ts,flipud(freq),log10(abs(cyx(:,2:round(freq_max/(seg_size/fs)))))');

h_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(h_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('Filtered cross EEG EMG Spectrum' ,'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)

% The coherence across the EEG and EMG signals
figure(4);
[h] = imagesc(Ts,flipud(freq),c_hat(:,2:round(freq_max/(seg_size/fs)))');

h_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(h_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('The Coherence', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)



