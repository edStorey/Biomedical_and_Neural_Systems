clear all
% If L and int are both 2^a and 2^b then they will be divisible and the
% algorithm will run quicker without errors
seg_size = 2^10;
fs = 1000;
NW = 2.5;

% Number of overlapping segments
Int = 2^7;

% Load ramp and hold data
load('jm201a.mat')

% Create slepian windows
[E,V] = dpss(seg_size,NW);
size_E = size(E);
t = 1/fs:1/fs:length(emg)*20/fs;

% The commented parts here are test data

% f1 = 50;
% f2 = 100;
% wave = sin((2*pi*t)*f1);
% wave2 = sin((2*pi*t)*f2);
%  x = wave';
%  y = wave2';
%  x = randn(length(emg),1);
%  y = x + randn(length(emg),1);
%  y = randn(length(emg),1)
%  y = x;

% set the inputs to be the eeg and emg data
x = eeg;
y = emg;

% full wave rectification
if x == emg    
    x = abs(x - mean(x));
elseif y == emg
    y = abs(y - mean(y));
end


ramp_dat = struct([]);
segment_totals = zeros(85,1);
ramps = 0;
h = waitbar(ramps/length(trig_start_hold),'Creating Tapered Data and Averaging');
for ramps = 1:length(trig_start_hold)
    
    % calculate segments
    segments = floor((samp_hold(ramps))/seg_size);
    
    % Amount of extra segments
    Int_segments = segments + (segments-1)*(Int-1);
    ramp_segs = 0;
    
    % here the data is taken out of the original set and ordered into
    % segments and trials
    for ramp_segs = 0:(Int_segments-1)
        dat_start = trig_start_hold(ramps)+ramp_segs*(seg_size/Int);
        dat_end = dat_start + seg_size -1;
        seg_temp = x(dat_start:dat_end);
        seg_temp2 = y(dat_start:dat_end);
        
        for a = 1:size_E(2)
            % tapering is done here
            tapered(:,a) = fft(seg_temp.*E(:,a))/seg_size;
            tapered2(:,a) = fft(seg_temp2.*E(:,a))/seg_size;
            
        end
        
        % Three different sets of taperes are shown, the last is the
        % cross-data so isn't squared or absed
        mtapered = (abs((tapered)).^2);
        mtapered2 = (abs((tapered2)).^2);
        mtapered12 = ((tapered.*conj(tapered2)));
        
        % smooth across tapered data
        ramp_data{ramps,ramp_segs+1} = mean(mtapered,2);
        ramp_data2{ramps,ramp_segs+1} = mean(mtapered2,2);
        ramp_data12{ramps,ramp_segs+1} = mean(mtapered12,2);
    end

   % record the number of segments here
   segment_totals(ramps,1) = Int_segments;
   segment_totals(ramps,2) = segments;
   
   waitbar(ramps/length(trig_start_hold), h,'Creating Tapered Data and Averaging');
end
close(h)

size_out = size(ramp_data);

% Currently this takes an average of all the segments in a trial, is there
% any point to this, does it make segmentsing pointless?
row_select = 0;
select = 1;
g = waitbar(row_select/size_out(2),'Creating Averaged Periodograms');
for row_select = 1:size_out(2)
    clear Row_dat Row_dat2 Row_dat12 Row_dat_coh
    clear IRx IRy IRxy
    
    % when averaging across segments size_out(1) becomes segment_totals(row_select)
    for move = 1:size_out(1)
        
        if segment_totals(move) >= row_select
            
            % take the tapered data out of the various structs
            Row_dat(select,:) = ramp_data{move,row_select};
            Row_dat2(select,:) = ramp_data2{move,row_select};
            Row_dat12(select,:) = ramp_data12{move,row_select};

            
            select = select+1;
        end
    end
    
    % Take the mean across all the trials for the current segment
    IRx = mean(Row_dat,1);
    IRy = mean(Row_dat2,1);
    IRxy = mean(Row_dat12,1);
    
    % Store the averaged values in a new matrix
    Row_Autocov_x(:,row_select) = IRx;
    
    Row_Autocov_y(:,row_select) = IRy;
    
    Row_Crosscov(:,row_select) = IRxy;
    select = 1;
    
    waitbar(row_select/size_out(2), g,'Creating Averaged Periodograms');
end

close(g)

% Coherence
chyx = (abs(Row_Crosscov).^2)./((Row_Autocov_x).*(Row_Autocov_y));
size_c = size(chyx);
freq = 0:fs/seg_size:100;
time = 1:1/Int:9;

% From here images of the various data collected are plotted with a white
% background and labelled axis.

% Plotting the coherence
figure(1)
[h] = imagesc(time,flipud(freq),flipud(chyx(1:round((seg_size/fs)*100),:)));
hold on
h_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(h_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('The Multi-Taper Coherence', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)
hold off

% The tapered EMG data
figure(2)
[g] = imagesc(time,flipud(freq),flipud(log10(Row_Autocov_y(1:round((fs\seg_size)*100),:))));
hold on
g_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(g_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('The Multi-Tapered EMG Data', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)
hold off

% The tapered EEG data
figure(3)
[p] = imagesc(time,flipud(freq),flipud(log10(Row_Autocov_x(1:round((fs\seg_size)*100),:))));
hold on
p_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(p_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
axis([1 9 1 100])
set(gcf,'color','w');
title('The Multi-Tapered EEG Data', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)
hold off

% The tapered EEG-EMG cross-spectrum
figure(4)
[q] = imagesc(time,flipud(freq),flipud(log10(abs(Row_Crosscov(1:round((fs\seg_size)*100),:)))));
hold on
q_plot=gca;
YTick_lab=get(gca,'YTickLabel');
YLim=get(gca,'YLim');
YTick_pos=get(gca,'YTick')';
set(q_plot,'YTick',YLim(2)-flipud(YTick_pos)+YLim(1),'YTickLabel',flipud(YTick_lab))
set(gcf,'color','w');
title('The Multi-Tapered Cross-Spectrum', 'FontSize', 20)
ylabel('Frequency', 'FontSize', 20)
xlabel('Time (1024 sample segments)', 'FontSize', 20)
hold off


