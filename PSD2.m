function [P, L, I] = PSD2(fxx,seg_size)
% PSD2 creates an auto-periodogram for the input data    
fxx = fxx - mean(fxx);    
L =  fix(length(fxx)/seg_size);

% Reshape the data into segments
fxx = fxx(1:L*seg_size);
fxx1 = reshape(fxx, seg_size, L);

% Create the periodograms
dx = fft(fxx1);
I = 1/(2*pi*seg_size)*abs(dx.*dx);

% average over the segments
P = mean(I,2);
P = P(2: length(P)/2+1);