function [CP I] = CPSD2(x, y, seg_size)

% CPSD2 creates the cross-periodogram for the two input datas
full = floor(length(x)/seg_size);
x = x(1:seg_size*full);
y = y(1:seg_size*full);

% Reshapes the data into segments
x = reshape(x, seg_size, full);
y = reshape(y, seg_size, full);

% Fourier Transform
dx = fft(x);
dy = fft(y);

% Cross-periodogram
I = (1/(2*pi*seg_size)*dy.*conj(dx));

% Average across the segments
CP = mean(I, 2);