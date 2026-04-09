function spectrograms(s,x,y,e,Fs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Plot the spectrograms, for curiosity
figure;
subplot(4,1,1)
spectrogram(s, hann(512), 256, 1024, Fs, 'yaxis')
title('Speech Spectrogram')

subplot(4,1,2)
spectrogram(x, hann(512), 256, 1024, Fs, 'yaxis')
title('Music Spectrogram')

subplot(4,1,3)
spectrogram(y, hann(512), 256, 1024, Fs, 'yaxis')
title('Output Spectrogram')

subplot(4,1,4)
spectrogram(e, hann(512), 256, 1024, Fs, 'yaxis')
title('Excitation Spectrogram')
end