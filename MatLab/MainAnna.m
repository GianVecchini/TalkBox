%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                             %
%             THE TALK BOX EFFECT             % 
%                                             %
%  Anna Chiara Melioli & Gianluigi Vecchini   %
%                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% IMPORT / CREATE SIGNALS
nCase = 1; % 1, 2, 3

switch nCase
    case 1
        % Observe that in this case Fs and Length is equal for the two
        % signals
        [s, Fs] = audioread('Toms_diner.wav');
        [x, ~] = audioread('moore_guitar.mp3');

        % Guitar is a stereo signal, we convert it to mono
        x = mean(x, 2); 
        sigLength = length(s);
    case 2
        [s, Fs] = audioread('Toms_diner.wav');
        [x,~] = audioread('tromba.mp3');

        sigLength = length(s);
    case 3
        % ?

end

% Such signals are non stationary, let's split them in short stationary
% frames. We will work separatly on those frame and then reconstruct the
% total signal using the OLA method.

%% DEFINE WINDOWING PARAMETERS 

% long frames --> high freq resolution and low time resolution
% short frames --> low freq resolution and high time resolution
wTimeLenght = 0.05; %con 50 invece di 80 sembra venire un po' meglio
%con 0.02 forse pure e l'envelope sembra meglio ma clippa
wLength = round(Fs*wTimeLenght)+1;

% A smaller hop gives an higher overlap. This results in having a smoother
% result but a havier computation.
hopSize = floor(wLength/2); 

% Type of window
w = hann(wLength,'periodic');  
% 'periodic' added only bc in the matlab documentation someone added it to check
% COLA condition through iscola() function

% Number of frames
nFrames = floor((sigLength-wLength)/hopSize) + 1;

% Overlap
ol =  (wLength - hopSize);

% Check of cola condition
if iscola(w,ol,'ola')
    disp('COLA condition satisfied');
else
    disp('COLA condition is not satisfied!!!')
end

%% TALK BOX FILTER

% Output signal
outLength = (nFrames - 1) * hopSize + wLength;
y = zeros(1, outLength);
e = zeros(1, outLength);

% Filter order: is the numer of poles in the filter and better is the
% approximation of the spectral envelope. Some observations:
%   - High P_s means that we have a better approximation of the spectra. 
%     If such value is too high then the result will appear so similar to
%     the original s signal and !almost not filtered
%   - Low P_s the words are less distinguishible
%   - High P_x means that the residual does not have any spectral content,
%     is less correlated and more similar to a white noise
%   - Low P_x gives a more rich excitation
P_s = 46; % Fs/1000 < p < Fs/1000+4 --> 44.1 < p < 48.1
P_x = 46;

% Compute for each frame
for i = 1:nFrames
    % Index of the frame
    idxFrame = (i-1) * hopSize + 1 : (i-1) * hopSize + wLength;
    
    % Extract frames using window 
    sFrame = s(idxFrame) .* w;
    xFrame = x(idxFrame) .* w;
    
    %--- FIND SHAPING FILTER: work with s ---%
    % Observe that we implement LPC using the the Levison-Durbin recursion. 
    % Such algorithm has a lighter computational cost

    % Auto-correlation computation of the frame
    [r_s, rlags_s] = xcorr(sFrame);

    % Extract the required autocorrelation coefficients
    central_lag_s = ceil(length(r_s) / 2);

    % Apply Levinson-Durbin recursion
    a_s = levinson(r_s(central_lag_s:(central_lag_s+P_s)), P_s);

    % Prediction error (not needed)
    e_s = conv(sFrame, a_s, 'same'); 

     %--- FIND EXCITATION: work with x ---%
    % Auto-correlation computation
    [r_x, rlags_x] = xcorr(xFrame);

    % Extract the required autocorrelation coefficients
    central_lag_x = ceil(length(r_x) / 2);

    % Apply Levinson-Durbin recursion
    a_x = levinson(r_x(central_lag_x:(central_lag_x+P_x)), P_x);
    
    % Excitation (prediction error)
    e_x = conv(xFrame, a_x, 'same'); 

    %--- RECONSTRUCTION: combine s' filter with x' excitation ---%
    yFrame = filter(1, a_s, e_x);

    % OLA
    y(idxFrame) = y(idxFrame) + (yFrame)';
    e(idxFrame) = e(idxFrame) + (e_x)';
    
    % For one frame let's plot spectrograms
    if i == 50
        Nfft = 1024; %2^(ceil(log2( L + M -1))); ok non centra forse

        S = abs(fft(sFrame, Nfft));
        S = S(1:Nfft/2+1);

        E_X = abs(fft(e_x, Nfft));
        E_X = E_X(1:Nfft/2+1);

        Y = abs(fft(yFrame, Nfft));
        Y = Y(1:Nfft/2+1);

        [H_s, w_s] = freqz(1, a_s, Nfft/2+1, Fs);
        
        S = 20*log10(S);
        E_X = 20*log10(E_X);
        H_s = 20*log10(abs(H_s));
        Y = 20*log10(Y);

        figure;
        plot(w_s, H_s - max(H_s), 'r', 'LineWidth', 2); hold on;
        plot(linspace(0, Fs/2, Nfft/2+1), S - max(S), 'k','LineWidth', 1); hold on;
        plot(linspace(0, Fs/2, Nfft/2+1), E_X - max(E_X), 'b', 'LineWidth', 1); hold on;
        plot(w_s, Y - max(Y), 'm', 'LineWidth', 1); 
        legend('LPC envelope','Spectrum of Speech', 'Excitation', 'Result');
        title('Speech Frame - LPC Envelope')
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')

        % Observe that we are plotting in dB, subtracting the maximum
        % means setting the peak to 0. 
    end
end

y = y ./ max(abs(y) + eps); % Normalize the output
sound(y, Fs)
audiowrite('effectOutput.wav', y, Fs);

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


%%
