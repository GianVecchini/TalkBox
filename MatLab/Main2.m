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
        [s, Fs] = audioread('giovannigiorgio.mp3');
        [x,~] = audioread('daftpunk.mp3');

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

wLength = round(Fs*wTimeLenght);
if rem(wLength,2) ~= 0
    wLength = wLength +1;
end

% A smaller hop gives an higher overlap. This results in having a smoother
% result but a havier computation.
hopSize = floor(wLength/2); 

% Type of window
w = hann(wLength,'periodic');  
% 'periodic' added only bc in the matlab documentation someone added it to check
% COLA condition through iscola() function

% Overlap
ol =  (wLength - hopSize);

% Check of cola condition
if iscola(w,ol,'ola')
    disp('COLA condition satisfied');
else
    disp('COLA condition is not satisfied!!!')
end

%% TALK BOX FILTER
% Filter order: is the numer of poles in the filter and better is the
% approximation of the spectral envelope. Some observations:
%   - High P_s means that we have a better approximation of the spectra. 
%     If such value is too high then the result will appear so similar to
%     the original s signal and !almost not filtered
%   - Low P_s the words are less distinguishible
%   - High P_x means that the residual does not have any spectral content,
%     is less correlated and more similar to a white noise
%   - Low P_x gives a more rich excitation 

% the rule for speech should be 
% Fs/1000 < p < Fs/1000+4 --> 44.1 < p < 48.1
P_s = 46; 
P_x = 70;
Nfft = 1024;

% doing it as a function so that we can try different parameters and all
% the parameters are out 
[y, e] = LPCandOLA(s, x, w, Fs, wLength, hopSize, P_s, P_x, Nfft);

y = y ./ max(abs(y) + eps); % Normalize the output
sound(y, Fs)
audiowrite('effectOutput.wav', y, Fs);

spectrograms(s,x,y,e,Fs);


% Experiments for the case 1! it changes a lot basing on soectrum 
%% other P_x example with some frequencies that are not shaped 
P_x = 50;
[y, e] = LPCandOLA(s, x, w, Fs, wLength, hopSize, P_s, P_x, Nfft);

y = y ./ max(abs(y) + eps); % Normalize the output
sound(y, Fs)
audiowrite('dinerLowp.wav', y, Fs);

spectrograms(s,x,y,e,Fs);

%% other P_x example too white --> we cannot hear the guitar anymore
P_x = 100;
[y, e] = LPCandOLA(s, x, w, Fs, wLength, hopSize, P_s, P_x, Nfft);

y = y ./ max(abs(y) + eps); % Normalize the output
sound(y, Fs)
audiowrite('effectOutput.wav', y, Fs);

spectrograms(s,x,y,e,Fs);

%% experiment: i divide the bandwidth in two regions with different lpc orders
P_x1 = 20;
P_x2 = 90;
Fcut = 4000;
[y,e] = LPCsplitted(s, x, w, Fs, wLength, hopSize, P_s, P_x1, P_x2, Fcut, Nfft);
y = y ./ max(abs(y) + eps); % Normalize the output
sound(y, Fs)
audiowrite('lpcsplitOutput.wav', y, Fs);

spectrograms(s,x,y,e,Fs);