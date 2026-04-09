function [y,e] = LPCsplitted(s, x, w, Fs, wLength, hopSize, P_s, P_x1, P_x2, Fcut, Nfft)
% INPUTS:
%   s: input signal speech
%   x: input signal sound
%   w: window 
%   Fs: frequncy of sampling
%   wLength: window length
%   hopSize: hop size of the windowing
%   P_s: p order used for LPC of the speech
%   P_x: p order used for LPC of the sound
%   Nfft: nfft for the plot of the frame example
%
% OUTPUTS:
%   y: the talking instrument
%
%   Detailed explanation goes here

sigLength = length(s);
nFrames = floor((sigLength-wLength)/hopSize) + 1;
outLength = (nFrames - 1) * hopSize + wLength;

% Output signal
y = zeros(1, outLength);
e = zeros(1, outLength);


for i = 1:nFrames % Compute for each frame
    
    idxFrame = (i-1) * hopSize + 1 : (i-1) * hopSize + wLength; % Index of the frame
    
    sFrame = s(idxFrame) .* w; % Windowing each signal
    xFrame = x(idxFrame) .* w;
    
    %--- FIND SHAPING FILTER: work with s ---%
    % Observe that we implement LPC using the the Levison-Durbin recursion. 
    % Such algorithm has a lighter computational cost

    [r_s, rlags_s] = xcorr(sFrame); % Auto-correlation computation of the frame

    % Extract the required autocorrelation coefficients
    central_lag_s = ceil(length(r_s) / 2);

    % Apply Levinson-Durbin recursion
    a_s = levinson(r_s(central_lag_s:(central_lag_s+P_s)), P_s);

    % Prediction error (not needed)
    e_s = conv(sFrame, a_s, 'same'); 

     %--- FIND EXCITATION: work with x ---%
    h_low = fir1(100, Fcut/(Fs/2), 'low');
    h_high = fir1(100, Fcut/(Fs/2), 'high');
    xLow = filter(h_low, 1, xFrame);
    xHigh = filter(h_high, 1, xFrame);

    [r_xLow, rlags_xL] = xcorr(xLow);
    [r_xHigh, rlags_xH] = xcorr(xHigh);

    central_lag_xL = ceil(length(r_xLow) / 2);

    a_x1 = levinson(r_xLow(central_lag_xL:(central_lag_xL+P_x1)), P_x1);
    e_x1 = conv(xLow, a_x1, 'same'); 

    central_lag_xH = ceil(length(r_xHigh) / 2);

    a_x2 = levinson(r_xHigh(central_lag_xH:(central_lag_xH+P_x2)), P_x2);
    e_x2 = conv(xHigh, a_x2, 'same'); 

    e_x = e_x1 + e_x2;

    %--- RECONSTRUCTION: combine s' filter with x' excitation ---%
    yFrame = filter(1, a_s, e_x);

    % OLA reconstruction
    y(idxFrame) = y(idxFrame) + (yFrame)';
    e(idxFrame) = e(idxFrame) + (e_x)';
    
    % For one frame let's plot spectrograms
    if i == 50
         %Nfft 2^(ceil(log2( L + M -1))); ok non centra forse

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


end