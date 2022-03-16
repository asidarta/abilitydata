
function [freq,xdft,psdx] = fft_matlab(mydata, Fs)
    % The purpose of this function is to perform FFT on the time-varying input signal
    % Inputs:  - time-series (mydata)
    %          - sampling freq (Fs)
    % Outputs: - Fourier magnitude
    %          - frequency 
    %          - power spectral density

N = length(mydata);               % Length of input data
xdft = fft(mydata);               % Perform FFT, output contains complex numbers
xdft = xdft(1:N/2+1);             % FFT is symmetrical, take the half of the content only!
xdft(2:end-1) = 2*xdft(2:end-1);  % But in order to conserve the total power, multiply by 2!
xdft = abs(xdft);                 % Get the magnitude only...

freq = 0:Fs/N:Fs/2;               % Define the discrete freq range
psdx = (1/(Fs*N))*abs(xdft).^2;   % Power spectral density

