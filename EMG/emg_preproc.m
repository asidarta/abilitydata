%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 3 Mar 2021. Last revision: 17 Mar 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mysegment, emg_clean] = emg_preproc (sampleinfo, ind, emgdata)
% This function performs preprocessing of EMG and MVC data, with inputs: 
%   sampleinfo: an array containing Fs and sf
%   ind: event markers, indices
%   emgdata: vector containing recording of one EMG sensor
% It then returns OUTPUT:
%   mysegment containing the EMG segment we want

% (1) EMG Data timing information, taken from 'sampleinfo' input
Fs = sampleinfo(1);       % retrieve the sampling frequency
sf = sampleinfo(2);       % Sampling factor (only for Qualisys)
T  = 1/Fs;                % Sampling period       
L  = length(emgdata);     % Length of signal
t  = (0:L-1)*T;           % Time vector
mysegment = [];           % initialize empty array for the output

% Note: for some subjects, the sampling rate was modified...
if Fs==1000,  Fs=2000;  end


%% (2) Load event marker data (index, only for Qualisys)
% Event marker data = take Frame (samples) not the Time (sec)!
if ~isempty(ind)
    mvcdata = false;   % Is this for MVC trials? 
    % Load event marker data (index, only for Qualisys)
    ind = struct2cell(ind);   % Convert struct type to cell.
    % Now we extract the Frame (samples) not the Time (sec)!
    ind = cell2mat(ind(2,:));
    % Multiply with sampling factor. Round up the value!
    ind = sf * ceil(ind);

    % Note: Due to EMG recording problem, sometimes there are missing data, 
    % but the event trigger was intact!
    while ( ind(length(ind))+1000 > length(emgdata))
        ind(length(ind)) = [];    % throw away from the last number!
    end
else
    mvcdata = true;  % this indicates it's for MVC
end


%% (3) Filter the EMG data with a bandpass filter!
lowf  = 40;         % Low freq cutoff (Hz)
highf = 150;        % High freq cutoff (Hz)
cutfreq = [lowf highf]/(Fs/2);      % provide cutoff freq for filter
[B,A] = butter(3,cutfreq);          % design butterworth filter type
emg_clean = filtfilt(B,A,emgdata);  % do filtering now!! 

figure(1)
plot(emgdata,'k'); hold on;
plot(emg_clean,'y'); 
title('Raw EMG (black) and filtered EMG data (yellow)');
% Draw multiple lines! Very nice, I found this feature online ---
if ~mvcdata, arrayfun(@(a)xline(a),ind); end;
hold off;

%% Perform FFT to see the freq content
Y = fft(emgdata);     % Perform FFT on time-series; contains complex numbers!
P2 = abs(Y/L);        % Take the magnitude of FFT
P1 = P2(1:L/2+1);     % FFT is symmetrical, take the half only!
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

figure(2)
subplot(1,2,1)
plot(f,P1);  title('Before bandpass')
ylabel('|Y(f)|');  xlabel('f (Hz)');
axis([0 500 0 1])

% (4) Again, FFT for the filtered emg
Y = fft(emg_clean);   % Perform FFT on time-series; contains complex numbers!
P2 = abs(Y/L);        % Take the magnitude of FFT
P1 = P2(1:L/2+1);     % FFT is symmetrical, take the half only!
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
subplot(1,2,2)
plot(f,P1);  title('After bandpass')
ylabel('|Y(f)|');  xlabel('f (Hz)');
axis([0 500 0 1])


%% (5) Next pre-processing steps
% Take the absolute value of the filtered emg data!
emg_clean_abs = abs(emg_clean);
% Envelope detector using Hilbert transform, connecting only the peaks!
emg_clean_envelope = envelope(emg_clean_abs,50,'peak');
figure(3)
plot(emg_clean_abs,'y'); hold on;
plot(emg_clean_envelope,'r'); 

% Now we define start/stop for MVC data only (not for trial data)! Use manual
% clicking to define 3 starts + 3 stops for 3 MVC trials. We just need X-axis info only!
if mvcdata, [ind, ~] = ginput(6); end
ind = round(ind);   % The ginput returns decimals, so we round the numbers

% Excitation detector....
baseline = mean(emg_clean_envelope(1:50));
baseline = mean(emg_clean_envelope(length(emgdata)-1:-1:length(emgdata)-50));
for i = 1:length(emg_clean_envelope)
    if (emg_clean_envelope(i) > 2*baseline)
        excite(i) = emg_clean_envelope(i);
    else
        excite(i) = NaN;
    end
end

%plot(excite, 'g.');
% Draw multiple lines! Very nice, I found this feature online!
arrayfun(@(a)xline(a),ind);
hold off;  


%% (6) Now transfer all the segments of interest into a new variable. 
% Note: For trial data, the ind is not empty but contains event marker. Use it to 
% define the start/stop of the window or segment we are interested in; which covers 
% 5 msec (100 frames) before and 500 msec (1000 frames) after the event trigger.
% For MVC data, the start/stop was defined manually.
if mvcdata  
    mystart = ind(1:2:5);     % Obtain the start points from the ind array
    myend   = ind(2:2:6);     % And... the end points 
else
    mystart = ind - 0.55*Fs;  % S01-S05: 0.05; S06-S11: 0.55
    myend   = ind + 0*Fs;     % S01-S05: 0.5; S06-S11: 0.00
end

% I keep adding new segments of interest into a big 2D matrix of data 
% belonging to the SAME muscle so I can take the average later. For manual
% selection (MVC) we should do time-interpolation to 2000 data points, so that 
% it'll be the same for all trials.
for j = 1:length(mystart)
    if mvcdata
        npoints = myend(j)-mystart(j)+1;  stretch = 2000;
        piece = interp1(1:npoints, emg_clean_envelope(mystart(j):myend(j)), ...
                    linspace(1,npoints,stretch));   % interpolate
        piece = piece';  % transpose (flip!)
    else
        piece = emg_clean_envelope(mystart(j):myend(j));
    end
    mysegment = [mysegment, piece];
end

% NOTE: Careful which processed EMG to use (filtered, absolute, or envelope)

%% (7) Plot segments of interest of each muscle
%          figure(4)
%          %plot(mysegment);
%          options.handle     = 4;
%          options.error      = 'sem';
%          options.color_area = [243 169 114]./255;    % Orange theme
%          options.color_line = [236 112  22]./255;
%          options.x_axis = -0.05*Fs: 2000/Fs: 0.5*Fs;
%          plot_areaerrorbar(mysegment',options);
%          ylim([-60,60]); xline(0);  % Indicate trigger line!
%          xlabel('Time (msec)');
%          ylabel('Amplitude (filtered)')
%          title(['Average EMG reading for Electrode-' num2str(mmm)]);

end