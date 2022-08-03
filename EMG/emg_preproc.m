%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 3 Mar 2021. Last revision: 17 Mar 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mysegment, emg_clean] = emg_preproc (sampleinfo, ind, emgdata, mvc, showplot)
% This function performs preprocessing of EMG and MVC data, with inputs: 
%    sampleinfo: an array containing Fs and sf
%    ind: event markers, indices
%    emgdata: vector containing recording of one EMG sensor
%    mvc: a flag to indicate if this is mvcdata
%    showplot: a flag to show figures or not
% It then returns OUTPUT:
%    mysegment: the preprocessed data we want
%    emg_clean: filtered EMG data

% (1) EMG Data timing information, taken from 'sampleinfo' input
Fs = sampleinfo(1);       % retrieve the sampling frequency
sf = sampleinfo(2);       % Sampling factor (only for Qualisys)
T  = 1/Fs;                % Sampling period       
L  = length(emgdata);     % Length of signal
t  = (0:L-1)*T;           % Time vector
eachtrial = [];           % initialize empty array for the output

% Note: for some subjects, the sampling rate was modified...
if Fs==1000,  Fs=2000;  end

% Check for NaN (undefined), then fill the gap by using interpolation
emgdata = fillmissing(emgdata, 'linear');   % >>> USEFUL function I found!!!
mystart = 5; myend = length(emgdata)-5;


%% (2) Load event marker data (index, only for Qualisys)
% Event marker data = take Frame (samples) not the Time (sec)!
if ~isempty(ind)
    % If not empty or there are > 2 event markers, then it's perturbation trials
    mvcdata = false;
    % Load event marker data (index, only for Qualisys)
    ind = struct2cell(ind);   % Convert struct type to cell first...
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
    mvcdata = mvc;    % this indicates it's for MVC
end


%% (3) Filter the EMG data with a bandpass filter!
lowf  = sampleinfo(3);               % Low freq cutoff (Hz)
highf = sampleinfo(4);               % High freq cutoff (Hz)
cutfreq = [lowf highf]/(Fs/2);       % provide cutoff freq for filter
[B,A] = butter(3,cutfreq);           % design butterworth filter type
emg_clean0 = filtfilt(B,A,emgdata);  % do filtering now!! 

% Recording near the trunk and shoulder may suffer from cardiac artifacts. This is 
% interesting but troublesome. Frequency content of normal PQRST is 8-50 Hz (Larisa GT, 2015). 
% You can either do a bandpass or perform Notch filter on 15Hz (heart rate)
%w0 = 20/(Fs/2);   % for 30Hz with a data sample rate of Fs  
%[num,den] = iirnotch(w0,w0/35);             % design notch filter
%emg_clean = filtfilt(num,den,emg_clean0);   % then do filtering 
emg_clean = emg_clean0;



%% (4) Perform FFT to see the freq content
if(showplot)
    figure(1)
    subplot(2,1,1)
    plot(emgdata,'k'); hold on;
    plot(emg_clean,'y'); 
    title('Raw EMG (black) and filtered EMG data (yellow)');
    % Draw multiple lines! Very nice, I found this feature online ---
    if ~mvcdata, arrayfun(@(a)xline(a),ind); end;
    hold off;
end

Y = fft(emgdata);     % Perform FFT on the original data
P2 = abs(Y/L);        % Take the magnitude of FFT
P1 = P2(1:L/2+1);     % FFT is symmetrical, take the half only!
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

if(showplot)
%     figure(2)
%     subplot(1,2,1)
%     plot(f,P1,'k');  title('Before bandpass')
%     ylabel('|Y(f)|');  xlabel('f (Hz)');
%     axis([0 500 0 500])

% Again, FFT for the filtered or cleaned emg
% Y = fft(emg_clean0);  % Perform FFT on bandpassed data
% P2 = abs(Y/L); 
% P1 = P2(1:L/2+1); 
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% subplot(1,2,2);  plot(f,P1,'b');  
% ylabel('|Y(f)|');  xlabel('f (Hz)');
% axis([0 500 0 500])
% title('After bandpass (blue), notch (yellow)')
% hold on;
% Y = fft(emg_clean);   % Perform FFT on the final cleaned data (after notch filter)
% P2 = abs(Y/L); 
% P1 = P2(1:L/2+1); 
% P1(2:end-1) = 2*P1(2:end-1);
% f = Fs*(0:(L/2))/L;
% subplot(1,2,2);  plot(f,P1,'y');
% hold off;
end


%% (5) Next pre-processing step: Hilbert transform
% Take the absolute value of the filtered emg data!
emg_clean_abs = abs(emg_clean);
% Envelope detector using Hilbert transform, connecting peaks for every nn = 50 samples
nn = 50;
emg_clean_envelope = envelope(emg_clean_abs,nn,'peak');

if(showplot)
     figure(1)
     subplot(2,1,2)
     plot(emg_clean_abs,'y'); hold on;
     plot(emg_clean_envelope,'r'); 
     if ~mvcdata, arrayfun(@(a)xline(a),ind); end;
     title('Linear envelope (red) and the rectified EMG data (yellow)')
     hold off;
end

% Now we define start/stop for MVC data only (not for trial data)! Use manual
% clicking to define 3 starts + 3 stops for 3 MVC trials. We just need X-axis info only!
%if mvcdata 
%    [ind, ~] = ginput(6);
%    ind = round(ind);   % The ginput returns decimals, so we round the numbers
%end



%% (6) Now transfer all trial segments of interest into a new variable. 
% Note: For trial data, the ind is not empty but contains event marker. Use it to 
% define the start/stop of the window or segment we are interested in; which coverts 
% 250 msec (500 frames) before and 1000 msec (2000 frames) after the event trigger.
% For MVC data, the start/stop was defined manually.
if mvcdata  
    mystart = 1; %ind(1:2:5);         % Obtain the start points from the ind array
    myend   = length(emgdata)-1;%ind(2:2:6);         % And... the end points 
else
    if(~isempty(ind))
        mystart = ind - 0.25*Fs;  % S01-S05: 0.05; S06-S11: 0.55
        myend   = ind + 1.0*Fs;   % S01-S05: 0.5; S06-S11: 0.00
    end
end

% I keep adding new trial segments of interest into a big 2D matrix of data 
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
        piece = emg_clean_envelope(mystart(j):myend(j));   % This is for perturbation!
    end
    eachtrial(:,j) = piece;  
end
%size(eachtrial)

% if mvcdata
     mysegment = eachtrial;
% else
%     mysegment = mean(eachtrial,2);
% end


end