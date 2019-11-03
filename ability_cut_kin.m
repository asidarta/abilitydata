
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code cuts kinematic trajectories. In addition to the
% textfile from the Visual3D, I also load marker trajectories from QTM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Define path and filename.
mypath = 'D:\extracted angle\SN010\';

% First, load the angle trajectories file......
myfile = 'SN010_0005_towel_L01.txt';
ntrial = 10;   % <<< We define number of trials in a task
mydata = readtable(strcat(mypath,myfile), 'HeaderLines',2);  
%mydata.Properties.VariableNames{'Var1'}='trial';
mydata.Properties.VariableNames={'trial','LElbowX','LElbowY','LElbowZ',...
    'LShoulderX','LShoulderY','LShoulderZ','LWristX','LWristY',...
    'LWristZ','HeadX','HeadY','HeadZ','PelvisX',...	
	'PelvisY','PelvisZ','RElbowX','RElbowY','RElbowZ','RShoulderX',...	
	'RShoulderY','RShoulderZ','RWristX','RWristY','RWristZ',...
    'TxX','TxY','TxZ'};
if(0)
% Load the right hand marker trajectories.......
myfile  = 'SN005_0008_towel_L02_Rhand.txt';
myRhand = readtable(strcat(mypath,myfile),'Delimiter','\t','HeaderLines',0);
myRhand.Rhand_path = sqrt(myRhand.X.^2 + myRhand.Y.^2 + myRhand.Z.^2);

% Load the right shoulder marker trajectories.......
myfile  = 'SN005_0008_towel_L02_Rshoulder.txt';
myRshoul= readtable(strcat(mypath,myfile),'Delimiter','\t','HeaderLines',0);  
myRshoul.Rshoul_path = sqrt(myRshoul.X.^2 + myRshoul.Y.^2 + myRshoul.Z.^2);

% Load the left hand marker trajectories.......
myfile  = 'SN005_0008_towel_L02_Lhand.txt';
myLhand = readtable(strcat(mypath,myfile),'Delimiter','\t','HeaderLines',0); 
myLhand.Lhand_path = sqrt(myLhand.X.^2 + myLhand.Y.^2 + myLhand.Z.^2);

% Load the left shoulder marker trajectories.......
myfile  = 'SN005_0008_towel_L02_Lshoulder.txt';
myLshoul= readtable(strcat(mypath,myfile),'Delimiter','\t','HeaderLines',0);  
myLshoul.Lshoul_path = sqrt(myLshoul.X.^2 + myLshoul.Y.^2 + myLshoul.Z.^2);

% Combine all table into one.
T = [mydata, table(myRhand.Y,  'VariableNames',{'Rhand_path'}), ... 
          table(myRshoul.Rshoul_path,'VariableNames',{'Rshoul_path'}), ...
          table(myLhand.Y,  'VariableNames',{'Lhand_path'}) , ...
          table(myLshoul.Lshoul_path,'VariableNames',{'Lshoul_path'}) ];
end
      
if(0)     
%%%%%%%%%%%% Do mini freq domain analyses using CSD %%%%%%%%%%%%%%%%
var1 = mydata{1:3000,9}; 
Fs = 200;                % Sampling frequency                    
T  = 1/Fs;               % Sampling period       
L  = length(var1);       % Length of signal
t  = (0:L-1)*T;          % Time vector
Y  = fft(squeeze(var1)); % Do FFT on the amplitude data
YL = abs(Y/L);
P1 = YL(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;      % Freq-range for X-axis, folded.
plot(f,P1)               % Plot the single-sided spectrum 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)'); ylabel('|P1(f)|')
xlim([0,15]); ylim([-2 20]);

hold on;

var1 = mydata{1:3000,18}; 
Fs = 200;                % Sampling frequency                    
T  = 1/Fs;               % Sampling period       
L  = length(var1);       % Length of signal
t  = (0:L-1)*T;          % Time vector
Y  = fft(squeeze(var1)); % Do FFT on the amplitude data
YL = abs(Y/L);
P2 = YL(1:L/2+1);
P2(2:end-1) = 2*P2(2:end-1);
f = Fs*(0:(L/2))/L;      % Freq-range for X-axis, folded.
plot(f,P2)               % Plot the single-sided spectrum 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)'); ylabel('|P2(f)|')
xlim([0,15]); ylim([-2 20]);

%%%% Cross spectrum analysis, measuring coherence %%%%%
numer = (abs(cpsd(P1,P2))).^2;
denom = abs(cpsd(P1,P1)).* abs(cpsd(P2,P2));
coherence = numer ./ denom; 
plot(coherence, 'b'); hold on;
%mscohere(P1,P2);  % or, using the function
end      

      
% Which joint angle is the cleanest? You may want to plot other
% trajectories too as comparison!!
plot(zscore(mydata.LElbowX)); hold on;
plot(zscore(mydata.RElbowX));
%plot(zscore(abcde.Z));
%plot(zscore(T.Lhand_path)); 
%plot(diff(T.Rhand_path)); hold on;
%plot(zscore(T.Rhand_path));
%plot(zscore(T.Lshoul_path));
%[x y] = ginput(ntrial*2);  

%plot(zscore(mydata.Index_pos_Z),'r'); hold on; % For HAND only!

%%%%%%%%%%%%%%%% I found this code online %%%%%%%%%%%%%%%%
ntrial = 10;%10
xlim([-100,5000]); X = []; Y = [];
title('Press Enter to Quit')
while 1
    [x,y,b] = ginput(1); % Get only 1 click
    if isempty(b)  % Press <ENTER> key to bail out
        break;
    elseif b==60
        ax = axis; width=ax(2)-ax(1); %height=ax(4)-ax(3);
        xlim([ax(1)-1000, ax(2)-1000]);
    elseif b==62
        ax = axis; width=ax(2)-ax(1); %sheight=ax(4)-ax(3);
        xlim([ax(1)+1000, ax(2)+1000]);
    else
        X=[X;x];
        line([x x],[-2 2],'Color','g');
        %Y=[Y;y];
    end;
end
close all;



%%%%%%%%%%%%%%%%%%%% PARSING THE ORIGINAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtain the indices in time (X-axis) for segmentation..
X = round(X);
normt     = 200; % How many time-points you want to scale down to?
cutdata   = [];
all_tnorm = [];
tscaled   = [];

% We define how many events we are interested in....
nevents = input('How many events per trial (2 or 3)?  ');

% NOTE: For the purpose of the folding towel task, I consider it to have
% 2 phases. Phase-1 involves bilateral arm, phase-2 is the lateral fold.
% Now estimated start, middle, and stop by visual inspection!!
if nevents == 3
     start  = X(1:nevents:length(X));
     middle = X(2:nevents:length(X));
     stop   = X(3:nevents:length(X));
     
     % Obtain the subset of the original data from start to stop ONLY...
     for k = 1:ntrial
        %%%% Phase 1 ---------------------------------------------------------
        tnorm_cut = [];  % initialize
        cutdata = mydata(start(k):middle(k),:);
        npoints = size(cutdata,1);
        for m = 1:size(mydata,2)-1
            % Perform time normalization using interp1 and linspace setting!
            tscaled  = interp1(1:npoints, cutdata{:,m+1}, linspace(1,npoints,normt));
            tnorm_cut = [tnorm_cut tscaled'];
        end
        tnorm_cut = [repmat(k,[normt 1]) tnorm_cut];  % Remember to add trial num
        tnorm_cut = [tnorm_cut repmat(1,[normt 1])];  % And add phase-1 too
        all_tnorm = [all_tnorm; tnorm_cut];
        %%%% Phase 2 ---------------------------------------------------------
        tnorm_cut = [];  % initialize
        cutdata = mydata(middle(k):stop(k),:);
        npoints = size(cutdata,1);
        for m = 1:size(mydata,2)-1
            % Perform time normalization using interp1 and linspace setting!
            tscaled  = interp1(1:npoints, cutdata{:,m+1}, linspace(1,npoints,normt));
            tnorm_cut = [tnorm_cut tscaled'];
        end
        tnorm_cut = [repmat(k,[normt 1]) tnorm_cut];
        tnorm_cut = [tnorm_cut repmat(2,[normt 1])];  % And this is phase-2
        all_tnorm = [all_tnorm; tnorm_cut];
    end

% Generally, we really just want START and STOP of the arm/hand!
elseif nevents == 2
     start  = X(1:nevents:length(X));
     stop   = X(2:nevents:length(X));
     
     % Obtain the subset of the original data from start to stop ONLY...
     for k = 1:ntrial
        tnorm_cut = [];  % initialize
        cutdata = mydata(start(k):stop(k),:);
        npoints = size(cutdata,1);
        for m = 1:size(mydata,2)-1
            % Perform time normalization using interp1 and linspace setting!
            tscaled  = interp1(1:npoints, cutdata{:,m+1}, linspace(1,npoints,normt));
            tnorm_cut = [tnorm_cut tscaled'];
        end
        tnorm_cut = [repmat(k,[normt 1]) tnorm_cut];  % Remember to add trial #
        tnorm_cut = [tnorm_cut repmat(1,[normt 1])];  % And add phase-1 too
        all_tnorm = [all_tnorm; tnorm_cut];
     end
end
    
% Convert the segmented array back to a table format!
all_tnorm = array2table(all_tnorm,'VariableNames', ...
                        [mydata.Properties.VariableNames 'phase' ] );

%plot(all_tnorm.Index_pos_Z);
plot(all_tnorm.LElbowX)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Save the cut trajectories to a textfile. Ask whether we want to save it.
fname = input('What is your filename you want to save? ','s');
fname = strcat(mypath,fname,'.txt');   % the same filepath
sprintf("Saving the segmented files")
writetable(all_tnorm, fname)

