
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code cuts kinematic trajectories. In addition to the
% textfile from the Visual3D, I also load marker trajectories from QTM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Define path and filename.
mypath = 'D:\extracted angle\';
mysubj = 'SN015';
mypath = strcat(mypath,mysubj);

% First, load the angle trajectories file......
myfile = '_towel_L';
ntrial = 10;   % <<< We define number of trials in a task
mydata = readtable(strcat(mypath,'\',mysubj,myfile,'.txt'), 'HeaderLines',1);  

mydata.Properties.VariableNames={'trial','LElbowX','LElbowY','LElbowZ',...
    'LShoulderX','LShoulderY','LShoulderZ','LWristX','LWristY',...
    'LWristZ','HeadX','HeadY','HeadZ','PelvisX',...	
	'PelvisY','PelvisZ','RElbowX','RElbowY','RElbowZ','RShoulderX',...	
	'RShoulderY','RShoulderZ','RWristX','RWristY','RWristZ',...
    'TxX','TxY','TxZ'};


if(1)
%%%%%%%%%%%%%  ANALOG SIGNAL FROM THE TABLE LOADCELLS  %%%%%%%%%%%%%%%%%%%%
analogfile = strcat(myfile,'_loadcells.txt');
loadcells  = readtable(strcat(mypath,'\',mysubj,analogfile), ...
                        'Delimiter','\t','HeaderLines', 0, 'ReadVariableNames', true);
loadcells = loadcells{:,[1,2,4]};  % Take only relevant columns

% Define and apply Butterworth filter for analog signals....
%Fc = 10; Fs_kin = 2000;     % 10 Hz cutoff frequency, 2000 Hz sampling rate.
%[B,A]=butter(2,Fc/(Fs_kin/2));
%loadcells = filtfilt(B,A,loadcells);

% Take note that the analog signal has higher sampling rate!
% Time-normalize to be the same as the kinematic data.................
loadcell2 = interp1(1:length(loadcells), ...  % Size of original data
                    loadcells, ...            % Name of original data
                    linspace(1,length(loadcells),size(mydata,1)));   % Downsizing!
end


if(0)
% Load the right hand marker trajectories.......
myfile  = 'SN005_0008_towel_L02_Rhand.txt';
myRhand = readtable(strcat(mypath,myfile),'Delimiter','\t','HeaderLines',0);
myRhand.Rhand_path = sqrt(myRhand.X.^2 + myRhand.Y.^2 + myRhand.Z.^2);

% Combine all table into one.
T = [mydata, table(myRhand.Y,  'VariableNames',{'Rhand_path'}), ... 
          table(myRshoul.Rshoul_path,'VariableNames',{'Rshoul_path'}), ...
          table(myLhand.Y,  'VariableNames',{'Lhand_path'}) , ...
          table(myLshoul.Lshoul_path,'VariableNames',{'Lshoul_path'}) ];
end
      



% Which joint angle is the cleanest? You may want to plot other
% trajectories too as comparison!!
plot(zscore(mydata.RElbowX), 'k'); hold on;
plot(zscore(mydata.LElbowX), 'b');
plot(zscore(loadcell2(:,2)), '.');
plot(zscore(loadcell2(:,3)), '.');


%plot(zscore(mydata.Index_pos_Z),'r'); hold on; % For HAND only!

%%%%%%%%%%%%%%%% I found this code online %%%%%%%%%%%%%%%%
ntrial = 10;
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
        line([x x],[-3 3],'Color','g');
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
     start  = X(1:nevents:length(X));  start = start - 50;
     middle = X(2:nevents:length(X));
     stop   = X(3:nevents:length(X));  stop = stop + 50;
     
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

% It's wise to plot the curves again with the segmentation indices for reference!!
figure(3)
  plot(zscore(mydata.RElbowX), 'k'); hold on;
  plot(zscore(mydata.LElbowX), 'b');
  xlim([-100,3000]);
  line([X X],[-3 3],'Color','g');
  savefig( strcat(mypath,'\',mysubj,myfile) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save the cut trajectories to a textfile. Ask whether we want to save it.
fname = strcat('cut',myfile,'.txt');
fname = strcat(mypath,'\',fname)   % the same filepath
sprintf("Saving the segmented files")
writetable(all_tnorm, fname)
save(strcat(mypath,'\',mysubj,myfile));

