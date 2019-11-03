
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code cuts kinematic trajectories. In addition to the
% textfile from the Visual3D specifically for GAIT.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Define path and filename.
%mypath = 'C:\Users\ananda.sidarta\Documents\extracted angle\SN6\';
mypath = 'D:\extracted angle\Stroke gait\';

% First, load the angle trajectories file......
%myfile = 'walk_RL4.csv';
myfile = 'S98.csv';

mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  
colNames = mydata.Properties.VariableNames;
%mydata.Properties.VariableNames{'Var1'}='trial';
%mydata.Properties.VariableNames={'trial','RhipX','RhipY','RhipZ','Rknee','RankleX','RankleY',...
%     'RankleZ','LhipX','LhipY','LhipZ','Lknee','LankleX','LankleY','LankleZ'};

% Load the left shoulder marker trajectories.......
%myfile  = 'SN005_0008_towel_L02_Lshoulder.txt';
%myLshoul= readtable(strcat(mypath,myfile),'Delimiter','\t','HeaderLines',0);  
%myLshoul.Lshoul_path = sqrt(myLshoul.X.^2 + myLshoul.Y.^2 + myLshoul.Z.^2);

T=mydata;
      
% Which joint angle is the cleanest? For gait, hip flexion/extension of R/L
% are alternating, so it is the clearest!
plot(zscore(T.hip_flexion_l)); hold on;
plot(zscore(T.ankle_angle_l));
%plot(zscore(T.RhipX));


%%%%%%%%%%%%%%%% I found this code online %%%%%%%%%%%%%%%%
X = []; Y = [];
hold on; xlim([0,1060]); 
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
normt     = 100; % How many time-points you want to scale down to?
cutdata   = [];
all_tnorm = [];
tscaled   = [];
     
% Obtain the subset of the original data from start to stop ONLY...
for k = 1 : (length(X)-1)
   tnorm_cut = [];  % initialize
   cutdata = mydata(X(k):X(k+1)-1,:);
   npoints = size(cutdata,1);
   for m = 1:size(mydata,2)-1
        % Perform time normalization using interp1 and linspace setting!
        tscaled  = interp1(1:npoints, cutdata{:,m+1}, linspace(1,npoints,normt));
        tnorm_cut = [tnorm_cut tscaled'];
   end
   tnorm_cut = [repmat(k,[normt 1]) tnorm_cut];  % Remember to add trial #
   all_tnorm = [all_tnorm; tnorm_cut];
end
    
% Convert the segmented array back to a table format!
all_tnorm = array2table(all_tnorm,'VariableNames',colNames );

%plot(all_tnorm.ankle_angle_r,all_tnorm.hip_flexion_r);
plot(all_tnorm.ankle_angle_l); hold on; plot(all_tnorm.hip_flexion_r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Save the cut trajectories to a textfile. Ask whether we want to save it.
fname = strcat(mypath,'new1',myfile);   % the same filepath
sprintf("Saving the segmented files");
writetable(all_tnorm, fname)

