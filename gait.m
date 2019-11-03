
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code examines grasping pattern of GAIT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all

mypath = 'D:\extracted angle\Normal gait\'
myfile = 'SN008 LL angles 10m.csv';

mydata = readtable(strcat(mypath,myfile), 'HeaderLines',2);

colNames = {'trial','ankle_angle_l','ankle_angle_lY','ankle_angle_lZ',...
         'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
         'knee_angle_l','knee_angle_lY','knee_angle_lZ',...
         'pelvis_tilt','pelvis_list','pelvis_rotation',...
         'ankle_angle_r','ankle_angle_rY','ankle_angle_rZ',...
         'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
         'knee_angle_r','knee_angle_rY','knee_angle_rZ'};       
        
[row, col] = size(mydata);
multiple = (col-1)/21;

myarray = mydata{:,2:col};

myarray2 = [myarray(:, 1:21); 
            myarray(:,22:42); 
            myarray(:,43:63);
            myarray(:,64:84)];

% Remove all NaN rows in the column dimension!
myarray2(any(isnan(myarray2), 2), :) = [];

% IMPORTANT: Let's cut based on the ANKLE flexion/extension angles
plot( zscore(myarray2(:,1)) ); 
hold on; 
plot( zscore(myarray2(:,4)) ); 


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
start  = X(1:2:length(X));
stop   = X(2:2:length(X));
     
% Obtain the subset of the original data from start to stop ONLY...
for k = 1 : length(start)
    tnorm_cut = [];  % initialize
    cutdata = myarray2(start(k):stop(k),:);   % we move back 5 frames!
    npoints = size(cutdata,1);
    for m = 1: size(cutdata,2)-1
        % Perform time normalization using interp1 and linspace setting!
        tscaled  = interp1(1:npoints, cutdata(:,m), linspace(1,npoints,normt));
        tnorm_cut = [tnorm_cut tscaled'];
    end
    tnorm_cut = [repmat(k,[normt 1]) tnorm_cut];  % Remember to add trial #
    tnorm_cut = [tnorm_cut repmat(1,[normt 1])];  % And add phase-1 too
    all_tnorm = [all_tnorm; tnorm_cut];
end

% Convert the segmented array back to a table format!
all_tnorm = array2table(all_tnorm,'VariableNames',colNames );

plot( all_tnorm{:,2} ); 
plot(all_tnorm.ankle_angle_r,all_tnorm.hip_flexion_r,'k.'); axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Save the cut trajectories to a textfile. Ask whether we want to save it.
newStr = extractBetween(myfile,1,5);
fname = strcat(mypath,'new_',newStr{1});   % the same filepath
sprintf("Saving the segmented files");
writetable(all_tnorm, fname)


