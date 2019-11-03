
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code examines grasping pattern of our upper limb tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; %close all

% Define path and filename.
mypath = 'D:\extracted angle\lateral\'

myfile = 'SN67_lateral_R_cube.txt';
mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  
%mydata = mydata([1:1800],:);
nframe = 200; %size(mydata,1)/ntrial;
ntrial = size(mydata,1)/nframe;         % <<< define number of trials per task

% Once loaded, we first compute the Euclidean distance between two finger
% markers to estimate the aperture size!!
aperture = sqrt( (mydata.Thumb_pos_X-mydata.Index_pos_X).^2 + ...
                 (mydata.Thumb_pos_Y-mydata.Index_pos_Y).^2 + ...
                 (mydata.Thumb_pos_Z-mydata.Index_pos_Z).^2 );
mydata.aperture = aperture;
%plot(aperture); 
%plot(mydata.Index_pos_X,'r-'); hold on; plot(mydata.Thumb_pos_X,'b-');
if(0)
plot3(  reshape(mydata.Thumb_pos_X, [200,10]), ...
        reshape(mydata.Thumb_pos_Y, [200,10]), ...
        reshape(mydata.Thumb_pos_Z, [200,10])); hold on;
    
plot3(  reshape(mydata.Index_pos_X, [200,10]), ...
        reshape(mydata.Index_pos_Y, [200,10]), ...
        reshape(mydata.Index_pos_Z, [200,10]), '.'); hold off;

plot3( (mydata.Thumb_pos_X-mydata.Index_pos_X), ...
       (mydata.Thumb_pos_Y-mydata.Index_pos_Y), ...
       (mydata.Thumb_pos_Z-mydata.Index_pos_Z) );
end
   
% We find mean trajectories and standard deviation of each trajectory!
avg = []; stdev = [];
for p = 1 : size(mydata,2)
    temp  = mydata{:,p};
    temp  = reshape(temp, [nframe, ntrial]);
    tempz = zscore(temp);
    %if p == 11
    %    plot(tempz); hold on;
    %end
    avg(:,p)   = mean(temp,2);
    stdev(:,p) = std(temp,0,2);
end

% Plot the graph with shading
j = 11;  % Which column you want to plot??
mymean = avg(:,j); myse = stdev(:,j);
lo = mymean - myse;
hi = mymean + myse;
nsample = length(lo);
mycolor = 'b';   % blue for Left -------------
hl = line(1:length(lo),mymean,'color',mycolor,'LineStyle','-.');
hp = patch([(1:nsample)'; (nsample:-1:1)'; 1], ...
           [lo; hi(nsample:-1:1); lo(1)], mycolor);
set(hp, 'facecolor', mycolor, 'edgecolor', 'none', 'FaceAlpha',.1);
title('Mean kinematics trajectories');
xlabel('Normalized time-points'); ylabel('Normalized distance');


hold on;
if(0)
j = 4;  % Which column you want to plot??
mymean = avg(:,j); myse = stdev(:,j);
lo = mymean - myse;
hi = mymean + myse;
nsample = length(lo);
mycolor = 'm';   % blue for Left -------------
hl = line(1:length(lo),mymean,'color', mycolor,'LineStyle','--');
hp = patch([(1:nsample)'; (nsample:-1:1)'; 1], ...
           [lo; hi(nsample:-1:1); lo(1)], mycolor);
set(hp, 'facecolor', mycolor, 'edgecolor', 'none', 'FaceAlpha',.1);
end

plot(reshape(aperture, [200,ntrial], '.'));
%ylim([50 150])