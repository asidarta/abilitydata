%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code compares kinematic trajectories of LEFT and RIGHT limbs
% using VECTOR CODING technique for Folding Towel........
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; 
mean_array=[]; std_array=[];

% Shall we perform while loop for all subjects?
subj = {'10'};
%subj = {'05','08','11','13'};

% Define filename for STROKE. Which task you want to see?
myfile = 'towel_cut_L.txt';


% READY to loop through all subjects!!!
for k = 1 : length(subj)

mypath = 'D:\extracted angle\SN0';
mypath = strcat(mypath, subj{k},'\');  k

% Skips the first three rows of data as table. Note this is similar to a open 
% table or dataframe in R.
mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  

% Data preparation: Loop across different joint angles. Our joint angles
% have been normalized to contain 200 x 2 phases = 400 time points.
nframe = 400;
ntrial = size(mydata,1)/nframe;  %<<<

% The column naming follows the existing textfile I received...
colNames = mydata.Properties.VariableNames;


%%%% This is the place where you have to define variables of interest! %%%%
var1 = mydata.RElbowX;  
var2 = mydata.LElbowX;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Can we produce Angle-angle Plot of two joints???
figure(2); hold on; 
pointsize = 5;
[unique_groups, ~, group_idx] = unique(mydata.phase);
num_groups = size(unique_groups, 1);
scatter(var1, var2, pointsize, mydata.phase);
axis square;
axis([-50,100,-50,100]);
cmap = jet(num_groups);    % or build a custom color map
colormap( cmap );          % apply the custom map
hold off;


%-------------------------------------------------------------------------
% Vector coding method: We extract phase information from the Angle-angle
% plot in the form of an angle; using atan between 2 consecutive points.
%-------------------------------------------------------------------------
% Shall we modify the data to be per trial? For stroke subjects the data
% is not normalized, so we plot everything.
var1 = reshape(var1, [nframe,ntrial]);
var2 = reshape(var2, [nframe,ntrial]);


% Compute the average trajectory for all trials....
movegui(figure(1),'east')
mean_r = mean(var1,2);   se_r = std(var1,1,2)/sqrt(ntrial);
mean_l = mean(var2,2);   se_l = std(var2,1,2)/sqrt(ntrial);
% Now plor both average curves (L/R) in one figure
%errorbar(mean_r,se_r,'bx'); hold on; errorbar(mean_l,se_l,'rx');
% Add standard error shading using patch function...
lo = mean_r - se_r; hi = mean_r + se_r;
nsample = length(lo);
mycolor = 'r';   % red for Right -------------
hl = line(1:length(lo),mean_r,'color',mycolor);
hold on;
hp = patch([(1:nsample)'; (nsample:-1:1)'; 1], ...
           [lo; hi(nsample:-1:1); lo(1)], mycolor);
set(hp, 'facecolor', mycolor, 'edgecolor', 'none', 'FaceAlpha',.1);
lo = mean_l - se_l;
hi = mean_l + se_l;
nsample = length(lo);
mycolor = 'b';   % blue for Left -------------
hl = line(1:length(lo),mean_l,'color', mycolor);
hp = patch([(1:nsample)'; (nsample:-1:1)'; 1], ...
           [lo; hi(nsample:-1:1); lo(1)], mycolor);
set(hp, 'facecolor', mycolor, 'edgecolor', 'none', 'FaceAlpha',.1);
title('Mean angular trajectories (degree)');
xlabel('Normalized time-points'); ylabel('Average angle');
line([200 200],[min(mean_l)-10 max(mean_l)+10], 'color', 'black'); 





%%%%%%%%%% Do you want to take the whole time frames? %%%%%%%%%%%%%%%%%%%
var1 = var1(1:200,:); var2=var2(1:200,:);    %%<<<<<<<<<
aaa  = [mean(var1,2) mean(var2,2)];

mysize = size(var1);

for col = 1:mysize(2)
    for row = 1:mysize(1)-1
        x(row,col) = var1(row+1,col)-var1(row,col);  % joint angle1
        y(row,col) = var2(row+1,col)-var2(row,col);  % joint angle2
        ang(row,col) = atan2(y(row,col),x(row,col));  % This is in radian
        ang(row,col) = ang(row,col) * 180/pi;        % convert to degree
        if (ang(row,col)<0 && ang(row,col)>-180)
            ang2(row,col) = ang(row,col) + 360;
        else
            ang2(row,col) = ang(row,col);
        end
    end
end

% After computing the vectors, let's plot the phase angles!
%figure(1);
%plot(ang2(:,3:4),'r-'); hold on;
%plot(ang(:,1),'b-');

% Compute circular mean (circular statistics) of the phase (Hamill, 2008).
% First find the mean cos and sin values for each time point. This operation 
% goes on for all rows!
y_ang = mean( sind(ang2),2 );
x_ang = mean( cosd(ang2),2 );

% Then find the mean, then obtain back the angle using trigonometry.
mean_ang = atan2(y_ang, x_ang);
mean_ang = mean_ang * 180/pi;   % convert rad to degree
if (mean_ang<0 & mean_ang>-180)
    mean_ang2 = mean_ang + 360;
else
    mean_ang2 = mean_ang;
end

%plot(mean_ang2 + 180, 'r.-'); %hold on;

ang2_circ = ang2* pi / 180;
mean_circ = circ_mean(ang2_circ');
mean_circ = mean_circ * 180 / pi;   % convert rad to degree
std_circ  = circ_std(ang2_circ');
std_circ  = std_circ * 180 / pi;

figure(3) % PLOT the mean and std using circular mean function
  plot(mean_circ + 180, 'b.'); hold on; 
  xlim([0,mysize(1)]); ylim([0 360])
  plot(std_circ);
  xlabel('Normalized time points'); ylabel('Phase angle (degree)');

% Create patches according to Chang et al, 2008, dividing areas where the
% angles are in-phase, out-of-phase, etc.
x  = ([0,mysize(1)]);
y1 = [ 22.5,  22.5];
y2 = [ 67.5,  67.5];
y3 = [112.5, 112.5];
y4 = [157.5, 157.5];
y5 = [202.5, 202.5];
y6 = [247.5, 247.5];
y7 = [292.5, 292.5];
y8 = [337.5, 337.5];
I  = patch([x fliplr(x)],[y1 fliplr(y2)], 'r'); alpha(0.1);  hold on
II = patch([x fliplr(x)],[y2 fliplr(y3)], 'g'); alpha(0.1);
III= patch([x fliplr(x)],[y3 fliplr(y4)], 'b'); alpha(0.1);
IV = patch([x fliplr(x)],[y4 fliplr(y5)], 'w'); alpha(0.1);
V  = patch([x fliplr(x)],[y5 fliplr(y6)], 'r'); alpha(0.1);
VI = patch([x fliplr(x)],[y6 fliplr(y7)], 'g'); alpha(0.1);
VII= patch([x fliplr(x)],[y7 fliplr(y8)], 'b'); alpha(0.1);
xlim([0,nframe]);

hold off
group = ordinal(mean_circ+180,{'w','r','g','b','w','r','g','b','w'},...
    [],[0,22.5,67.5,112.5,157.5,202.5,247.5,292.5,337.5,360]);

% Compute the number of elements of a specific value!
in_phase   = sum(group == 'r')/mysize(1)
knee_phase = sum(group == 'g')/mysize(1)
anti_phase = sum(group == 'b')/mysize(1)
hip_phase  = sum(group == 'w')/mysize(1)



% IMPORTANT: Let's create an array of mean and variability for different subjects,
% using the circular statistics respectively.
mean_array = [mean_array; mean_circ];
std_array  = [std_array;  std_circ ];

normal_vect_std = circ_rad2ang(circ_std(circ_ang2rad(mean_array)));
normal_vect     = circ_rad2ang(circ_mean(circ_ang2rad(mean_array)));

figure(3) % PLOT the mean and std using circular mean function
  plot(normal_vect + 180, 'b.'); hold on; 
  xlim([0,mysize(1)]); ylim([0 360])
  plot(normal_vect_std);
  xlabel('Normalized time points'); ylabel('Phase angle (degree)');


% This is to compute average phase components (in-phase, etc) for all subjects.
in_phase=[]; anti_phase=[]; left_phase=[]; right_phase=[];
for i = 1:size(mean_array,1)
    temp = mean_array(i,:); %i
    group = ordinal(temp+180,{'w','r','g','b','w','r','g','b','w'},...
    [],[0,22.5,67.5,112.5,157.5,202.5,247.5,292.5,337.5,360]);
    % Compute the number of elements within each ordinal range!
    in_phase    = [in_phase;    sum(group == 'r')/size(mean_array,2)]; 
    right_phase = [right_phase; sum(group == 'g')/size(mean_array,2)]; %y
    anti_phase  = [anti_phase;  sum(group == 'b')/size(mean_array,2)];
    left_phase  = [left_phase;  sum(group == 'w')/size(mean_array,2)]; %x
end

% Now, plot bar-graph with error bar to compare different phases averaged
% across n=8 subjects!
phases = [mean(in_phase);mean(anti_phase);mean(right_phase);mean(left_phase)];
stderr = [std(in_phase);std(anti_phase);std(right_phase);std(left_phase)]/sqrt(8);

if k == length(subj)
    figure(4)
    bar(phases); hold on; 
    ylim([0,0.9]); ylabel('Proportion');
    er = errorbar(1:4,phases,stderr,stderr); 
    er.Color = [0 0 0];                 
    er.LineStyle = 'none';
    xticklabels({'in-phase','anti-phase','LEFT-lead','RIGHT-lead'})
    view([90 -90]);   % Nice function to rotate plot!!!
    hold off;
end


if(0) % This is to plot average plots
    affected = mean_array([1,10,3,4,5,6,7,16],:);
    normal   = mean_array([9,2,11,12,13,14,15,8],:);
    affected_std  = circ_rad2ang(circ_std(circ_ang2rad(affected)));
    affected_mean = circ_rad2ang(circ_mean(circ_ang2rad(affected)));
    normal_std  = circ_rad2ang(circ_std(circ_ang2rad(normal)));
    normal_mean = circ_rad2ang(circ_mean(circ_ang2rad(normal)));
    
    figure(4) % PLOT the mean and std using circular mean function
    plot(normal_mean + 180, 'b.-'); hold on; 
    xlim([0,mysize(1)]); ylim([0 360])
    plot(normal_std);
    xlabel('Normalized time points'); ylabel('Phase angle (degree)'); 
    plot(affected_mean + 180, 'r.-');
    plot(affected_std);
    
    knee_normal   = mean_array([9,2,11,12,13,14,15,8],:); 
    
    affected_phases = [mean(in_phase([1,10,3,4,5,6,7,16],:));
                       mean(anti_phase([1,10,3,4,5,6,7,16],:));
                       mean(hip_phase([1,10,3,4,5,6,7,16],:));
                       mean(knee_phase([1,10,3,4,5,6,7,16],:))];
                   
    affected_stderr = [std(in_phase([1,10,3,4,5,6,7,16],:));
                       std(anti_phase([1,10,3,4,5,6,7,16],:));
                       std(hip_phase([1,10,3,4,5,6,7,16],:));
                       std(knee_phase([1,10,3,4,5,6,7,16],:))]/sqrt(8); 
                   
    normal_phases   = [mean(in_phase([9,2,11,12,13,14,15,8],:));
                       mean(anti_phase([9,2,11,12,13,14,15,8],:));
                       mean(hip_phase([9,2,11,12,13,14,15,8],:));
                       mean(knee_phase([9,2,11,12,13,14,15,8],:))];
                   
    normal_stderr   = [std(in_phase([9,2,11,12,13,14,15,8],:));
                       std(anti_phase([9,2,11,12,13,14,15,8],:));
                       std(hip_phase([9,2,11,12,13,14,15,8],:));
                       std(knee_phase([9,2,11,12,13,14,15,8],:))]/sqrt(8);    

    subplot(1,2,1)
    bar(normal_phases); hold on;
    ylim([0,0.9]); ylabel('Proportion');
    er = errorbar(1:4,normal_phases,normal_stderr,normal_stderr); 
    er.Color = [0 0 0];                 
    er.LineStyle = 'none';
    xticklabels({'in-phase','anti-phase','RIGHT-lead','LEFT-lead'})
    view([90 -90]);   % Nice function to rotate plot!!!
    subplot(1,2,2)
    bar(affected_phases); hold on;
    ylim([0,0.9]); ylabel('Proportion');
    er = errorbar(1:4,affected_phases,affected_stderr,affected_stderr); 
    er.Color = [0 0 0];                 
    er.LineStyle = 'none';
    view([90 -90]);   % Nice function to rotate plot!!!
    hold off; 
end

pause  % Pause for each subject and press any key to continue!

end