%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code compares kinematic trajectories of LEFT and RIGHT limbs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all; 
mean_array=[]; std_array=[];

% Define path and filename for STROKE.
mypath = 'D:\extracted angle\Stroke gait\';
myfile = 'newS142.csv';
%98,108,110,128,129,136,137,142
% R,  L,  R,  R,  R,  R,  R,  L

% Define path and filename for normal gait.
%mypath = 'C:\Users\ananda.sidarta\Documents\extracted angle\Normal gait\';
%myfile = 'new_SN035.txt';

% Skips the first three rows of data as table. Note this is similar to a open 
% table or dataframe in R.
mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  
colNames = mydata.Properties.VariableNames;
nframe = 100;
ntrial = size(mydata,1)/nframe;  %<<<

%%%%% This is the place where you have to define variables of interest!!!
var1 = mydata.hip_flexion_r;  
var2 = mydata.hip_flexion_l;  

%plot(var1); hold on; plot(var2);

% Can we produce Angle-angle Plot of two joints???
hold on; figure(1)
plot(var1, var2, 'b.');
axis square;
axis([-70,70,-70,70]);
%xlabel('Left Hip (flex, deg)'); ylabel('Right Hip (flex, deg)');
%idx = DBSCAN([var1 var2],2,1); gscatter(var1,var2,idx);
hold off;


if(1)
%-------------------------------------------------------------------------
% Vector coding method: We extract phase information from the Angle-angle
% plot in the form of an angle; using atan between 2 consecutive points.
%-------------------------------------------------------------------------
% Shall we modify the data to be per trial? For stroke subjects the data
% is not normalized, so we plot everything.
var1 = reshape(var1, [nframe,ntrial]);
var2 = reshape(var2, [nframe,ntrial]);
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

%plot(mean_ang2 + 180, 'r.-'); 
%hold on;

ang2_circ = ang2* pi / 180;
mean_circ = circ_mean(ang2_circ');
mean_circ = mean_circ * 180 / pi;   % convert rad to degree
std_circ  = circ_std(ang2_circ');
std_circ  = std_circ * 180 / pi;


figure(2) % PLOT the mean and std using circular mean function
plot(mean_circ + 180, 'b-'); hold on; 
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
in_phase   = sum(group == 'r')/mysize(1)%*ntrial
knee_phase = sum(group == 'g')/mysize(1)%*ntrial
anti_phase = sum(group == 'b')/mysize(1)%*ntrial
hip_phase  = sum(group == 'w')/mysize(1)%*ntrial

end


% IMPORTANT: Let's create an array of mean and variability for different subjects,
% using the circular statistics respectively.
mean_array = [mean_array; mean_circ];
std_array  = [std_array;  std_circ ];

normal_gait_std = circ_rad2ang(circ_std(circ_ang2rad(mean_array)));
normal_gait     = circ_rad2ang(circ_mean(circ_ang2rad(mean_array)));



% This is to compute average phase components (in-phase, etc) per subject!!
if(1) 
in_phase=[]; anti_phase=[]; hip_phase=[]; knee_phase=[];
for i = 1:size(mean_array,1)
    temp = mean_array(i,:); i
    group = ordinal(temp+180,{'w','r','g','b','w','r','g','b','w'},...
    [],[0,22.5,67.5,112.5,157.5,202.5,247.5,292.5,337.5,360]);
    % Compute the number of elements within each ordinal range!
    in_phase   = [in_phase;   sum(group == 'r')/size(mean_array,2)]; 
    knee_phase = [knee_phase; sum(group == 'g')/size(mean_array,2)]; %y
    anti_phase = [anti_phase; sum(group == 'b')/size(mean_array,2)];
    hip_phase  = [hip_phase;  sum(group == 'w')/size(mean_array,2)]; %x
end

% Now, plot bar-graph with error bar to compare different phases averaged
% across n=8 subjects!
phases = [mean(in_phase);mean(anti_phase);mean(hip_phase);mean(knee_phase)];
stderr = [std(in_phase);std(anti_phase);std(hip_phase);std(knee_phase)]/sqrt(8);

figure(3)
bar(phases); hold on; 
ylim([0,0.9]); ylabel('Proportion');
er = errorbar(1:4,phases,stderr,stderr); 
er.Color = [0 0 0];                 
er.LineStyle = 'none';
xticklabels({'in-phase','anti-phase','RIGHT-lead','LEFT-lead'})
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

