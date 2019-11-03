
%%% Load normative dataset then clinical dataset

clear all; clc; close all
nframe = 100;

%mypath = 'C:\Users\ananda.sidarta\Documents\extracted angle\Normal gait\';
%myfile = 'new_SN006.txt';
%mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  
%var1 = reshape(mydata.hip_flexion_r, [nframe, size(mydata,1)/nframe]);

%mypath = 'C:\Users\ananda.sidarta\Documents\extracted angle\Stroke gait\';
%myfile = 'newS110.csv';
%mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0); 
%var2 = reshape(mydata.hip_flexion_r, [nframe, size(mydata,1)/nframe]);

load vector_coding_bilateral_hip_trajectories.mat
var1 = Right;

load vector_coding_bilateral_hip_trajectories_stroke.mat
var2 = Affected;

%diff = (var1(:,1:12)-var2);
%var1 = var1(11:100,:); var2 = var2(1:90,:);
subplot(1,2,1);
plot(var1,'g.'); hold on; plot(var2,'r');
ylabel('Joint angle (deg)'); xlabel('Normalized time');




% Using spm1D to understand temporal differences
%%%%%%%%%%%%%%%%%%%%%%%%%%
%var1 = zscore(var1); var2 = zscore(var2);
t  = spm1d.stats.ttest2(var2',var1');
ti = t.inference(0.05, true);
subplot(1,2,2);
ti.plot();
ti.plot_threshold_label();
ti.plot_p_values();


figure(3)
plot(mean(var1'),'b*'); 
hold on; 
plot(mean(var2'),'r*');