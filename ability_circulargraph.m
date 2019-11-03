%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a small piece of analysis using Connectogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear; clc; 

% Define path and filename.
mypath = 'D:\extracted angle\SN005\';
%myfile = 'SN005_0006_Towel_5Jul19'
myfile = 'towel_cut_L.txt';

ntrial = 7;  %%%%%<<<<<<<
ntime  = 400; %%%%%<<<<<<<

% Skips the first three rows of data as table.
% Note this is similar to a open table or dataframe in R.
mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  
% Remove the first and second column
mydata.trial = [];
mydata.phase=[];
colNames = mydata.Properties.VariableNames;

% Option 1: Here, I compute mean trajectories of all trials first,
% then estimate the correlation among these mean trajectories.
avg=[];
for i = 1:size(mydata,2)
    temp1  = mydata{:,i};
    temp2  = reshape(temp1, [ntime ntrial]);
    avg(:,i) = mean(temp2(1:200,:),2);
end

avg_table = array2table(avg);
avg_table.Properties.VariableNames = colNames;
% Here I'm confused whether I should keep the variable as table!!
x_corr = corr( table2array(avg_table) );
x_corr(find(abs(x_corr) < 0.95)) = 0;
%x_corr = x_corr + 1;

% Now plot the circular graph!
%figure(1)
%circularGraph(abs(x_corr), 'Label', colNames);



% Option 2: Why don't we take the whole full trajectories and estimate the
% correlation, rather than looking into per trial basis.
avg = table2array(mydata);
% Compute the mean of correlation
corr_mean = corr(avg);
corr_mean(find(abs(corr_mean) < 0.85)) = 0;
% Now plot the circular graph!
figure(2)
circularGraph(abs(corr_mean), 'Label', colNames);




% Option 3: But we can also compute the correlations between trajectories
% of each trial then compute the mean correlation values.
avg=[];
enddata = size(mydata,1)/ntrial;
for j = 1:ntrial
    temp3 = table2array(mydata);
    temp4 = corr( temp3( (j-1)*enddata+201:j*enddata,:) );
    avg(:,:,j) = temp4;
end

% Compute the mean of correlation
corr_mean = mean(avg,3);
corr_mean(find(abs(corr_mean) < 0.8)) = 0;
% Now plot the circular graph!
%figure(3)
%circularGraph(abs(corr_mean), 'Label', colNames);




% Option 4: Using Dynamic Time Warping (DTW)
avg = table2array(mydata);
similar = [];
for i = 1:size(avg,2)
    for j = 1:size(avg,2)
        similar(i,j) = dtw(avg(1:200,i),avg(1:200,j));
    end
end
corr_mean = similar;

% Now plot the circular graph!
corr_mean(find(abs(corr_mean) < 3e+4)) = 0;
figure(4)
circularGraph(abs(corr_mean), 'Label', colNames);