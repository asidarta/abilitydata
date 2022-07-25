clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code is specifically written to analyze MVC data from the
% the Fall FYP project (Last revision: 17 Mar 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% (1) Load the textfile
fpath  = 'C:\Users\ananda.sidarta\Documents\MATLAB\'; % Main folder for the data
mysubj  = 'S04';  % Subject number
myfile = strcat(fpath,'\',mysubj,'_','MVC','.txt')

% Wee Kiat's method is better because we are loading a dataset with column header (string)
% By treating the data into a table, we can manipulate easily directly
% Also, column headers will be readily available so you won't get confused
mvc = readtable(myfile, 'Delimiter', '\t');
mvcfile_name = mvc.Properties.VariableNames;
mvc = mvc{:, 2:end};    % Convert table into 2D array (matrix) format
mvc = 0.05*mvc;    % Amplify the data 1000x

Fs  = 2000;  sf  = 10;

%% (2) Loop over each subject
% Loop over all electrodes, 12 for first two subjects
    for elect = 1:12
        each_mvc = mvc(:,elect);
        mymvc = [];  % the variable to contain MVC segment of interest
    
        % Try-catch function will handle error situation where EMG data are missing or cannot be loaded. If error, warning msg will be called.
        try
            % Check if there are NaN undefined data, then throw them away!
            each_mvc(isnan(each_mvc)) = [];  
            % Tilda can be used in the output part when calling the function
            % but cannot use it to define the inputs of the function.
            [mymvc] = emg_preproc([Fs, sf], [], each_mvc);

        catch
            warning('--- Error encountered while processing the data!');
        end
        
        pause
        
    % For MVC, compute the mean excitation of particular muscle. Some papers computed 
    % root-mean square instead of mean. That doesn't differ much!
    all_mean{elect} = mean(mean(mymvc,2));    % Each electrode will give 1 number
    fprintf('Finished processing subject "%s", Electrode-%s\n', num2str(mysubj), num2str(elect));
    
    end  % finishing all muscle electrodes...

    
%% (3) NOTE: 'Save' function will save the variables of interest into a file
% save(strcat(fpath,mysubj),'all_mean','all_std');