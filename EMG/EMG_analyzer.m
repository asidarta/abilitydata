%clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code is specifically written to analyze EMG data from the
% the Fall FYP project (Last revision: 17 Mar 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Before this, you first have to export QTM into MAT file
% (1) Define main folder directory to suit different dataset accordingly;
mypath  = 'C:\Users\ananda.sidarta\Documents\MATLAB\'; % Main folder for your data
mysubj  = 'S03';              % Which subject number?
myblock = {'FF'}; % Four block conditions: F, FF, P, S
amp     = 500;                % Amplication factor to raw EMG data; Starting from S06, use lower number 0.5

%% (2) Loop over each experimental condition ('F' block,... etc), then within each condition we loop over
% each muscle (there are up to 14 electrodes). Within it, we again loop over several sub-blocks cos we 
% divided the trials for a better data saving process in QTM (Mocap).

for blk = 1:length(myblock) % condition blocks
    for mmm = 1:14  % electrode number
    mysegment = [];  % Declare a variable for EMG segments for all blocks!  
        for sub = 1:7   % sub-blocks...
        % Try-catch function will handle error situation where EMG data are missing or cannot be loaded. If error, warning msg will be called.
        try
            % Once filepath is defined, we will now load subject's MAT file
            myfile = strcat(mypath,mysubj,'\',mysubj,'_',myblock{blk},'_',num2str(sub),'.mat');
            mymatfile = matfile(myfile);              % Load the MAT file
            tempfield = fieldnames(mymatfile);        % Obtain the fieldname of matfile
            mymatfile = mymatfile.(tempfield{2});        % Grab the content of the 2nd item
            emg = mymatfile.Analog(2).Data(mmm,:);       % Extract EMG data from the MAT file; change to Analog (1) for S05
            emg = amp*emg';                              % Amplify the data, if needed?
            
            % Check if there are NaN undefined data, then throw them away!
            emg(isnan(emg)) = [];  
            % Take the first index of the analog data for S05
            Fs = mymatfile.Analog(2).Frequency;          % Sampling freq, change to (1) for S05                    
            sf = mymatfile.Analog(2).SamplingFactor;     % Sampling factor, change to (1) for S05
            myevents = mymatfile.Events;                 % Note: data type is a struct
            
            % Now, call the function to retrieve the segments we want!
            [mysegment, emgclean] = emg_preproc([Fs, sf], myevents, emg);
        catch
            warning('--- Error encountered while loading and processing the data!');
        end
        end  % finishing all sub-blocks...
    
        %% Save the mean and SD of each electrode
        all_mean{blk,mmm} = mean(mysegment,2); % Row (1) x column (2)
        all_std{blk,mmm}  = std(mysegment,0,2)/sqrt(size(mysegment,2));
        fprintf('Finished processing trials "%s", Electrode-%s\n', myblock{blk}, num2str(mmm));
        pause  % Pause allows me to check each EMG data one by one

    end  % finishing all muscle electrodes...
    %pause  % Pause allows me to check each condition one by one

end % finishing all trial conditions...

%% (3) NOTE: 'Save' function will save the variables of interest into a file
% save(strcat(mypath,mysubj),'all_mean','all_std');


%% (4) Now plot the curves. Here, each figure corresponds to each Muscle-pair! 
% Inside each figure, you will see 4 different conditions!
x_axis = -0.05*Fs: 2000/Fs: 0.5*Fs; 
for i = 0:6
    figure(i+1);
    subplot(2,2,1)
      plot(x_axis,all_mean{1,2*i+1}); hold on;
      plot(x_axis,all_mean{1,2*i+2}); 
      axis([-100,1000,-10,100]); xline(0);
      xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
      title('Average for "F"'); legend('Left','Right');
    subplot(2,2,2)
      plot(x_axis,all_mean{2,2*i+1}); hold on;
      plot(x_axis,all_mean{2,2*i+2});
      axis([-100,1000,-10,100]); xline(0);
      xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
      title('Average for "FF"'); legend('Left','Right');
    subplot(2,2,3)
      plot(x_axis,all_mean{3,2*i+1}); hold on;
      plot(x_axis,all_mean{3,2*i+2}); 
      axis([-100,1000,-10,100]); xline(0);
      xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
      title('Average for "P"'); legend('Left','Right');
    subplot(2,2,4)
      plot(x_axis,all_mean{4,2*i+1}); hold on;
      plot(x_axis,all_mean{4,2*i+2}); 
      axis([-100,1000,-10,100]); xline(0);
      xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
      title('Average for "S"'); legend('Left','Right');
end