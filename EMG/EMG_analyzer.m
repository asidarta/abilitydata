
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code is specifically written to analyze EMG data from the
% the Balance & Fall FYP project (Revised: July 2022)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
fprintf('Hello. This is the MATLAB program to analyze EMG muscle activity!\n\r')

%% Before this, you first have to export QTM into MAT file
% (0) Some basic variables with pre-assigned value
amp     = 500;            % Default:
Nelectr = 16;             % There are 16 EMG electrodes
emgactive = [];           % This is for subject-level emg-activity
showplot  = input('Do you want to see any graph? (0=No; 1=Yes)');
% (For Enoch: amplication factor to raw EMG data; Starting from S06, use lower number 0.5)
% (For Enoch: there were a total of 14 electrodes)


% (1) Point to project folder with UIgetdir function to open a folder dialog window
fprintf('First, please point to the directory of interest.\n');
while 1
    mypath  = uigetdir('C:\Users\ananda.sidarta\Desktop\MVC');  
    pause(0.5)
    if (length(mypath)>1)
        mypath 
        break
    else
        fprintf('!! No directory selected. Please select the directory!\n'); 
    end
end

% (2) Then select and load all relevant files using UIgetfile function
fprintf('Now select and load ALL of the necessary .MAT files.\n');
while 1
    mypath  = strcat(mypath, '\Raw data\');
    listfiles = uigetfile(mypath,'Select ALL input data file(s)','MultiSelect','on');
    if(size(listfiles,2)>0)
        fprintf('Total files selected: %d\n\r',size(listfiles,2)); 
        break
    else
        fprintf('!! Warning! Please select at least 2 file to proceed!\n\r'); 
    end
    pause(1)
end

disp('Ensure the folder and file list are correct. Press <ENTER> to continue!'); pause


%% (3) Loop over each file in the list, then within each file we do preprocessing of each muscle
mysegment = cell({}); 
emgclean  = cell({}); 

for n = 1:size(listfiles,2)
    myfile = strcat(mypath,listfiles{n})       % Point to the file in the list one by one
    mymatfile = matfile(myfile);               % Load the MAT file
    tempfield = fieldnames(mymatfile);         % Obtain the fieldname of matfile
    mymatfile = mymatfile.(tempfield{2});      % Grab the content of the 2nd item

    try
        emg = mymatfile.Analog(2).Data;            % Extract all EMG data from the MAT file 
        emg = amp*emg';                            % Amplify the data, if needed?            
        Fs = mymatfile.Analog(2).Frequency;        % Sampling freq, change to (1) for S05                    
        sf = mymatfile.Analog(2).SamplingFactor;   % Sampling factor, change to (1) for S05
        myevents0 = mymatfile.Events;              % Note: data type is a struct
        % NOTE: Take the first index of the analog data for S05 for Enoch's subject
        % NOTE: change to Analog (1) for S05 for Enoch's subject

        % I realized that there are event markers which carry different name...
        myevents = get_subset(myevents0, 'Camera Sync Unit Event');

        if(isempty(myevents)||length(myevents)<=3)
            disp('No event detected. This means this is not a perturbation trial!');
            mvcdata = true;
        else
            disp('Multiple event detected. This is a perturbation trial!');
            mvcdata = false;
        end
    catch
        warning('!! Error encountered while loading the EMG file!');
    end

    try           
        for mmm = 1:Nelectr     % for each electrode number
            % This is to bpf filter EMG data contaminated by ventricular artifacts!!
            if (mmm==7)||(mmm==8)||(mmm>=13),  lofreq = 60; hifreq = 400;
            else,                              lofreq = 10; hifreq = 400;
            end
            fprintf('Processing EMG data for: %s \n', mymatfile.Analog(2).Labels{mmm});
            [mysegment_file, emgclean_file] = ...
                        emg_preproc([Fs,sf,lofreq,hifreq], myevents, emg(:,mmm), mvcdata, showplot);
                    
            % Update variables where each row is each muscle and each column is data from each file
            mysegment{mmm,n} = mysegment_file;
            emgclean{mmm,n}  = emgclean_file;
            if(showplot), pause; end
        end
    catch
        warning('!! Error encountered while processing the EMG data!');
    end
    
    fprintf('Finished processing the file-%d.......\n',n);
    pause(1.5)
end  %end of for-loop


%% (4) Compute root mean square value as measure of EMG activity
% First reshape each muscle data to be a long vector, then compute RMS value
for mmm = 1:Nelectr
    temp = cell2mat(mysegment(mmm,:));
    mysegment_reshaped = reshape(temp, [size(temp,1)*size(temp,2) 1]);
    emgactive(mmm,1) = rms( mysegment_reshaped );
    emgactive(mmm,2) = std( mysegment_reshaped );
end
        


%% (5) Plot segments of interest of each muscle together with the average activity
fprintf('\n\rPreparing average muscle activity. Please wait........\n\r');

if(isempty(myevents)||length(myevents)<=3)
    xdata = 1:2000;
else
    xdata = [-0.25*Fs:1.0*Fs];
end

for mmm = 1:Nelectr
    figure(mmm)
    meanxxx = [];
    for f = 1:size(listfiles,2)
        % each column is each file in the list, each row is each muscle
        toplot = mysegment(mmm,f);
        % Plot can only display array, not cells
        plot( xdata, cell2mat(toplot), 'Color', [0.8 0.8 0.8] ); 
        hold on
        if(length(myevents)>1)
            %xticklabels([-250:250:1000]);
            xline(0);    % Indicate trigger line!
            xlabel('Time (msec)');   ylabel('Amplitude (filtered)')
        else
            xlabel('Samples');   ylabel('Amplitude (filtered)')
        end
        title(['Average EMG data (black) for ' mymatfile.Analog(2).Labels{mmm}]);
        meanxxx = [meanxxx mean(cell2mat(toplot),2)];
    end
    mean_to_plot = mean(meanxxx,2);
    plot( xdata, mean_to_plot ,'k');
    latency(mmm) = get_latency(mean_to_plot, 0.25*Fs ,Fs);  % Extract latency of muscle activation
    xline(latency(mmm)*Fs,'r');    % Indicate latency with a line!
    hold off
end


%% (6) Lastly, save the processed file
pathout = strsplit(mypath,'Raw data');
outname = input('Give a name of the output file ("process_###"): ', 's');
save( strcat(pathout{1},'\Processed data\',outname,'.mat'), ...    % output path
      'mysegment','emgactive','latency');               % variables of interest
fprintf('Output file has been saved\n\r');
  
  
  

%%
% %% (2) Loop over each experimental condition ('F' block,... etc), then within each condition we loop over
% % each muscle (there are up to 14 electrodes). Within it, we again loop over several sub-blocks cos we 
% % divided the trials for a better data saving process in QTM (Mocap).
% 
% for blk = 1:length(myblock)    % condition blocks
%     for mmm = 1:Nelectrode     % electrode number
%     mysegment = [];    % Declare a variable for EMG segments for all blocks!  
%         for sub = 1:7    % sub-blocks...
%         % Try-catch function will handle error situation where EMG data are missing or cannot be loaded. If error, warning msg will be called.
%         try
%             % Once filepath is defined, we will now load subject's MAT file
%             myfile = strcat(mypath,'\',mysubj,'\',mysubj,'_',myblock{blk},'.mat');
%             mymatfile = matfile(myfile);               % Load the MAT file
%             tempfield = fieldnames(mymatfile);         % Obtain the fieldname of matfile
%             mymatfile = mymatfile.(tempfield{2});      % Grab the content of the 2nd item
%             emg = mymatfile.Analog(2).Data(mmm,:);     % Extract EMG data from the MAT file; change to Analog (1) for S05
%             emg = amp*emg';                            % Amplify the data, if needed?
%             
%             % Check if there are NaN undefined data, then throw them away!
%             emg(isnan(emg)) = [];  
%             % Take the first index of the analog data for S05
%             Fs = mymatfile.Analog(2).Frequency;          % Sampling freq, change to (1) for S05                    
%             sf = mymatfile.Analog(2).SamplingFactor;     % Sampling factor, change to (1) for S05
%             myevents = mymatfile.Events;                 % Note: data type is a struct
%             
%             % Now, call the function to retrieve the segments we want!
%             [mysegment, emgclean] = emg_preproc([Fs, sf], myevents, emg);
%         catch
%             warning('--- Error encountered while loading and processing the data!');
%         end
%         end  % finishing all sub-blocks...
%     
%         %% Save the mean and SD of each electrode
%         all_mean{blk,mmm} = mean(mysegment,2); % Row (1) x column (2)
%         all_std{blk,mmm}  = std(mysegment,0,2)/sqrt(size(mysegment,2));
%         fprintf('Finished processing trials "%s", Electrode-%s\n', myblock{blk}, num2str(mmm));
%         pause  % Pause allows me to check each EMG data one by one
% 
%     end  % finishing all muscle electrodes...
%     %pause  % Pause allows me to check each condition one by one
% 
% end % finishing all trial conditions...
% 
% %% (3) NOTE: 'Save' function will save the variables of interest into a file
% % save(strcat(mypath,mysubj),'all_mean','all_std');
% 
% 
% %% (4) Now plot the curves. Here, each figure corresponds to each Muscle-pair! 
% % Inside each figure, you will see 4 different conditions!
% x_axis = -0.05*Fs: 2000/Fs: 0.5*Fs; 
% for i = 0:6
%     figure(i+1);
%     subplot(2,2,1)
%       plot(x_axis,all_mean{1,2*i+1}); hold on;
%       plot(x_axis,all_mean{1,2*i+2}); 
%       axis([-100,1000,-10,100]); xline(0);
%       xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
%       title('Average for "F"'); legend('Left','Right');
%     subplot(2,2,2)
%       plot(x_axis,all_mean{2,2*i+1}); hold on;
%       plot(x_axis,all_mean{2,2*i+2});
%       axis([-100,1000,-10,100]); xline(0);
%       xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
%       title('Average for "FF"'); legend('Left','Right');
%     subplot(2,2,3)
%       plot(x_axis,all_mean{3,2*i+1}); hold on;
%       plot(x_axis,all_mean{3,2*i+2}); 
%       axis([-100,1000,-10,100]); xline(0);
%       xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
%       title('Average for "P"'); legend('Left','Right');
%     subplot(2,2,4)
%       plot(x_axis,all_mean{4,2*i+1}); hold on;
%       plot(x_axis,all_mean{4,2*i+2}); 
%       axis([-100,1000,-10,100]); xline(0);
%       xlabel('Time (msec)'); ylabel('Amplitude (filtered)')
%       title('Average for "S"'); legend('Left','Right');
% end