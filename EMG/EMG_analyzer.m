
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code is specifically written to analyze EMG data from the
% the Balance & Fall FYP project
%      Rev July 2022 to include quiet standing trials (new project batch)
%      Rev Aug 2022 to retrieve latency and peak500, plot EMG data,
%      and make use of event markers to take relevant period
%      Rev to add function to combine left and right muscles data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;
fprintf('Hello. This is the MATLAB program to analyze EMG muscle activity!\n\r')

%% Before this, you first have to export QTM into MAT file
% (0) Some basic variables with pre-assigned value
amp = 1;   % Amplification factor <I decided to remove this to just use raw voltage>
emgactive = [];           % This is for subject-level emg-activity
showplot  = input('Do you want to see any graph? (0=No; 1=Yes)');


% (1) Point to project folder with UIgetdir function to open a folder dialog window
fprintf('First, please point to the directory of interest.\n');
while 1
    mypath  = uigetdir('/Users/sngqiwen/Desktop/M4/scholarly project/matlab/codes/data/');  
    pause(0.5)
    if (length(mypath)>1)
        fprintf('   %s\n',mypath);
        break
    else
        fprintf('!! No directory selected. Please select the directory!\n'); 
    end
end

% NOTE: Extract info if this is the Part-1 (quiet standing) or Part-2 (perturbation)
temppath = strsplit(mypath,'EMG');
temppath = temppath{1};
mypart   = str2num(temppath(length(temppath)-1));

% Check operating system of the current machine 
if ismac
    fprintf('   MacOS platform detected\n\r');   
    mypath = strcat(mypath, '/Raw data/');
elseif ispc
    fprintf('   Windows platform detected\n\r'); 
    mypath = strcat(mypath, '\Raw data\');
else
    fprintf('   Oops, platform not supported. Please Quit MATLAB.\n\r')
end


% (2) Then select and load all relevant files using UIgetfile function
fprintf('Now select and load ALL relevant files.\n');
while 1
    listfiles = uigetfile(mypath,'Select ALL input data file(s)','MultiSelect','on');
    if( ~iscell(listfiles) )   % To cater if you load 1 file only
        listfiles = mat2cell(listfiles,1);
        fprintf('   Total files selected: 1\n');  break
    else
        if(size(listfiles,2)>0)
            fprintf('   Total files selected: %d\n',size(listfiles,2)); 
            break
        else
            fprintf('!! Warning! Please select at least 1 file to proceed!\n'); 
        end
    end
    pause(1)
end

disp('   Ensure the folder and file list are correct. Press <Enter> to continue!'); pause


%% (3) Loop over each file in the list, then within each file we do preprocessing of each muscle
mysegment = cell({}); 
mysegment_out = cell({});
emgclean  = cell({}); 

for n = 1:size(listfiles,2)
    myfile = strcat(mypath,listfiles{n})       % Point to the file in the list one by one
    mymatfile = matfile(myfile);               % Load the MAT file
    tempfield = fieldnames(mymatfile);         % Obtain the fieldname of matfile
    mymatfile = mymatfile.(tempfield{2});      % Grab the content of the 2nd item containing EMG data

    try
        emg = mymatfile.Analog(2).Data;            % Extract all EMG data from the MAT file
        emg = amp*emg';                            % Amplify the data, if needed?            
        Nelectr = size(emg,2);                     % Define how many EMG electrodes
        Fs = mymatfile.Analog(2).Frequency;        % Sampling freq, change to (1) for S05                    
        sf = mymatfile.Analog(2).SamplingFactor;   % Sampling factor, change to (1) for S05
        myevents0 = mymatfile.Events;              % Note: data type is a struct

        % Important: I realized that there are event markers with different name!
        if (mypart==1)
            myevents  = get_subset(myevents0, 'Software trigger event');  % Part-1, quiet standing
        else  
            myevents  = get_subset(myevents0, 'Camera Sync Unit Event');  % Part-2, perturbation  
        end
        if(isempty(myevents))
            disp('No event detected. This means this is not a perturbation trial!');
        else
            disp('Multiple event detected. This is a perturbation trial!');
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
                        emg_preproc([Fs,sf,lofreq,hifreq], myevents, emg(:,mmm), mypart, showplot);
                    
            % Update variables where each row is each muscle and each column is data from each file
            mysegment{mmm,n} = mysegment_file;
            emgclean{mmm,n}  = emgclean_file;
            if(showplot), pause; end  % Pause for each muscle
        end
    catch
        warning('!! Error encountered while processing the EMG data!');
    end
    
    %% (4a) Remove bad trials. Important especially for perturbation trials
    % You write it inside a bracket the trial number of the set/block to be removed
    toremove = input('\nSelect the trial to remove (press <Enter> if none!) ');
    for mmm = 1:Nelectr
        temp = cell2mat(mysegment(mmm,:));  % Obtain the matrix of interest
        temp(:,toremove)=[];   % Then remove the columns of bad trials, if any 
        mysegment_out{mmm,n} = temp;
    end
    pause(1)
    fprintf('Finished processing the file-%d.......\n',n);
    pause(1)
    
end  %end of loop to run through each EMG file


%% (4b) Compute root mean square value as measure of EMG activity
% To compute EMG parameters, first reshape each muscle data to be a long vector, then compute 
% RMS value + variability
for mmm = 1:Nelectr
    temp = cell2mat(mysegment_out(mmm,:));
    mysegment_reshaped = reshape(temp, [size(temp,1)*size(temp,2) 1]);
    emgactive(mmm,1) = rms( mysegment_reshaped );
    emgactive(mmm,2) = std( mysegment_reshaped );
end


%% (5) Plot segments of interest of each muscle together with the average activity
fprintf('\nPreparing average muscle activity. Please wait! \n'); 
pause(0.5);
combine = input('Do you want to combine left and right side? (0=No; 1=Yes)');
xdata = [-0.25*Fs:1.0*Fs];      % Don't change this, because you have to change a few things!

if (combine == 1)   % If you choose to COMBINE left and right electrodes
    for mmm = 1:2:Nelectr
        figure(mmm)
        meanxxx = [];
        for f = 1:size(listfiles,2)
            % Each column is each file in the list, each row is each muscle. Use mysegment_out 
            % instead of mysegment!
            toplot = [mysegment_out(mmm,f); mysegment_out(mmm+1,f)];
            % Then, we combine left and right muscle data
            toplot = [cell2mat(toplot(1,1)) cell2mat(toplot(2,1))];
            % Plot can only display array, not cells
            plot( xdata, toplot, 'Color', [0.7 0.7 0.7] ); 
            hold on
            meanxxx = [meanxxx mean(toplot,2)];
        end
        xticklabels([-250:250:1000]);
        xline(0);    % Indicate trigger line!
        xlabel('Time (msec)');   ylabel('Amplitude (filtered)')
        figlabel = strsplit(mymatfile.Analog(2).Labels{mmm}, ' ');
        title(['Average EMG data (black) for ' string(figlabel(2)) ]);
        mean_to_plot(mmm,:) = mean(meanxxx,2);
        plot( xdata, mean_to_plot(mmm,:) ,'k');

        % Extract latency and peak of muscle activation for each trial
        latency(mmm) = get_latency( mean_to_plot(mmm,:), 0.25*Fs ,Fs, figlabel(2) );
        
        % Indicate latency with a line! Latency is meaningful only for perturbation trials
        if(~isempty(latency(mmm).Latency_sample)) && (mypart==2)
            xline(latency(mmm).Latency_sample,'r');    
        end
        ylim([0 1200]);
        hold off

        % Remove empty row. This is very useful function from the web........
        latency = latency( all(~cellfun(@isempty,struct2cell(latency))) );
    end

else   % If you choose to see the data SEPARATELY
    for mmm = 1:Nelectr
        figure(mmm)
        meanxxx = [];
        for f = 1:size(listfiles,2)
            % Each column is each file in the list, each row is each muscle. Use mysegment_out 
            % instead of mysegment!
            toplot = mysegment_out(mmm,f);
            % Plot can only display array, not cells
            plot( xdata, cell2mat(toplot), 'Color', [0.7 0.7 0.7] ); 
            hold on
            meanxxx = [meanxxx mean(cell2mat(toplot),2)];
        end
        xticklabels([-250:250:1000]);
        xline(0);    % Indicate trigger line!
        xlabel('Time (msec)');   ylabel('Amplitude (filtered)')
        title(['Average EMG data (black) for ' mymatfile.Analog(2).Labels{mmm}]);
        mean_to_plot(mmm,:) = mean(meanxxx,2);
        plot( xdata, mean_to_plot(mmm,:) ,'k');

        % Extract latency and peak of muscle activation for each trial
        latency(mmm) = get_latency(mean_to_plot(mmm,:), 0.25*Fs ,Fs, mymatfile.Analog(2).Labels{mmm});
        % Indicate latency with a line! Latency is meaningful only for perturbation trials
        if(~isempty(latency(mmm).Latency_sample)) && (mypart==2)
            xline(latency(mmm).Latency_sample,'r');    
        end
        ylim([0 1200]);
        hold off
    end
end

% Showing the average of the whole 16 EMG electrodes
figure(mmm+1)
    plot(mean_to_plot'); ylim([0 1000]);
    legend; % Show the legend so we know which muscle it is



%% (6) Lastly, save the processed file into a separate folder
pathout = strsplit(mypath,'Raw data');
outname = input('Give a name of the output file ("process_###"): ', 's');
fprintf('   Saving all necessary output files, please wait..... \n');

% Saving the figures....
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = num2str(get(FigHandle, 'Number'));
  set(0, 'CurrentFigure', FigHandle);
  savefig(fullfile(strcat(pathout{1},'/Processed data/'), [strcat(FigName,outname) '.fig']));
end

fprintf('   All figures of preprocessed data have been saved! \n');

if ismac
    outfile = strcat(pathout{1},'/Processed data/',outname);
else
    outfile = strcat(pathout{1},'\Processed data\',outname);
end

% Save the variables of interest
save( strcat(outfile,'.mat'),'Fs','sf',...
         'mysegment_out','emgactive',...
         'mean_to_plot','latency' );
  
fprintf('   Output file has been saved. Analysis completed! \n\r');


%------------------------------