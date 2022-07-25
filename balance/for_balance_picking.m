
clear; clc; close all


%%
% Define the working directory. 
wd = 'C:\Users\ananda.sidarta\Documents\MATLAB\balance\cop_balance\';
cd(wd);

COP_each_leg = [];   % COP metrics for each leg separately
COP_contrast = [];   % COP metrics to contrast L/R side
meanSpect = [];
showPlot  = 1;       % Want to show the figure?

age = [3,6,2,2,4,3,5,2,5,2,6,6,5,2,3,6,5,5,6,3,3,4,5,3,4,4,4,5,2,6,5,5,3,2,3,6,4,3,5,5,3,5,5,7,5,6,5,6,...
       5,2,5,4,6,7,6,5,5,6,7,6,4,6,3,5,7,4,5,2,7,5,7,6,6,6,6,5,6,7,7,5,5,5,3,2,6,5,5,6,5,5,6,6,6,6,6,6,7,5,2,6];

   
   
for subj = 1:1      % Looping for EVERY subject (numbered 1 - 100)
for task = 1:1      % Looping for EVERY pick-up task (left=1, right=2)

    if task == 1,   taskStr = 'left';   % String handler for task
    else,  taskStr = 'right';    
    end
    
    % Now open the relevant file based on subject and task type
    myfile = strcat('balance_segmented_cofp_',taskStr,'_',num2str(subj) );
    
    % Then load the textfile as Table. Please skip the first FIVE rows!
    T = readtable(myfile,'HeaderLine',5,'ReadVariableNames',0);%,'Format','%f%f%f%f%f%f%f');
    T(:,1)=[];              % Delete the first column!
    d = table2array(T);     % Convert to numeric ARRAY
    d = 1000*d;             % Convert to millimeter(mm) unit.
    Fs = 2000;              % Sampling frequency of the recording
    
    % Design a 4th order Butterworth Low Pass Filter, cutoff freq = 5Hz
    cutfreq = 5/(Fs/2); 
    [B,A] = butter(4,cutfreq);
    
    for trial = 0:2      % Looping for EVERY repetition, rep 1 to 3
        %% Deal with the LEFT side data first ---------------------------------------------------------
        side = 1;
        dd = d(:,1+trial*3:2+trial*3);     % automatically take the correct columns for left FP!!
        dd_clean_L = filtfilt(B, A, dd);   % filter the input data, subscript L for left!
    
        % Demean the cleaned data
        demeanedX_L = dd_clean_L(:,1)-mean(dd_clean_L(:,1));
        demeanedY_L = dd_clean_L(:,2)-mean(dd_clean_L(:,2));
        % Flip the axis??
        demeanedX_L = -demeanedX_L; demeanedY_L = -demeanedY_L;
    
        % Compute time domain parameters for x-axis (LR) and y-axis (AP)
        [pk2pkx,meanDx,rmsDx,meanVelx,mfreqx,slopex] = time_domain(demeanedX_L, Fs);
        [pk2pky,meanDy,rmsDy,meanVely,mfreqy,slopey] = time_domain(demeanedY_L, Fs);
    
        % Compute 95% conf ellipse of COP
        [ellip,x0,y0,area] = error_ellipse([demeanedX_L demeanedY_L]);
    
        % Compute frequency domain parameters for x-axis (LR) and y-axis (AP)
        [powerX,p50X,p95X,spm1X,spm2X,cfreqX,freqdX,~,~] = freq_domain(demeanedX_L, Fs);
        [powerY,p50Y,p95Y,spm1Y,spm2Y,cfreqY,freqdY,f,Y] = freq_domain(demeanedY_L, Fs);
 
        % Compile ALL variables for quantifying COP
        COP_each_leg = [COP_each_leg; [subj,task,side,age(subj),meanDx,meanDy,rmsDx,rmsDy,...
                               pk2pkx,pk2pky,meanVelx,meanVely,slopex,slopey,mfreqx,mfreqy,area, ...
                               powerX,powerY,p50X,p50Y,p95X,p95Y,cfreqX,cfreqY,freqdX,freqdY ] ];
    
        % Visualize the data! Separate reaching downward and return upward portion
        if(showPlot) 
        fff = figure('WindowState','maximized');
        fff.Name = strcat('COP_trajectory_Eyes_',taskStr);
        subplot(2,4,2)
            plot(demeanedX_L(1:3600),   demeanedY_L(1:3600),   'b.'); hold on;  
            plot(demeanedX_L(3601:7200),demeanedY_L(3601:7200),'r.');
            title('LEFT FOOT');
            xlabel('Left-Right'); ylabel('Anterior-Posterior');
            plot(ellip(:,1)+x0, ellip(:,2)+y0);
            axis([-70 70 -70 70]); 
            hold off;
        subplot(2,4,1)
            plot(f,Y,'-'); xlim([0 10]);
            xlabel('Frequency'); ylabel('Amplitude');
        subplot(2,4,5)
            plot(1:3600,   demeanedX_L(1:3600),   'b.'); hold on; 
            plot(3601:7200,demeanedX_L(3601:7200),'r.'); hold off; 
            xlabel('Time'); ylabel('Left-Right'); 
            ylim([-70,70]);
        subplot(2,4,6)
            plot(1:3600,   demeanedY_L(1:3600),   'b.'); hold on; 
            plot(3601:7200,demeanedY_L(3601:7200),'r.'); hold off; 
            %hold on; plot(mdl.Fitted); hold off
            xlabel('Time'); ylabel('Anterior-Posterior');
            legend('off');  ylim([-70,70]);
        end
    
        %% Then, with the RIGHT side data -------------------------------------------------------------
        side = 2;
        dd = d(:,10+trial*3:11+trial*3);
        dd_clean_R = filtfilt(B, A, dd);   % filter the input data
    
        % Demean the cleaned data
        demeanedX_R = dd_clean_R(:,1)-mean(dd_clean_R(:,1));
        demeanedY_R = dd_clean_R(:,2)-mean(dd_clean_R(:,2));
        % Flip the axis??
        demeanedX_R = -demeanedX_R; demeanedY_R = -demeanedY_R;
    
        % Compute time domain parameters for x-axis (LR) and y-axis (AP)
        [pk2pkx,meanDx,rmsDx,meanVelx,mfreqx,slopex] = time_domain(demeanedX_R, Fs);
        [pk2pky,meanDy,rmsDy,meanVely,mfreqy,slopey] = time_domain(demeanedY_R, Fs);
    
        % Compute 95% conf ellipse of COP
        [ellip,x0,y0,area] = error_ellipse([demeanedX_R demeanedY_R]);
    
        % Compute frequency domain parameters for x-axis (LR) and y-axis (AP)
        [powerX,p50X,p95X,spm1X,spm2X,cfreqX,freqdX,~,~] = freq_domain(demeanedX_R, Fs);
        [powerY,p50Y,p95Y,spm1Y,spm2Y,cfreqY,freqdY,f,Y] = freq_domain(demeanedY_R, Fs);
 
        % Adding those variables for quantifying COP
        COP_each_leg = [COP_each_leg; [subj,task,side,age(subj),meanDx,meanDy,rmsDx,rmsDy,...
                               pk2pkx,pk2pky,meanVelx,meanVely,slopex,slopey,mfreqx,mfreqy,area, ...
                               powerX,powerY,p50X,p50Y,p95X,p95Y,cfreqX,cfreqY,freqdX,freqdY ] ];
                       
        % Visualize the data on the same figure
        if(showPlot)
        subplot(2,4,3)
            plot(demeanedX_R(1:3600),   demeanedY_R(1:3600),   'b.'); hold on;  
            plot(demeanedX_R(3601:7200),demeanedY_R(3601:7200),'r.');
            title('RIGHT FOOT');
            xlabel('Left-Right'); ylabel('Anterior-Posterior');
            plot(ellip(:,1)+x0, ellip(:,2)+y0);
            axis([-70 70 -70 70]); 
            hold off;
        subplot(2,4,4)
            plot(f,Y,'-'); xlim([0 10]);
            xlabel('Frequency'); ylabel('Amplitude');
        subplot(2,4,7)
            plot(1:3600,   demeanedX_R(1:3600),   'b.'); hold on; 
            plot(3601:7200,demeanedX_R(3601:7200),'r.'); hold off; 
            xlabel('Time'); ylabel('Left-Right'); 
            ylim([-70,70]);
        subplot(2,4,8)
            plot(1:3600,   demeanedY_R(1:3600),   'b.'); hold on; 
            plot(3601:7200,demeanedY_R(3601:7200),'r.'); hold off; 
            %hold on; plot(mdl.Fitted); hold off
            xlabel('Time'); ylabel('Anterior-Posterior');
            legend('off');  ylim([-70,70]);
        end
    
        %% Left-Right time-series comparison: DTW, In-phase Synch, etc. --------------------------------
        dtw_x = dtw(demeanedX_L,demeanedX_R);
        dtw_y = dtw(demeanedY_L,demeanedY_R);
        dtw_R = dtw(sqrt(demeanedX_L.^2+demeanedY_L.^2),sqrt(demeanedX_R.^2+demeanedY_R.^2));
        mean_sync_x = instant_phase(demeanedX_L,demeanedX_R);
        mean_sync_y = instant_phase(demeanedY_L,demeanedY_R);
        mean_sync_R = instant_phase(sqrt(demeanedX_L.^2+demeanedY_L.^2),sqrt(demeanedX_R.^2+demeanedY_R.^2));
    
        COP_contrast = [COP_contrast; [subj,task,age(subj),dtw_x,dtw_y,dtw_R,...
                                       mean_sync_x,mean_sync_y,mean_sync_R] ];
    
    end
    pause;   % Pls comment this if you don't wanna wait!   

end
end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Convert matrix into table
myHeader = {'subj','task','side','age','MDx','MDy','RMSx','RMSy','pktopkX','pktopkY','Velx','Vely',...
            'slopeX','slopeY','mfreqx','mfreqy','area','powerX','powerY','p50X','p50Y','p95X',...
            'p95Y','cfreqX','cfreqY','freqdX','freqdY' };


%hist(COPtoSave.powerY, 60);
%COP_each_leg(COP_each_leg(:,19) > 300, 1:2)
%COP_each_leg(COP_each_leg(:,1) == 53, :)=[];

ssss = mean(COP_each_leg,1); ssss = array2table(ssss,'VariableNames',myHeader)

% Finally, save the subject-level metrics!
%COPtoSave = array2table(COP_each_leg,'VariableNames',myHeader);
%writetable(COPtoSave, 'C:\Users\ananda.sidarta\Documents\MATLAB\balance\balancestatic2.csv')


cd('C:\Users\ananda.sidarta\Documents\MATLAB')

    