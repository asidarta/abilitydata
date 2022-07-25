
clear; clc; close all
%%% Subject #11,#12,#53 is bad!!!

%%
% Define the working directory. 
wd = 'C:\Users\ananda.sidarta\Documents\MATLAB\balance\cop_balance_static\';
cd(wd);

COP_each_leg = [];   % COP metrics for each leg separately
COP_contrast = [];   % COP metrics to contrast L/R side
meanSpect = [];
showPlot  = 1;       % Want to show the figure???? ------------------

age = [3,6,2,2,4,3,5,2,5,2,6,6,5,2,3,6,5,5,6,3,3,4,5,3,4,4,4,5,2,6,5,5,3,2,3,6,4,3,5,5,3,5,5,7,5,6,5,6,...
       5,2,5,4,6,7,6,5,5,6,7,6,4,6,3,5,7,4,5,2,7,5,7,6,6,6,6,5,6,7,7,5,5,5,3,2,6,5,5,6,5,5,6,6,6,6,6,6,7,5,2,6];

   
   
for subj = 42:42    % Looping for EVERY subject (numbered 1 - 100)
close all; subj   % Which subject are we at now?
    
for task = 1:2      % Looping for EVERY task (closed=1, open=2)

    if task == 1,   taskStr = 'close';   % String handler for task
    else,  taskStr = 'open';    
    end
    figure(task);  % Open a figure according to task EC/EO
    
    % Now open the relevant file based on subject and task type
    myfile = strcat(taskStr, 'eyes_cofp_',  num2str(subj) );
    
    % Then load the textfile as Table. Please skip the first FIVE rows!
    T = readtable(myfile,'HeaderLine',5,'ReadVariableNames',0);%,'Format','%f%f%f%f%f%f%f');
    T(:,1)=[];              % Delete the first column!
    T = T(:,1:6);           % If Force-plate 3 exists, delete the data! Not used.
    d = table2array(T);     % Convert to numeric ARRAY
    d = -1000*d;            % Convert to millimeter(mm) unit. INVERSE??
    Fs = 7200/30;           % Sampling frequency of the recording
    
    % Do you want to focus on transition between tasks???????? (closed=1, open=2)
    %if task == 1,   d = d(1:2400,:);   % Close-eye, early part
    %else,  d = d(4800:7200,:);         % Open-eye, last part    
    %end
    
    % Design a 4th order Butterworth Low Pass Filter, cutoff freq = 5Hz
    cutfreq = 5/(Fs/2); 
    [B,A] = butter(4,cutfreq);
    
    %% Deal with the LEFT side data first ---------------------------------------------------------
    side = 1;
    dd = d(:,1:2);
    dd_clean_L = filtfilt(B, A, dd);   % filter the input data, subscript L for left!
    
    % Demean the cleaned data
    demeanedX_L = dd_clean_L(:,1)-mean(dd_clean_L(:,1));
    demeanedY_L = dd_clean_L(:,2)-mean(dd_clean_L(:,2));
    
    % Compute time domain parameters for x-axis (LR) and y-axis (AP)
    [pk2pkx,meanDx,rmsDx,meanVelx,mfreqx,totexX,slopex,zcrossX] = time_domain(demeanedX_L, Fs);
    [pk2pky,meanDy,rmsDy,meanVely,mfreqy,totexY,slopey,zcrossY] = time_domain(demeanedY_L, Fs);
    
    % Compute fractal dimension
    totex  = sqrt(totexX^2 + totexY^2);  % total excursion
    fractd = log(size(dd,1))/log(size(dd,1)*max(pk2pkx,pk2pky)/totex);
    
    % Compute 95% conf ellipse of COP
    [ellip,x0,y0,area] = error_ellipse([demeanedX_L demeanedY_L]);
    
    % Compute frequency domain parameters for x-axis (LR) and y-axis (AP)
    [powerX,p50X,p95X,spm1X,spm2X,cfreqX,freqdX,~,~] = freq_domain(demeanedX_L, Fs);
    [powerY,p50Y,p95Y,spm1Y,spm2Y,cfreqY,freqdY,f,Y] = freq_domain(demeanedY_L, Fs);
 
    % Compile ALL variables for quantifying COP
    COP_each_leg = [COP_each_leg; [subj,task,side,age(subj),meanDx,meanDy,rmsDx,rmsDy,pk2pkx,pk2pky,...
                               meanVelx,meanVely,slopex,slopey,mfreqx,mfreqy,totex,zcrossX,zcrossY,area,...
                               fractd,powerX,powerY,p50X,p50Y,p95X,p95Y,cfreqX,cfreqY,freqdX,freqdY ] ];
                           
    
    % Visualize the data, for LEFT foot
    if(showPlot)
    fff = figure(task);
    fff.WindowState = 'maximized';
    fff.Name = strcat('COP_trajectory_Eyes_',taskStr);
    subplot(2,4,2)
        plot(demeanedX_L,demeanedY_L,'k.');
        hold on;  title('LEFT FOOT');
        xlabel('Left-Right'); ylabel('Anterior-Posterior');
        plot(ellip(:,1)+x0, ellip(:,2)+y0);
        axis([-20 20 -20 20]); 
        hold off;
    subplot(2,4,1)
        plot(f,Y,'-'); xlim([0 3]);
        xlabel('Frequency'); ylabel('Amplitude');
    subplot(2,4,5)
        plot(demeanedX_L,'b.'); 
        xlabel('Time'); ylabel('Left-Right'); 
        ylim([-30,30]);
    subplot(2,4,6)
        plot(demeanedY_L,'b.'); %hold on; plot(mdl.Fitted); hold off
        xlabel('Time'); ylabel('Anterior-Posterior');
        legend('off');  ylim([-30,30]);
    end
    
    %% Then, with the RIGHT side data -------------------------------------------------------------
    side = 2;
    dd = d(:,4:5);
    dd_clean_R = filtfilt(B, A, dd);   % filter the input data, subscript R for left!
    
    % Demean the cleaned data
    demeanedX_R = dd_clean_R(:,1)-mean(dd_clean_R(:,1));
    demeanedY_R = dd_clean_R(:,2)-mean(dd_clean_R(:,2));
    
    % Compute time domain parameters for x-axis (LR) and y-axis (AP)
    [pk2pkx,meanDx,rmsDx,meanVelx,mfreqx,totexX,slopex,zcrossX] = time_domain(demeanedX_R, Fs);
    [pk2pky,meanDy,rmsDy,meanVely,mfreqy,totexY,slopey,zcrossY] = time_domain(demeanedY_R, Fs);
    
    % Compute fractal dimension
    totex  = sqrt(totexX^2 + totexY^2);  % total excursion
    fractd = log(size(dd,1))/log(size(dd,1)*max(pk2pkx,pk2pky)/totex);
    
    % Compute 95% conf ellipse of COP
    [ellip,x0,y0,area] = error_ellipse([demeanedX_R demeanedY_R]);

    % Compute frequency domain parameters for x-axis (LR) and y-axis (AP)
    [powerX,p50X,p95X,spm1X,spm2X,cfreqX,freqdX,~,~] = freq_domain(demeanedX_R, Fs);
    [powerY,p50Y,p95Y,spm1Y,spm2Y,cfreqY,freqdY,f,Y] = freq_domain(demeanedY_R, Fs);
 
    % Adding those variables for quantifying COP
    COP_each_leg = [COP_each_leg; [subj,task,side,age(subj),meanDx,meanDy,rmsDx,rmsDy,pk2pkx,pk2pky,...
                               meanVelx,meanVely,slopex,slopey,mfreqx,mfreqy,totex,zcrossX,zcrossY,area,...
                               fractd,powerX,powerY,p50X,p50Y,p95X,p95Y,cfreqX,cfreqY,freqdX,freqdY ] ];
                       
    % Visualize the data on the same figure for RIGHT foot
    if(showPlot)
    subplot(2,4,3)
        plot(demeanedX_R,demeanedY_R,'k.'); title('RIGHT FOOT');
        hold on;
        xlabel('Left-Right'); ylabel('Anterior-Posterior');
        plot(ellip(:,1)+x0, ellip(:,2)+y0);
        axis([-20 20 -20 20]); 
        hold off;
    subplot(2,4,4)
        plot(f,Y,'-'); xlim([0 3]);
        xlabel('Frequency'); ylabel('Amplitude');
    subplot(2,4,7)
        plot(demeanedX_R,'b.'); 
        xlabel('Time'); ylabel('Left-Right'); 
        ylim([-30,30]);
    subplot(2,4,8)
        plot(demeanedY_R,'b.'); %hold on; plot(mdl.Fitted); hold off
        xlabel('Time'); ylabel('Anterior-Posterior');
        legend('off');  ylim([-30,30]);
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
%pause;   % Pls comment this if you don't wanna wait!
end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Convert matrix into table
myHeader = {'subj','task','side','age','MDx','MDy','RMSx','RMSy','pktopkX','pktopkY','Velx','Vely',...
            'slopeX','slopeY','mfreqx','mfreqy','totex','zcrossX','zcrossY','area','fractd','powerX','powerY',...
            'p50X','p50Y','p95X','p95Y','cfreqX','cfreqY','freqdX','freqdY' };


%hist(COPtoSave.area, 60);
%COP_each_leg(COP_each_leg(:,19) > 300, 1:2)
%COP_each_leg(COP_each_leg(:,1) == 53, :)=[];

%ssss = mean(COP_each_leg,1); ssss = array2table(ssss,'VariableNames',myHeader)

% Finally, save the subject-level metrics!
COPtoSave = array2table(COP_each_leg,'VariableNames',myHeader);
%writetable(COPtoSave, 'C:\Users\ananda.sidarta\Documents\MATLAB\balance\balancestatic.csv')


cd('C:\Users\ananda.sidarta\Documents\MATLAB')

  