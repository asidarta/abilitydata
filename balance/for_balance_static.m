
clear; clc; close all
%%% Subject #53 is bad!!!

%%
% Define the working directory. 
wd = 'C:\Users\ananda.sidarta\Documents\MATLAB\balance\cop_balance_static\';
cd(wd);
COPmetrics  = [];
meanSpect = [];
age = [3,6,2,2,4,3,5,2,5,2,6,6,5,2,3,6,5,5,6,3,3,4,5,3,4,4,4,5,2,6,5,5,3,2,3,6,4,3,5,5,3,5,5,7,5,6,5,6,...
       5,2,5,4,6,7,6,5,5,6,7,6,4,6,3,5,7,4,5,2,7,5,7,6,6,6,6,5,6,7,7,5,5,5,3,2,6,5,5,6,5,5,6,6,6,6,6,6,7,5,2,6];

   
   
for subj = 1:100   % Looping for EVERY subject (numbered 1 - 100)
for task = 1:2     % Looping for EVERY task (closed=1, open=2)
for side = 1:2     % Looping for EACH forceplate (left=1, right=2)

    if task == 1,   taskStr = 'close';   % String handler for task
    else,  taskStr = 'open';    
    end
    if side == 1,  sideStr = 'LEFT';    % String handler for FP side
    else,  sideStr = 'RIGHT';    
    end
    
    % Now open the relevant file based on subject and task type
    myfile = strcat(taskStr, 'eyes_cofp_',  num2str(subj) );
    
    % Then load the textfile as Table. Please skip the first FIVE rows!
    T = readtable(myfile,'HeaderLine',5,'ReadVariableNames',0);%,'Format','%f%f%f%f%f%f%f');
    T = T(:,1:6);           % If Force-plate 3 exists, delete the data! Not used.
    T(:,1)=[];              % Delete the first column!
    d = table2array(T);     % Convert to numeric ARRAY
    d = 1000*d;             % Convert to millimeter(mm) unit.
    Fs = 2000;              % Sampling frequency of the recording

    % Decide which Force Plate (Left or right). 
    if side == 1,  dd = d(:,1:2); 
    else,          dd = d(:,4:5); 
    end
    
    % Design a Low Pass Filter
    cutfreq = 5/(Fs/2);            % provide cutoff freq = 5Hz
    [B,A] = butter(4,cutfreq);     % design 4th order butterworth filter
    dd_clean = filtfilt(B,A,dd);   % do filtering now!! 
    
    % Demean the cleaned data
    demeanedX = dd_clean(:,1)-mean(dd_clean(:,1));
    demeanedY = dd_clean(:,2)-mean(dd_clean(:,2));
    
    % Compute time domain parameters for x-axis (LR) and y-axis (AP)
    [pk2pkx,meanDx,rmsDx,meanVelx,mfreqx,slopex] = time_domain(demeanedX, Fs);
    [pk2pky,meanDy,rmsDy,meanVely,mfreqy,slopey] = time_domain(demeanedY, Fs);
    
    % Compute 95% conf ellipse of COP
    [ellip,x0,y0,area] = error_ellipse([demeanedX demeanedY]);
    
    % Compute frequency domain parameters for x-axis (LR) and y-axis (AP)
    [powerX,p50X,p95X,spm1X,spm2X,cfreqX,freqdX,~,~] = freq_domain(demeanedX, Fs);
    [powerY,p50Y,p95Y,spm1Y,spm2Y,cfreqY,freqdY,f,Y] = freq_domain(demeanedY, Fs);
    
    %meanSpect = [meanSpect; b(1:500)']; 
  
    % Compile ALL variables for quantifying COP
    COPmetrics = [COPmetrics; [subj,task,side,age(subj),meanDx,meanDy,rmsDx,rmsDy,...
                           pk2pkx,pk2pky,meanVelx,meanVely,slopex,slopey,mfreqx,mfreqy,area, ...
                           powerX,powerY,p50X,p50Y,p95X,p95Y,cfreqX,cfreqY,freqdX,freqdY ] ];

    % Visualize the data!
    if(0) fff = figure(side);
    fff.Name = strcat('COP from the_',sideStr,'_Force Plate; Eyes_',taskStr);
    subplot(2,2,1)
        plot(demeanedX,demeanedY,'k.'); hold on;
        xlabel('Left-Right'); ylabel('Anterior-Posterior');
        plot(ellip(:,1)+x0, ellip(:,2)+y0);
        axis([-30 30 -30 30]); 
        hold off;
    subplot(2,2,2)
        plot(f,Y,'-'); axis([0 15 0 50]);
        xlabel('Frequency'); ylabel('Amplitude');
    subplot(2,2,3)
        plot(demeanedX,'b.'); 
        xlabel('Time'); ylabel('Left-Right'); 
        ylim([-30,30]);
    subplot(2,2,4)
        plot(demeanedY,'b.'); %hold on; plot(mdl.Fitted); hold off
        xlabel('Time'); ylabel('Anterior-Posterior');
        legend('off');
        ylim([-30,30]);
    end
end
end
%pause;   % Pls comment this if you don't wanna wait!
end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Convert matrix into table
myHeader = {'subj','task','side','age','MDx','MDy','RMSx','RMSy','pktopkX','pktopkY','Velx','Vely',...
            'slopeX','slopeY','mfreqx','mfreqy','area','powerX','powerY','p50X','p50Y','p95X',...
            'p95Y','cfreqX','cfreqY','freqdX','freqdY' };


%hist(COPtoSave.powerY, 60);
COPmetrics(COPmetrics(:,19) > 300, 1:2)
COPmetrics(COPmetrics(:,1) == 53, :)=[];
%ssss = mean(COPmetrics,1); COPtoSave = array2table(ssss,'VariableNames',myHeader)

% Finally, save the subject-level metrics!
%COPtoSave = array2table(COPmetrics,'VariableNames',myHeader);
%writetable(COPtoSave, 'C:\Users\ananda.sidarta\Documents\MATLAB\balance\balancestatic2.csv')


cd('C:\Users\ananda.sidarta\Documents\MATLAB')

  