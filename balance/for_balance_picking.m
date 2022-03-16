
clear; clc;
%close all

%%
% Define the working directory. 
wd = 'C:\Users\ananda.sidarta\OneDrive - Nanyang Technological University\ananda\For Lester\Ananda\outputV3D\bySubject\balance\force';
cd(wd);
time_COP  = [];
meanSpect = [];

%for subj = 1:60    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ensure all files are of type csv (text delimited) so that we can read table as Numeric.
myfile = strcat('balance_R_cofp_15');%, num2str(subj) );

% Then load the textfile as Table. Please skip the first FIVE rows!
T = readtable(myfile,'HeaderLine',5,'ReadVariableNames',0);
T(:,1)=[];              % Delete the first column!
d = table2array(T);     % Convert to numeric ARRAY
Fs = 2000;              % sampling frequency

% Decide which Force Plate (Left or right). Convert to Centimeter unit.
dd = 100*d(:,1:2); 

% Compute summary statistics (Mean, std) of the data
meanX = mean(dd(:,1)); 
meanY = mean(dd(:,2));
demeanedX = dd(:,1)-meanX;
demeanedY = dd(:,2)-meanY;
pktopkX = max(dd(:,1))-min(dd(:,1));
pktopkY = max(dd(:,2))-min(dd(:,2));

% Compute Time-domain variables to quantify COP -------------------------------
meandistX = sum(abs(demeanedX))/length(demeanedX);  % Mean distance X (MDx)
meandistY = sum(abs(demeanedY))/length(demeanedY);  % Mean distance Y (MDy)
rmsdistX  = std(demeanedX);
rmsdistY  = std(demeanedY);
velX = gradient(demeanedX,1/Fs);
velY = gradient(demeanedY,1/Fs);
meanVelX  = mean(velX);    % Mean velocity X (MVx)
meanVelY  = mean(velY);    % Mean velocity Y (MVy)
mfreqX = meanVelX / (4*sqrt(2)*meandistX);
mfreqY = meanVelY / (4*sqrt(2)*meandistY);

% Compute 95% conf ellipse of COP
[ellip,x0,y0,area] = error_ellipse(dd);

% Compile the time-domain variables for quantifying COP
time_COP = [time_COP; [meandistX meandistY rmsdistX rmsdistY pktopkX pktopkY ...
                       meanVelX meanVelY mfreqX mfreqY area]];

% Compute Freq-domain variables to quantify COP -------------------------------

% Obtain the spectrum(AP axis), Fs = sampling freq
[a,b] = fft_matlab(demeanedY,Fs);
meanSpect = [meanSpect; b(1:500)']; 


%end     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean(time_COP,1)

% See the data!
figure(1)
subplot(2,2,1)
  plot(dd(:,1),dd(:,2),'k.'); hold on;
  xlabel('Left-Right'); ylabel('Anterior-Posterior');
  plot(ellip(:,1)+x0, ellip(:,2)+y0);
  axis([meanX-2,meanX+2 meanY-2,meanY+2]);
subplot(2,2,2)
  plot(a,b,'-');
  axis([0 30 0 0.5]); 
subplot(2,2,3)
  plot(dd(:,1)); ylabel('Left-Right'); 
  hold on;
  ylim([meanX-1,meanX+1]);
subplot(2,2,4)
  plot(dd(:,2)); ylabel('Anterior-Posterior');
  ylim([meanY-1,meanY+1]);
  
figure(2)
  plot(mean(meanSpect,1))

  