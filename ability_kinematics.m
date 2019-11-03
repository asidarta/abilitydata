
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code compares kinematic trajectories of LEFT and RIGHT limbs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear; clc; 

% Define path and filename.
mypath = 'D:\extracted angle\SN005\';
%myfile = 'SN005_0006_Towel_5Jul19'
myfile = 'towel_cut_R.txt';
ntrial = 10;  %%%%%<<<<<<<<<<

% Skips the first three rows of data as table.
% Note this is similar to a open table or dataframe in R.
mydata = readtable(strcat(mypath,myfile), 'HeaderLines',0);  
% T.x = []   %%% Remove the first column
%mydata.time = [];
colNames = mydata.Properties.VariableNames;

% Data preparation: Loop across different joint angles. Our joint angles
% have been normalized to contain 200 time points per phase.
% The column naming follows the existing textfile I received...
for j = 1 : 9
  % I put some tweaks because of directional convention of Visual3D.
  if j == 1 || j == 4 || j == 7 
      R = mydata{:,13+j};   % Only for flex/extend, no need any tweaks!!
  else
      R = -mydata{:,13+j}; 
  end
  L = mydata{:,4+j};
  % Question: Is amplitude normalization within the same limb required 
  % here?? We are comparing L/R. It's good to see the plot of both first!
  %R = (R - mean(R))/std(R);
  %L = (L - mean(L))/std(L);ddc
  
  if(0)
      temp = strsplit( string(colNames(4+j)), '_');  % for plot naming  
      figure(gcf);%movegui(figure(1),'east')
      plot(R);  hold on;  plot(L);
  end
  
  % Reshape the time series for each elbow angle so that the column
  % represent different trials. This is to cater spm1d syntax.
  % Below, the arrays are in the form of [time x trial x joint]
  R_ang_reshaped(:,:,j) = reshape(R, [400 ntrial]);
  L_ang_reshaped(:,:,j) = reshape(L, [400 ntrial]);
  
  % We also get velocity as first derivative of position data.
  R_vel_reshaped(:,:,j) = reshape(diff([R;0]), [400 ntrial]);
  L_vel_reshaped(:,:,j) = reshape(diff([L;0]), [400 ntrial]);

end

% Compute mean and std angular trajectories across trials (2nd dim).
meanR_all = mean(R_ang_reshaped,2);
stdR_all  = std(R_ang_reshaped,[],2)/sqrt(ntrial);
meanL_all = mean(L_ang_reshaped,2);
stdL_all  = std(L_ang_reshaped,[],2)/sqrt(ntrial);

% Let's PLOT the angle you want (start from index 1,2,3 for elbow,
% 4,5,6 for shoulder, and 7,8,9 for wrist)
j = 6;  
movegui(figure(1),'east')
mean_r = meanR_all(:,:,j);  se_r   = stdR_all(:,:,j);
mean_l = meanL_all(:,:,j);  se_l   = stdL_all(:,:,j);
% Now plor both average curves (L/R) in one figure
%errorbar(mean_r,se_r,'bx'); hold on; errorbar(mean_l,se_l,'rx');
% Add standard error shading using patch function...
lo = mean_r - se_r;
hi = mean_r + se_r;
nsample = length(lo);
mycolor = 'r';   % red for Right -------------
hl = line(1:length(lo),mean_r,'color',mycolor);
hold on;
hp = patch([(1:nsample)'; (nsample:-1:1)'; 1], ...
           [lo; hi(nsample:-1:1); lo(1)], mycolor);
set(hp, 'facecolor', mycolor, 'edgecolor', 'none', 'FaceAlpha',.1);
lo = mean_l - se_l;
hi = mean_l + se_l;
nsample = length(lo);
mycolor = 'b';   % blue for Left -------------
hl = line(1:length(lo),mean_l,'color', mycolor);
hp = patch([(1:nsample)'; (nsample:-1:1)'; 1], ...
           [lo; hi(nsample:-1:1); lo(1)], mycolor);
set(hp, 'facecolor', mycolor, 'edgecolor', 'none', 'FaceAlpha',.1);
title('Mean angular trajectories (degree)');
xlabel('Normalized time-points'); ylabel('Average angle');
line([200 200],[min(mean_l)-10 max(mean_l)+10], 'color', 'black'); 



% PART 1: Using spm1d tools to compare................................
% Use spm1d to compare both ANGLE time-series.
for j = 1 : 9
  temp = strsplit( string(colNames(4+j)), 'L');   % for plot naming
  movegui(figure(2),'west')
  figure(gcf)
  sgtitle("Comparison of Angle (Top) & Velocity (Bottom) of RIGHT, LEFT")
  subplot(2,9,j);
  % Use spm1d to compare both angular time-series. Note input arrays format.
  t  = spm1d.stats.ttest_paired(L_ang_reshaped(:,:,j)', ...
                                R_ang_reshaped(:,:,j)');
  ti = t.inference(0.02, true);
  ti.plot();
  % Use spm1d also to compare both VELOCITY time-series.. 
  subplot(2,9,9+j)
  title(temp(2));
  t  = spm1d.stats.ttest_paired(R_vel_reshaped(:,:,j)', ...
                                L_vel_reshaped(:,:,j)');
  ti = t.inference(0.05, true);
  ti.plot()
end



% PART 2: Joint-space plotting. The idea is to plot joint angles in X,Y,Z
% system as given by the Virtual3D.
clist = colormap(hsv(400));
plot3(mydata.LElbowX, mydata.LElbowY, mydata.LElbowZ,'r.'); hold on;
plot3(mydata.LShoulderX, mydata.LShoulderY, mydata.LShoulderZ,'m.');
plot3(mydata.LWristX, mydata.LWristY, mydata.LWristZ,'b.');
xlabel('Flexion/Extension'); 
ylabel('Abduction/Adduction'); 
zlabel('Rotation');





% PART 3: Cross-correlation analyses...................................
% Do it on joint angle trajectory data. Here we test similarity between 
% 2 time-series (L/R) as a function of lag or delay.
%figure(j+1); figure(gcf)
%[c,lags] = xcorr(R, L, '0', 'coeff');
%title(temp(1));
%stem(lags, c);
  
% Note: Are tweaks and normalization steps still needed here?!
% One possible gotcha in cross correlation is that a DC bias in the waveform 
% will corrupt the result.

% Can we plot cross-correlation diagram between LEFT and RIGHT kinematics?
% We can use the mean angular trajectories. Note: Use squeeze function!!

x_corr = corr( squeeze(meanL_all(200:400,:)), ...
               squeeze(meanR_all(200:400,:)) );

x_corr(find(abs(x_corr) < 0.6)) = 0;
myLabel = colNames(:,5:22);  % what are the joint names?
imagesc(x_corr)  % produce the heatmap
set(gca, 'XTick', 1:length(myLabel)); % center x-axis ticks on bins
set(gca, 'YTick', 1:length(myLabel)); % center y-axis ticks on bins
set(gca, 'XTickLabel', colNames(:, 5:13)); % set x-axis labels
set(gca, 'YTickLabel', colNames(:,14:22)); % set y-axis labels
xtickangle(45);
axis square; colorbar
title('Cross-corr between Left and Right', 'FontSize', 14); % set title
colormap('jet'); % set the colorscheme

% Cross corr between Right & Left of the same joint angles
x_corr = corr( squeeze(meanL_all(1:60,:)), squeeze(meanR_all(1:60,:)) )




% PART 4: Frequency spectrum analysis of time series as whole-trial!
% This may not give more information than the temporal ones.
% NOTE = in this case I believe it is not good to use the cut time-series
% as you want to capture variability in time also.
var1 = mydata{:,9}; 
Fs = 200;                % Sampling frequency                    
T  = 1/Fs;               % Sampling period       
L  = length(var1);       % Length of signal
t  = (0:L-1)*T;          % Time vector
Y  = fft(squeeze(var1)); % Do FFT on the amplitude data
YL = abs(Y/L);
P1 = YL(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;      % Freq-range for X-axis, folded.
plot(f,P1)               % Plot the single-sided spectrum 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)'); ylabel('|P1(f)|')
xlim([0,15]); ylim([-2 20]);

hold on;

var1 = mydata{:,12}; 
Fs = 200;                % Sampling frequency                    
T  = 1/Fs;               % Sampling period       
L  = length(var1);       % Length of signal
t  = (0:L-1)*T;          % Time vector
Y  = fft(squeeze(var1)); % Do FFT on the amplitude data
YL = abs(Y/L);
P2 = YL(1:L/2+1);
P2(2:end-1) = 2*P2(2:end-1);
f = Fs*(0:(L/2))/L;      % Freq-range for X-axis, folded.
plot(f,P2)               % Plot the single-sided spectrum 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)'); ylabel('|P2(f)|')
xlim([0,15]); ylim([-2 20]);

%%%% Cross spectrum analysis, measuring coherence %%%%%
% PSD: distribution of power along the frequency axis.
% CSD: similar to cross-correlation, so can find the power shared by a 
% given frequency for the two signals using its squared module, and the 
% phase shift between the two signals at that frequencyt.
%hpsd = dspdata.psd(P2,'Fs',Fs); plot(hpsd);
% MSC: magnitude squared coherence, indicate how well two signals are 
% correlated at each frequency. It's measured by the magnitude of
% cross spectral density between two signals over spectral density of 
% the individual signal.
% Compute MSC by doing element-wise operation below.
numer = (abs(cpsd(P1,P2))).^2;
denom = abs(cpsd(P1,P1)).* abs(cpsd(P2,P2));
coherence = numer ./ denom; 
plot(coherence, 'b'); hold on;
mscohere(P1,P2);  % or, using the function


% PART 5: Angle-angle curve. This is to visualize coordination between
% two joint angles........................................
var1 = mydata{:,14};  % Left elbowX
%var1 = mydata{:,8};
var2 = mydata{:,5}; % Right elbowX
%var1 = mydata{:,17};
var1 = reshape(var1, [400 ntrial]);
var2 = reshape(var2, [400 ntrial]);
plot(var1, '.'); hold on; plot(var2, '*');
plot(var1(:,1), '.'); hold on; plot(var2(:,2), '*');

clist = colormap(hsv(400));
figure(3); hold on;
title('Spatial coordination between Left and Right');
for trial = 1:ntrial   
    for time = 1:400
      xlim([0,100]); ylim([0,100]);
      plot(var1(time,trial), var2(time,trial), '.', 'color', clist(time,:));
      %pause(0.05)
    end
    %close all;
end






% PART 5: Using radar plot. First I compute the correlation of joint angles
% between two adjacent joints, e.g. wrist and elbow, elbow and shoulder.
% Then I plot the multiple correlation values using radar plot.

% LET US START BY LOOKING AT ELBOW AND SHOULDER. 
% Start by looking at the first folding stage!
to_compare = cat(3, R_ang_reshaped(1:100,:,1:6), ...
                    L_ang_reshaped(1:100,:,1:6));
for angle = 1:6
    for j = 1:12
        for trial = 1:ntrial
            corr_trial(trial) = corr(to_compare(:,trial,angle), ... 
                                     to_compare(:,trial,j));
            corr_mean(angle,j) = mean(corr_trial);
        end
    end
    % Plot according to this joint sequence: ELBOW, SHOULDER
    figure(angle)
    radarplot( 1+corr_mean(angle,:), ...
              {'R-ElbowX','R-ElbowY','R-ElbowZ',...
               'R-ShoulderX','R-ShoulderY','R-ShoulderZ',...
               'L-ElbowX','L-ElbowY','L-ElbowZ',...
               'L-ShoulderX','L-ShoulderY','L-ShoulderZ'});
end

% Then the SECOND folding stage
to_compare = cat(3, R_ang_reshaped(101:200,:,1:6), ...
                    L_ang_reshaped(101:200,:,1:6));
for angle = 1:6
    for j = 1:12
        for trial = 1:ntrial
            corr_trial(trial) = corr(to_compare(:,trial,angle), ... 
                                     to_compare(:,trial,j));
            corr_mean(angle,j) = mean(corr_trial);
        end
    end
    % Plot according to this joint sequence: ELBOW, SHOULDER
    figure(angle)
    radarplot( 1+corr_mean(angle,:), ...
              {'R-ElbowX','R-ElbowY','R-ElbowZ',...
               'R-ShoulderX','R-ShoulderY','R-ShoulderZ',...
               'L-ElbowX','L-ElbowY','L-ElbowZ',...
               'L-ShoulderX','L-ShoulderY','L-ShoulderZ'});
end


% Now between ELBOW and WRIST --------------------------------------------
to_compare = cat(3, R_ang_reshaped(1:100,:,[1:3,7:9]), ...
                    L_ang_reshaped(1:100,:,[1:3,7:9]));
for angle = 1:6
    for j = 1:12
        for trial = 1:ntrial
            corr_trial(trial) = corr(to_compare(:,trial,angle), ... 
                                     to_compare(:,trial,j));
            corr_mean(angle,j) = mean(corr_trial);
        end
    end
    % Plot according to this joint sequence: ELBOW, SHOULDER
    figure(angle)
    radarplot( 1+corr_mean(angle,:), ...
              {'R-ElbowX','R-ElbowY','R-ElbowZ',...
               'R-WristX','R-WristY','R-WristZ',...
               'L-ElbowX','L-ElbowY','L-ElbowZ',...
               'L-WristX','L-WristY','L-WristZ'});
end

to_compare = cat(3, R_ang_reshaped(101:200,:,[1:3,7:9]), ...
                    L_ang_reshaped(101:200,:,[1:3,7:9]));
for angle = 1:6
    for j = 1:12
        for trial = 1:ntrial
            corr_trial(trial) = corr(to_compare(:,trial,angle), ... 
                                     to_compare(:,trial,j));
            corr_mean(angle,j) = mean(corr_trial);
        end
    end
    % Plot according to this joint sequence: ELBOW, SHOULDER
    figure(angle)
    radarplot( 1+corr_mean(angle,:), ...
              {'R-ElbowX','R-ElbowY','R-ElbowZ',...
               'R-WristX','R-WristY','R-WristZ',...
               'L-ElbowX','L-ElbowY','L-ElbowZ',...
               'L-WristX','L-WristY','L-WristZ'});
end




