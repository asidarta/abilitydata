

function [pk2pk,meanD,rmsD,meanVel,mfreq,totex,slope,zcross] = time_domain(tdata,Fs)

    % Compute time domain variables to quantify centre of pressure (COP) in a particular direction
    % Inputs   : time-series (tdata) and sampling freq (Fs)
    % Outputs  : several parameters in the time domain

    pk2pk = max(tdata)-min(tdata);          % Peak-to-peak magnitude
    meanD = sum(abs(tdata))/length(tdata);  % Mean distance (MD)
    rmsD  = std(tdata);                     % RMS distance
    velx  = abs(gradient(tdata,1/Fs));      % Note the absolute here, take the magnitude
    meanVel = mean(velx);                   % Mean velocity (MV)
    mfreq = meanVel / (4*pi*meanD);         % Median frequency (MFREQ)
    totex = sum( abs(diff(tdata)) );        % Total excursion

    mdl = fitlm(1:length(tdata),tdata);              % fit a linear reg line
    slope = 1000*abs(mdl.Coefficients.Estimate(2));  % get the abs linear slope
    
    zcross = sum(diff(sign(tdata))~= 0);    % Counting how many mean/zero crossing
    
