

function [powerX,p50,p95,spm1,spm2,cfreq,freqd,f,Y] = freq_domain(tdata, Fs)

    % Compute Freq-domain variables to quantify time-varying centre of pressure (COP)
    % Inputs   : time-series (tdata) and sampling freq (Fs)
    % Outputs  : several parameters in the frequency domain

    [f,Y,psdX] = fft_matlab(tdata,Fs);         % Obtain LR spectrum (X axis), sample freq Fs
    
    endd = find(f>=5,1);
    dff = f(3:endd);                           % Disregarding 1st,2nd,3rd freq content, up to 5Hz
    powerX = sum(psdX(3:endd));                % Total power = sum of PSD, up to 5Hz only...   
    spm1 = sum(psdX(3:endd).*(dff)');          % 1st Spectral moment, up to 5Hz
    spm2 = sum(psdX(3:endd).*(dff.^2)');       % 2nd Spectral moment, up to 5Hz
    cfreq = sqrt(spm2/powerX);                 % Centroidal frequency (Hz)
    freqd = sqrt(1-(spm1^2/(powerX*spm2)));    % Frequency dispersion, unitless

    p50 = f( find(cumsum(psdX)>0.50*powerX,1) );    % Freq, below which 50% power is contained    
    p95 = f( find(cumsum(psdX)>0.95*powerX,1) );    % Freq, below which 95% power is contained

        
    