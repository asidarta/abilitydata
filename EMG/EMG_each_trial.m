
% Do you want to see latency for each trial?

Nelectr = 16;  % How many EMG electrodes?
Fs  = 2000;
ref = 0.25*Fs;

ncol = size(mysegment_out,2);

for mmm = 1:Nelectr    % For each electrode
    for f = 1:ncol
        % For each column is each file in the list, each row is each muscle. 
        % Use mysegment_out instead of mysegment!
        temp = mysegment_out(mmm,f);
        temp = temp{:};   % Convert cell to array

        % For each trial of each muscle electrode....
        for zz = 1: size(temp,2)
            baseline = mean(temp(1:500,zz));
            datapeak = max(temp(:,zz)) - baseline;
        
            % Start of muscle excitation is w.r.t to the perturbation onset, e.g. excite = 230 
            % is the 230th element post perturbation. Note I also added an offset, such that 
            % the start of excitation at least > 20 msec post perturbation onset.
            excite = find(temp((ref+40):length(temp(:,zz)),zz) > baseline+0.1*datapeak,1);

            % Extract latency and peak
            latency_in_sample(mmm,zz) = excite+40;
            latency_in_time(mmm,zz) = latency_in_sample(mmm,zz)/Fs;
            peak500(mmm,zz) = max(temp(excite+40: excite+40+ 0.5*Fs)) - baseline;
            latency_in_time(mmm,zz) = latency_in_sample(mmm,zz)/Fs;
        end 
    end        
end