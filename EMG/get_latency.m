%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 2 Aug 2022. Last revision: N/A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function latencyout = get_latency (mydata, ref, Fs, emglabel)
% This function extracts latency of an EMG dataset w.r.t perturbation onset!
% The input is the cleaned EMG data and a reference sample point (perturbation onset!)
% Fs = sample rate;  peak500 = max amplitude within 500 msec post-perturbation


% Excitation detector: We assume activity that occurs way before the perturbation to be 
% the baseline. If muscle activity is higher than baseline, then the muscle is excited.
% Excitation = activity > baseline + 2*SD of baseline;
% Or excitation = baseline + 10 % of (maximum peak - baseline)

try
    % Baseline activity is taken as average activity prior to perturbation
    baseline = mean(mydata(1:500));
    datapeak = max(mydata) - baseline;
    % Start of muscle excitation is w.r.t to the perturbation onset, e.g. excite = 230 
    % is the 230th element post perturbation. Note I also added an offset, such that 
    % the start of excitation at least > 20 msec post perturbation onset.
    excite = find(mydata((ref+40):length(mydata)) > baseline+0.1*datapeak,1);

    % Extract latency and peak
    latency_in_sample = excite+40;
    latency_in_time = latency_in_sample/Fs;
    peak500 = max(mydata(latency_in_sample: latency_in_sample+ 0.5*Fs)) - baseline;
    latency_in_time = latency_in_sample/Fs;
catch
    latency_in_time = 0;
    latency_in_sample = 0;
    peak500 = 0;
    warning('!! No latency found. This is not perturbation trial');
end

% Pack the output information as struct
latencyout = struct('Muscle', emglabel,...
                    'Latency', latency_in_time,...
                    'Latency_sample', latency_in_sample,...
                    'Peak500', peak500);

end
