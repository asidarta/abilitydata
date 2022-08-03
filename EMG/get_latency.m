%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created: 2 Aug 2022. Last revision: N/A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function latency_in_time = get_latency (mydata, ref, Fs)
% This function extracts latency of an EMG dataset w.r.t perturbation onset!
% The input is the cleaned EMG data and a reference sample point (perturbation onset!)
% Fs = sample rate


% Excitation detector: We assume activity that occurs way before the perturbation to be 
% the baseline. If muscle activity is higher than baseline, then the muscle is excited.
% Excitation = activity > baseline + 2*SD of baseline;
% Or excitation = baseline + 10 % of (maximum peak - baseline)

try
    baseline = mean(mydata(1:100));
    datapeak = max(mydata) - baseline;
    excite = find(mydata(ref:length(mydata)) > baseline+0.1*datapeak );

    excite(1);
    latency_in_samples = excite(1);
    latency_in_time = latency_in_samples/Fs;
catch
    latency_in_time = 0;
    warning('!! No latency found. This is not perturbation trial');
end

end
