% x = cos(pi/4*(0:99));
% y = hilbert(x);
% sigphase = ((angle(y)))';
% plot(x); hold on; plot(sigphase);

function mean_synch = instant_phase(x,y)
    % The following function computes the instantaneous phase between two input signals 
    % x(t) and y(t). First, we do the hilbert transform of each data, then estimate 
    % the angle (phase). Differences in their phases will be subsequently computed.
    
    x2 = hilbert(x);
    y2 = hilbert(y);
    sigphase_x2 = angle(x2)';   % Angle to be between {-pi,+pi}
    sigphase_y2 = angle(y2)';
    
    phase_synch = 1 - sin(abs(sigphase_y2-sigphase_x2)/2);
    %plot(phase_synch);hold on; plot(sigphase_x2);plot(sigphase_y2);
    mean_synch = mean(phase_synch);
    