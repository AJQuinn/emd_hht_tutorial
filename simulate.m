function obj = simulate(S)
%%function obj = cfc.simulate(S)
% Generate a simulated signal with a known cross frequency coupling
%
% Abridged version of the function available in the CFC toolbox
% https://github.com/AJQuinn/cfc
% https://github.com/AJQuinn/cfc/blob/master/%2Bcfc/simulate.m
%
% Requried fields in S
% 
% S.seconds            - Duration of signal
% S.sample_rate        - Sample rate of signal
% S.modulating_freq    - Frequency of low frequency modulating signal (in Hz)
% S.modulating_amp     - Amplitude of low frequency modulating signal (default = 2)
% S.modulated_freq     - Frequency of high frequency modulating signal (in Hz)
% S.modulated_amp      - Amplitude of high frequency modulating signal (default = 1)
% S.phase_lag          - Phase lag between modulating signal and modulated envelope (default = 0)
% S.method             - Method for generating signal (Default = 'basic')
% S.noise_ratio        - Noise in dB (default = -40)

%% set defaults
if ~isfield(S,'modulating_amp')
   S.modulating_amp = 2; 
end

if ~isfield(S,'modulated_amp')
   S.modulated_amp = 1; 
end

if ~isfield(S,'phase_lag')
   S.phase_lag = 0; % Lock to peak of theta wave
end

if ~isfield(S,'noise_level')
    S.noise_level = -50;
end

if ~isfield(S,'method')
   S.method = 'basic';
end

%% Housekeeping
obj.time_vect = [0:1/S.sample_rate:S.seconds];
obj.sr = S.sample_rate;

if isfield(S, 'switching_freq')
    obj.switching_freq = S.switching_freq;
    state_switching = square(sin(2*pi*S.switching_freq*obj.time_vect));
else
    S.switching_freq = 0;
    state_switching = zeros(size(obj.time_vect));
end

if ~isfield(S,'randseed')
    rng('shuffle');
end

%% Development note
%
% Each function below should populate obj with a modulated and a modulating
% time-series along with any relevant state switching. These will be combined
% along with any scaled noise at the end.
%
% obj.modulating_ts
% obj.modulated_ts
%
% Any additional arguments required for the function should be commented here

%% Main loop
if strcmp(S.method,'none')
        % two uncoupled signals, tonic low and high frequency oscillations

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        % ensure this is tonic
        state_switching = ones(size(modulating_ts,1),size(modulating_ts,2));

        obj.modulating_ts = obj.modulating_ts.*state_switching;

elseif strcmp(S.method,'basic')
        % A simple phase amplitude coupling between low and high frequency signals

        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);
        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'Abreu2010')
        % An analytical formulation of a nonlinear near-bed wave
        %
        % The degree of nonlinearity is specified by S.nonlin_deg which is strictly [-1 < x < 1]
        %
        % Analytical approximate wave form for asymmetric waves - Abreu et al 2010
        % http://doi.org/10.1016/j.coastaleng.2010.02.005
        % Equation 7.

        if ~isfield(S,'nonlin_deg')
            error('Please make sure S.nonlin_deg is defined to use this simulation')
        end

        if abs(S.nonlin_deg) >= 1
           warning(sprintf('Please ensure that S.nonlin_deg is [-1 < x < 1], current setting is S.nonlin_deg = %d',S.nonlin_deg))
        end

        if ~isfield(S,'nonlin_phi')
            S.nonlin_phi = -pi/2;
        end

        % create nonlinear low frequency signal
        factor = sqrt(1-S.nonlin_deg.^2);
        num = sin(2*pi*S.modulating_freq*obj.time_vect) + ( S.nonlin_deg*sin(S.nonlin_phi) / (1+sqrt(1-S.nonlin_deg^2)) );
        denom = 1- S.nonlin_deg*cos(2*pi*S.modulating_freq.*obj.time_vect + S.nonlin_phi);
        obj.modulating_ts = S.modulating_amp.*factor.*( num ./ denom );

        % create simple high frequency signal
        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        am = sin(2*pi*S.modulating_freq*obj.time_vect + S.phase_lag) + 1;
        am(am < 0) = 0;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'Abreu2010am')
        % An analytical formulation of a nonlinear near-bed wave
        %
        % The degree of nonlinearity is specified by S.nonlin_deg which is strictly [-1 < x < 1]
        %
        % The direction of skewness is determined by nonlin_phi, pi/2 is a
        % central.
        %
        % Analytical approximate wave form for asymmetric waves - Abreu et al 2010
        % http://doi.org/10.1016/j.coastaleng.2010.02.005
        % Equation 7.

        if ~isfield(S,'nonlin_deg')
            error('Please make sure S.nonlin_deg is defined to use this simulation')
        end

        if abs(S.nonlin_deg) >= 1
           warning(sprintf('Please ensure that S.nonlin_deg is [-1 < x < 1], current setting is S.nonlin_deg = %d',S.nonlin_deg))
        end

        if ~isfield(S,'nonlin_phi');
            S.nonlin_phi = -pi/2;
        end

        obj.modulating_ts = S.modulating_amp*sin(2*pi*S.modulating_freq*obj.time_vect);

        % create simple high frequency signal
        obj.modulated_ts = S.modulated_amp*sin(2*pi*S.modulated_freq*obj.time_vect);

        state_switching = ( 1 + square(sin(2*pi*S.switching_freq*obj.time_vect)) ) / 2;

        % Create modulation
        % create nonlinear amplitude modulation signal
        factor = sqrt(1-S.nonlin_deg.^2);
        num = sin(2*pi*S.modulating_freq*obj.time_vect) + ( S.nonlin_deg*sin(S.nonlin_phi) / (1+sqrt(1-S.nonlin_deg^2)) );
        denom = 1- S.nonlin_deg*cos(2*pi*S.modulating_freq.*obj.time_vect + S.nonlin_phi);
        am = -factor.*( num ./ denom );

        am = am + abs(min(am)).*S.modulated_amp;
        am(am < 0) = 0;
        obj.am = am;
        am = ( am ./ max(am) ) .* state_switching;
        am(state_switching == 0) = 1;

        % modulate signal
        obj.modulated_ts = (obj.modulated_ts.*am);

elseif strcmp(S.method,'modal')
        % Generate a phase-amplitude-coupling signal between an autoregressive
        % oscillator and a high frequency signal. The low frequency signal has
        % a bandwidth introducing more realistic features such as envelope
        % dynamics.

        % Define AR oscillator
        r = .995;
        wr = 2*pi*S.modulating_freq/S.sample_rate;
        a1 = [1 -2*r*cos(wr) (r^2)];

        time_vect = 0:1/S.sample_rate:S.seconds;
        obj.modulating_ts = randn(1,length(time_vect));

        % Generate signal from model and parameters
        for idx = 4:length(time_vect)
            for ilag = 2:3
                obj.modulating_ts(idx) = obj.modulating_ts(idx) - squeeze(a1(ilag))*obj.modulating_ts(idx-ilag+1);
            end
        end

        % Normalise
        obj.modulating_ts = ( obj.modulating_ts - mean(obj.modulating_ts) ) / std(obj.modulating_ts);

        % Generate modulated signal
        obj.modulated_ts = S.modulated_amp.*sin(2*pi*S.modulated_freq*obj.time_vect);

        state_switching = ( square(sin(2*pi*S.switching_freq*obj.time_vect))+1 )/ 2;

        % Modulate high frequency amplitude by magnitude of low frequency
        t = obj.modulating_ts + 2;
        t(t<0) = 0;
        t = t./max(t);
        obj.modulated_ts = (obj.modulated_ts.*t.*state_switching);
        
else
    warning('S.method not recognised');
end


%% Add noise scaled to the modulating time-series and return
obj.noise = scalesignal(randn(size(obj.modulated_ts,1),size(obj.modulated_ts,2)),...
                             S.noise_level,...
                             obj.modulated_ts+zscore(obj.modulating_ts));

obj.signal = obj.modulated_ts + zscore(obj.modulating_ts) + obj.noise;

obj.state_switching = state_switching;

end


function [signal,resig_db,accuracy] = scalesignal(signal, diff_db, ref_signal);

    % Get original signal power
    sig_db = 20 * log10( sqrt ( sum ( power(signal,2) ) ) );

    % Calculate actual dB scaling if passed a reference signal
    if nargin == 3
        ref_db = 20 * log10( sqrt ( sum ( power(ref_signal,2) ) ) );
        diff_db = ref_db - sig_db + diff_db;
    end

    % Calculate scaling factor
    scale_factor = 10 ^ ( diff_db / 20);

    % Scale signal
    signal = scale_factor * signal;

    % Calculate scaled signal db
    resig_db = 20 * log10( sqrt ( sum ( power(signal,2) ) ) );

    % Estimate accuracy
    accuracy = (sig_db - diff_db) / resig_db;
end

