%% Time-Frequency analysis with Empirical Mode Decomposition
%
% This tutorial introduces two concepts:
% 1) High-resolution time-frequency analysis using Instantaneous Frequency
% estimation
% 2) Separation of a complicated signal into simpler oscillatory
% components.
%
% Taken together, these steps allow us to compute time-frequency plots at a
% much higher resolution than available with most standard methods.
%
% Here we will analyse at some simulations and some example data using the
% fieldtrip toolbox

%% Getting started
%
% To run this tutorial you will need to  be running matlab R2018a or more
% recent with the signal processing toolbox.
% The  Empirical Mode Decomposition functions used here are not available
% in older versions.
%
%----------
%
% Next, you will need a recent version of Fieldtrip on the  matlab path.
% Fieldtrip and a host of tutorials are available here:
%
% http://www.fieldtriptoolbox.org
%
% Fieldtrip is used to compute a Hanning-Taper time-frequency transform for
% comparison with the Empirical Mode Decomposition.
%
%----------
%
% Finally, the simulate.m function provided here is an abridged
% version of the simulation functions in the CFC toolbox. For more details see:
%
% https://github.com/AJQuinn/cfc

%% Getting Started
%
% Add relevant toolboxes to path

addpath /path/to/fieldtrip % Update to match your folder structure
addpath /path/to/emd_hht_tutorial % Update to match your folder structure

%% Sanity Checking
%
% If you can run this cell with no warnings you should be good to go!
%
% If you get warnings, they should help identify what you are missing.

% Is matlab recent enough?
if verLessThan('matlab','9.4')
    warning('Matlab is older than R2018a - the tutorial might not run...');
end

% Do we have fieldtrip?
if isempty( which('ft_freqanalysis') )
    warning('Fieldtrip function is missing!');
end

% Do we have the signal processing toolbox and EMD functions?
if ~license( 'test', 'Signal_Toolbox' )
    warning('Signal Processing Toolbox not found - the tutorial might not run...');
end

if isempty( which('emd') )
    warning('emd function is missing!');
end
if isempty( which('hht') )
    warning('hht function is missing!');
    found = false;
end


%% Some simulations
%
% Here we will generate some simulated time-series to look at in the
% frequency domain.
%
% simple_signal contains a simple sinusoidal 10Hz modulating signal with a
% 53Hz amplitude modulated signal
%
% nonlin_signal is similar to simple_signal but the 10Hz oscillation has a
% non-sinusoidal waveform shape
%
% nonstat_signal is similar to simple_signal but the 10Hz oscillation has
% quickly changing dynamics

% Define a simple signal
S = [];
S.sample_rate = 1024;
S.seconds = 1;
S.method = 'basic';
S.modulating_freq = 10;
S.modulated_freq = 53;
simple_signal = simulate(S);

% Define a signal with a non-sinusoidal low frequency component
S.method = 'Abreu2010';
S.nonlin_deg = .5;
nonlin_signal = simulate(S);

% Define a signal with a non-stationary low frequency component
% this signal will be different each time unless a random seed is set.
S.method = 'modal';
S.nonlin_deg = .5;
nonstat_signal = simulate(S);

% Concatenate the low and high frequency components of each signal
imf1 = [simple_signal.modulated_ts; simple_signal.modulating_ts]';
imf2 = [nonlin_signal.modulated_ts; nonlin_signal.modulating_ts]';
imf3 = [nonstat_signal.modulated_ts; nonstat_signal.modulating_ts]';

%% Hilbert-Huang Transform
%
% Here we will compute and plot the Hilbert-Huang Transform for each signal
% using the plot_hht function. This calls hht internally and creates a
% simple visualisation. Note, these plots show amplitude rather than
% power (power = amplitude.^2)
%
% The HHT plot is a sparse distribution of the instantaneous amplitude and
% frequency of each component in the signal. The instantaneous amplitude (IA)
% is computed from the abs of the Hilbert Transform. The instantaneous
% frequency (IF) is computed from the derivative of the instataneous phase.

% We compute a value for each signal at every time-point. Next we build a
% sparse matrix with an entry in every column for each signal component
% (here we have 2 components). We define a set of frequency bins and enter
% each IA value into the bin corresponding to its IF value.
%
%
% The first signal shows constant amplitude at 10Hz and modulated power at 53Hz
%
% The second signal shows rapid frequency variation in 10Hz component. The
% HHT has sufficient resolution to resolve frequency variation within
% single cycles. The modulation in frequency within cycles relates to which
% parts of the cycle are progressing faster or slower. This signal has a
% narrow/fast peak with a frequency around 12Hz and a wide/slow trough with
% a frequency around 8Hz.
%
% The third signal shows rapid changes in time. Here both the frequency and
% amplitude of both signals shows rapid changes. Note that the frequency
% estimates are very noisy when the signal is not a relatively simple
% mono-component oscillation.

figure('Position',[100 100 1024 512])
subplot(3,3,[1,4])
plot_hht(imf1,S.sample_rate,gca)
title('Standard PAC')
subplot(3,3,7)
plot(simple_signal.time_vect,imf1)
xlabel('Time (seconds)')
grid on

subplot(3,3,[2,5])
plot_hht(imf2,S.sample_rate,gca)
title('Non-sinusoidal PAC')
subplot(3,3,8)
plot(nonlin_signal.time_vect,imf2)
xlabel('Time (seconds)')
grid on

subplot(3,3,[3,6])
plot_hht(imf3,S.sample_rate,gca)
title('Non-stationary PAC')
subplot(3,3,9)
plot(nonstat_signal.time_vect,imf3)
xlabel('Time (seconds)')
grid on

%% Empirical Mode Decomposition
%
% The Hilbert-Huang Transform allows us to compute a high-resolution
% time-frequency transform BUT requires a set of well-behaved
% relatively-monocomponent set of signals to operate on. The Empirical Mode
% Decomposition is able to split a complicated signal into such a set of
% components using the Sift algorithm
%
% Here we look at an example from the matlab emd documentation. A composite
% signal is defined from a number of different signals which change over time.
%
% This signal is then split into a set of Intrinsic Mode Functions (IMFs)
% using the emd. Note that the IMFs are relatively clear oscillatory
% signals but they retain the rapid dynamics in the signal. These IMFs are
% appropriate inputs into the Hilbert-Huang Transform.
%
% The first IMF contains the fastest dynamics in the signal. Each
% successive IMF isolates slower and slower dynamics until the final IMF
% contains only a non-oscillatory trend.
%
% The top of the figure output contains the mixed signal and the second
% plot contains first 4 Intrinsic Mode Functions. The bottom figure
% contains the HHT.

% Define our signal
fs = S.sample_rate;
t = (0:1/S.sample_rate:4)';
x1 = sin(2*pi*50*t) + sin(2*pi*200*t) + 0.1*randn(length(t),1);
x2 = sin(2*pi*50/2*t) + sin(2*pi*200/2*t) + sin(2*pi*250*t) + 0.1*randn(length(t),1);
x = [x1;x2];
time_vect = linspace(0,8,size(x,1));

% Compute Intrinsic Mode Functions with Empirical Mode Decomposition
imf = emd(x);

% Plot the simulated signal and the IMFs
figure
subplot(7,1,1)
plot(time_vect,x)
title('Time Series')
xlim([3 5])
subplot(7,1,[2 3 4])
n_imfs = 4;
plot(time_vect,imf(:,1:n_imfs) - 2*repmat(1:n_imfs,size(imf,1),1))
yticks(-2 * (7:-1:1))
yticklabels({'IMF7','IMF6','IMF5','IMF4','IMF3','IMF2','IMF1'})
title('Intrinsic Mode Functions')
xlim([3 5])
ax = subplot(7,1,[5 6 7]);
plot_hht( imf, S.sample_rate, ax, 350);
title('Hilbert-Huang Transform')
xlim([3 5])

%% Fieldtrip Example Data
%
% Here we load in some fieldtrip example data from this tutorial:
%
% http://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis
%
% Please download this file and make sure it is in the current folder
%
% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/timefrequencyanalysis/dataFIC.mat
%
% We take a signal trial from one channel from this data and apply the EMD
% to extract IMFs.
%
% Again, note that the IMFs are relatively narrow-band oscillations with
% the fastest dynamics in the first IMF. They are highly dynamic within an
% IMF as the EMD does not supress non-linearity or non-stationarity in a
% signal.

% Load data
load dataFIC.mat
channel = 9; % MLC24

% Compute IMFs with EMD.
[imf, residual, info] = emd(dataFIC.trial{1}(channel,:));

% Summary figure
figure
subplot(4,1,1)
plot(dataFIC.time{1},dataFIC.trial{1}(channel,:))
title('Time Series')
grid on

subplot(4,1,[2 3 4])
plot(dataFIC.time{1},imf - 1e-12*repmat(1:5,size(imf,1),1))
yticks(-1e-12 * (5:-1:1))
yticklabels({'IMF5','IMF4','IMF3','IMF2','IMF1'})
grid on
title('Intrinsic Mode Functions')


%% Compute a Hanning Taper Time-Frequency Analysis
%
% This is a standard Hanning-Taper analysis based on this tutorial
%
% http://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis/
%
% Please check this page for a more detailed discussion of Hanning-Taper TF
% estimation.
%
% We will use this as a comparison for the EMD & HHT approach.
load dataFIC.mat

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:1:40;                         % analysis 2 to 30 Hz in steps of 2 Hz
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -1:0.05:2;                      % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg.pad          = 'nextpow2';                     % Faster FFT computation
TFRhann = ft_freqanalysis(cfg, dataFIC);

%% Compute a EMD + HHT Time Frequency analysis
%
% Here we compute a EMD+HHT analysis on the example data. This cell
% performs several operations
%
% 1) Preallocate arrays
% 2) Define channel and smoothing filter
% 3) Loop through trials
% 3.1) Compute IMFs from trial using EMD
% 3.2) Compute instantaneous frequency and energy from IMFs for this trial
% 3.3) Loop through IMFs
% 3.3.1) Create sparse HHT array
% 3.3.2) Smooth sparse HHT array


% 1) Preallocate arrays
%
% we have:
% * 77 trials
% * 900 time-points
% * 10 IMFs
% * 100 Frequencies

erf = zeros(77,900);          % Evoked Field
imfs = zeros(77,900,10);       % Intrinsic Mode Functions
hhts = zeros(77,100,900,10);   % Hilbert-Huang Transform

% 2) Define channel and smoothing filter
channel = 9; % MLC24
h = fspecial('gauss', 12,2);

% 3) Loop through trials
for ii = 1:77

    % 3.1) Compute IMFs from trial using EMD
    erf(ii,:) = dataFIC.trial{ii}(channel,:);
    tmp = emd(erf(ii,:),'Interpolation','pchip',...
                        'SiftRelativeTolerance',.02,...
                        'Display',0);
    imfs(ii,:,1:size(tmp,2)) = tmp;

    % 3.2) Compute instantaneous frequency and energy from IMFs for this trial
    [P1,F1,T1,insf,inse] = hht(tmp,dataFIC.fsample,'FrequencyLimits',[0,75],'FrequencyResolution',.5);

    % 3.3) Loop through IMFs
    for jj = 1:size(tmp,2)
        % 3.3.1) Create sparse HHT array
        freq_bins = discretize(insf(:,jj),0:100);
        good_bins = find(isfinite(freq_bins));
        s = sparse(good_bins,freq_bins(good_bins),inse(good_bins,jj));
        % 3.3.2) Smooth sparse HHT array
        s = filter2(h,full(s)');
        hhts(ii,1:size(s,1),1:size(s,2),jj) = s;
    end

end

%% Method comparison figure
%
% Here create a figure containing:
% * The Evoked Field
% * The trial-averaged Hanning-Taper time-frequency representation
% * The trial-averaged Hilbert-Huang time-frequency representation
%
% Note that both transforms pick up on similar features, though the EMD+HHT
% approach has less blurring across time and frequency. For instance, the
% increase in beta power around 1 second is quite gradual in the
% Hanning-Taper but almost instantanous for the HHT.

figure
subplot(511)
plot(dataFIC.time{1},mean(erf),'k')
title('Evoked Field')
grid on;

subplot(5,1,[2 3])
contourf(TFRhann.time,TFRhann.freq, squeeze(TFRhann.powspctrm(channel,:,:)),32,'linestyle','none')
title('Hanning Tapers')
ylabel('Frequency (Hz)');
grid on

subplot(5,1,[4 5])
toplot = squeeze(sum(sum(hhts(:,:,:,1:5),1),4));
contourf(dataFIC.time{1},1:100,toplot,linspace(0,10e-25,36),'linestyle','none');
title('HHT')
grid on
xlabel('Time (seconds')
ylabel('Frequency (Hz)');
ylim([0 40])



%% Method comparison figure
%
% Here create a figure containing:
% * The Evoked Field
% * The baseline-corrected trial-averaged Hanning-Taper time-frequency representation
% * The baseline-corrected Hilbert-Huang time-frequency representation
% using a baseline of -.5 to -.1 seconds
%
% The picture is similar to above, again note the increased sharpness of
% the beta onset around 1 second

figure
subplot(511)
plot(dataFIC.time{1},mean(erf),'k')
grid on;

baseline_inds = TFRhann.time > -.5 & TFRhann.time < -.1;
ax = subplot(5,1,[2 3]);
toplot = squeeze(TFRhann.powspctrm(channel,:,:));
toplot = toplot - mean(toplot(:,baseline_inds),2);
contourf(TFRhann.time,TFRhann.freq, toplot,linspace(-3e-27,3e-27,36),'linestyle','none')
grid on;
grid(gca,'minor')
colorbar('NorthOutside')
caxis([-3e-27,3e-27])
set(ax,'XMinorTick','on');
ax.XAxis.MinorTickValues = TFRhann.time;
title('Hanning Tapers')
xlabel('Time (seconds')
ylabel('Frequency (Hz)')

baseline_inds = dataFIC.time{1} > -.5 & dataFIC.time{1} < -.1;
subplot(5,1,[4 5])
toplot = squeeze(mean(sum(hhts(:,:,:,1:7),4),1));
toplot = toplot - mean(toplot(:,baseline_inds),2);
contourf(dataFIC.time{1},1:100,toplot,linspace(-3e-27,3e-27,36),'linestyle','none');
title('HHT')
grid on;caxis([-3e-27,3e-27])
ylim([0 40])
colorbar('NorthOutside')
xlabel('Time (seconds')
ylabel('Frequency (Hz)')

%% Advanced points
%
% There are two disadvantages of the EMD+HHT you may have noticed. These
% can both be greatly in other implementations of the EMD. We have used the
% matlab version here for simplicity but strongly recommmend considering these
% factors if using the EMD in your own future applications.
%
% Firstly,is that the instantanous frequency estimation can be quite noisy.
% This arises as it is computed from the differential of the instantaneous
% phase. Differentiation is challenging in a noisy signal and can lead to
% increased noise levels in the differentiated signal. This can be greatly
% improved by using smoothing filters on the phase prior to computing
% frequency.
%
% Secondly, mode-mixing can be a major problem in the EMD. this occurs when
% a signal component isn't cleanly isolated into a single IMF. This is
% particularly common for transient signals which come on and off quickly
% in a noisy signal. This problem can be reduced by using the Ensemble EMD
% or the Masked EMD to ensure a cleaner separation of signals. There is no
% clear way to guarantee that signals will be separated a priori for a new signal
%
% We strongly recommend that each EMD project should explore the masked and
% ensemble EMD options to make sure that decompositions are of a high
% quality

%% Further Reading
%
% The original EMD paper - 92 pages and 12,000 citations...
%
% Huang, N. E., Shen, Z., Long, S. R., Wu, M. C., Shih, H. H., Zheng,
%        Q., ? Liu, H. H. (1998). The empirical mode decomposition and the Hilbert
%        spectrum for nonlinear and non-stationary time series analysis. Proceedings
%        of the Royal Society of London. Series A: Mathematical, Physical and
%        Engineering Sciences, 454(1971), 903?995.
%        https://doi.org/10.1098/rspa.1998.0193
%
%
% Some more detailed discussion on instantaneous frequency and how to
% estimate it accurately
%
% Huang, N. E., Wu, Z., Long, S. R., Arnold, K. C., Chen, X., & Blank,
%        K. (2009). On Instantaneous Frequency. Advances in Adaptive Data Analysis,
%        1(2), 177?229. https://doi.org/10.1142/s1793536909000096
%
%
% A discussion on mode-mixing and how masked EMD can reduce the problem
%
% Ryan Deering, & James F. Kaiser. (2005). The Use of a Masking Signal
%        to Improve Empirical Mode Decomposition. In Proceedings. (ICASSP ?05). IEEE
%        International Conference on Acoustics, Speech, and Signal Processing, 2005.
%        IEEE. https://doi.org/10.1109/icassp.2005.1416051
%
% A frank discussion of some problems with the EMD and how you might solve
% them. Very practical paper focused on implementation details
%
% Rato, R. T., Ortigueira, M. D., & Batista, A. G. (2008). On the HHT,
%        its problems, and some solutions. Mechanical Systems and Signal
%        Processing, 22(6), 1374?1394. https://doi.org/10.1016/j.ymssp.2007.11.028
