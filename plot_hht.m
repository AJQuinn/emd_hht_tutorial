function plot_hht( imf, sample_rate, ax, max_freq )
%% function plot_hht( imf, sample_rate, ax, max_freq )
%
% A function to compute and plot Hilbert-Huang Transforms
% based on signalwavelet.internal.guis.plot.hhtPlot

%%
% Housekeeping variables
if nargin < 4  || isempty(max_freq)
    max_freq = 75;
end

if nargin < 3 || isempty(ax)
    figure
    subplot(111);
    ax = gca;
    hold on
end

seconds = size(imf,1)/sample_rate;
time = linspace(0,seconds,seconds*sample_rate);

%%
% compute instantaneous stats
[~,~,~,insf,inse] = hht(imf,sample_rate,'FrequencyLimits',[1,max_freq],'FrequencyResolution',.5);

%%
% Plot HHT
for ii = 1:size(insf,2)
    patch( ax,[nan;time';nan],[nan;insf(:,ii);nan],[nan;sqrt(inse(:,ii));nan],...
                        'EdgeColor','interp','EdgeAlpha','interp',...
                        'FaceColor', 'none', 'FaceVertexAlphaData',[nan;sqrt(inse(:,ii));nan],...
                        'LineWidth', 2, 'FaceAlpha', 'interp');
end

ylim(ax,[0 max_freq]);
xlim(ax,[0,seconds]);
xlabel('Time (seconds)')
ylabel('Frequency (Hz)')
colorbar('North')
grid on
end