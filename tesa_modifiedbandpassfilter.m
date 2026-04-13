% tesa_modifiedbandpassfilter() - do a variant of band-pass filtering, 
%                                 specialized to avoid propagating large-amplitude 
%                                 low-latency stimulation artifact to surrounding times
%                                 Auto-regressive prediction is used to pad
%                                 the time periods around where the TMS
%                                 artifact has been removed. Based on
%                                 function from the AARATEP pipeline by
%                                 Chris Cline. 
%
% Usage:
%   >>  EEG = tesa_modifiedbandpassfilter( EEG );
%   >>  EEG = tesa_modifiedbandpassfilter( EEG, 'key1', value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs:
%   'lowCutoff',num   - num defines the low band cut-off frequency for a high-pass filter. 
%                   num is a number in Hz. 
%                   default = 1
% 
%   'highCutoff',num   - num defines the high band cut-off frequency for a low-pass filter. 
%                   num is a number in Hz. 
%                   default = []
% 
%   'filterMethod','str'- 'str' is a filtering method. 'butterworth' | 'eegfiltnew' 
%                   default = 'butterworth'
%
%   'artifactTimespan', [vec] - vec is a vector defining the start and end
%                   time of TMS-pulse artifact removal in seconds. If tesa_removedata
%                   has been run, the most recent cut values will be used
%                   and multiplied by 3. 
%                   default = [] 
%
%   'doPiecewise','logical' - true|false. Do autoregressive prediction forward from negative times and backward from positive times,
%                   and apply filter in piecewise fashion to avoid large post-pulse artifact from "leaking"
%                   into pre-pulse time periods
%                   default = true
%  
%   'pieceWiseTimeToExtend',num - num defines the length of the time window for the 
%                   autoregressive prediction in seconds.
%                   default = 0.5
%
%   'prePostExtrapolationDurations',[vec] - Regardless of center artifact timespan, you may want to extrapolate the signal on
%                   either end of available data to reduce filter edge effects.
%                   Specify positive prePostExtrapolationDurations values to do so.
%                   default = []
%
%   'interpolationArgs',{cell} - optional, non-default interpolation args
%                   default = {}
%
%   'doDebug',logical - true|false. Do additional debug plots. 
%                   default = false
%
%   'filtOrder', num - defines the filter order (for Butterworth filters)
%                   default = 4
%
%   'filtType', 'str' - defines the filter type. 'auto' | 'bandpass'| 'bandstop' | 'highpass' | 'lowpass'
%                   auto = selects filter type based on lowCutoff and
%                   highCutoff inputs. If both are defined then a band-pass
%                   is used. If only one is defined then either high-pass
%                   or low-pass filter is used accordingly.
%                   default = 'auto'
%
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples
%   EEG = tesa_modifiedbandpassfilter( EEG ); % default use
%   EEG = tesa_modifiedbandpassfilter( EEG, 'lowCutoff', 2, 'pieceWiseTimeToExtend', 0.9 ); %user defined input
%
% See also:
%   tesa_filtbutter

% This script was adapted by Nigel Rogasch for the TESA toolbox. Original
% code is available from:
% https://github.com/chriscline/AARATEPPipeline

% MIT License
% 
% Copyright (c) 2021 Chris Cline
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function EEG = tesa_modifiedbandpassfilter(EEG,varargin)

if nargin < 1
	error('Not enough input arguments.');
end

p = inputParser();
% p.addRequired('EEG', @isstruct);
p.addParameter('lowCutoff', 1, @(x) isempty(x) || isscalar(x));
p.addParameter('highCutoff',[], @(x) isempty(x) || isscalar(x));
p.addParameter('filterMethod', 'butterworth', @ischar);
p.addParameter('artifactTimespan', [], @c_isSpan); % this can be wider than typical timespan to be sure to remove large artifact signals 
												   %  (only used for interpolating a temporary signal for filtering, not the output)
p.addParameter('doPiecewise', true, @islogical);
p.addParameter('piecewiseTimeToExtend', 0.5, @isscalar);
p.addParameter('prePostExtrapolationDurations', [0 0], @(x) isvector(x) && all(x >= 0) && length(x)==2);
													% regardless of center artifact timespan, you may want to extrapolate the signal on
													% either end of available data to reduce filter edge effects.
													% Specify positive prePostExtrapolationDurations values to do so.
p.addParameter('interpolationArgs', {}, @iscell); % optional, non-default interpolation args
p.addParameter('doDebug', false, @islogical);
p.addParameter('filtOrder', 4, @(x) isempty(x) || isscalar(x));
p.addParameter('filtType', 'auto', @ischar); 
p.parse(varargin{:});
s = p.Results;
s.EEG = EEG;

assert(any(arrayfun(@(x) ~isempty(x) && x>0, [s.lowCutoff, s.highCutoff])), 'Must specify at least one of lowCutoff and highCutoff');

% Check if tmscut exists and if it does, use 3 times these values
if isempty(s.artifactTimespan)
    if isfield(EEG,'tmscut')
        lastcut = size(EEG.tmscut,2);
        artifactTimespan = EEG.tmscut(lastcut).cutTimesTMS * 0.001;
        s.artifactTimespan = artifactTimespan*3;
    end
end

assert(~isempty(s.artifactTimespan), 'Must specify artifact timespan');

if s.doDebug
	iTr = 1;
	iCh = 1;
end

if s.doPiecewise
	% do autoregressive prediction forward from negative times and backward from positive times,
	%  and apply filter in piecewise fashion to avoid large post-pulse artifact from "leaking"
	%  into pre-pulse time periods
	interpArgs = c_cellToStruct(s.interpolationArgs);
	if ~isfield(interpArgs, 'method')
		interpArgs.method = 'ARExtrapolation';
	else
		assert(isequal(interpArgs.method, 'ARExtrapolation'), 'Only one interp method supported for piecewise modification');
	end
	
	if ~c_isFieldAndNonEmpty(interpArgs, 'prePostFitDurations')
		% use a longer pre duration since signals are expected to be more stationary there
		prePostFitDurations = [100 30]*1e-3;
	else
		prePostFitDurations = interpArgs.prePostFitDurations;
	end
	
	assert(EEG.xmin - s.prePostExtrapolationDurations(1) <= s.artifactTimespan(2) - s.piecewiseTimeToExtend);
	assert(EEG.xmax + s.prePostExtrapolationDurations(2) >= s.artifactTimespan(1) + s.piecewiseTimeToExtend);
	
	interpArgs = c_structToCell(interpArgs);
	EEG_pre = c_EEG_ReplaceEpochTimeSegment(EEG,...
		'timespanToReplace', [s.artifactTimespan(1), s.artifactTimespan(1) + s.piecewiseTimeToExtend],...
		interpArgs{:},...
		'prePostFitDurations', [prePostFitDurations(1) 0]);
	
	if s.prePostExtrapolationDurations(1) > 0
		EEG_pre = c_EEG_ReplaceEpochTimeSegment(EEG_pre,...
			'timespanToReplace', [EEG.xmin - s.prePostExtrapolationDurations(1), EEG.xmin],...
			interpArgs{:},...
			'prePostFitDurations', [0, s.artifactTimespan(1) - EEG.xmin]);
		if s.prePostExtrapolationDurations(2) > 0
			EEG_pre = c_EEG_grow(EEG_pre,...
				'timespan', [EEG_pre.xmin, EEG_pre.xmax + s.prePostExtrapolationDurations(2)],...
				'padWith', NaN);
		end
	end

	EEG_post = c_EEG_ReplaceEpochTimeSegment(EEG,...
		'timespanToReplace', [s.artifactTimespan(2) - s.piecewiseTimeToExtend, s.artifactTimespan(2)],...
		interpArgs{:},...
		'prePostFitDurations', [0 prePostFitDurations(2)]);

	if s.prePostExtrapolationDurations(2) > 0
		EEG_post = c_EEG_ReplaceEpochTimeSegment(EEG_post,...
			'timespanToReplace', [EEG.xmax, EEG.xmax + s.prePostExtrapolationDurations(2)],...
			interpArgs{:},...
			'prePostFitDurations', [EEG.xmax - s.artifactTimespan(2), 0]);
		if s.prePostExtrapolationDurations(1) > 0
			EEG_post = c_EEG_grow(EEG_post,...
				'timespan', [EEG_post.xmin - s.prePostExtrapolationDurations(1), EEG_post.xmax],...
				'padWith', NaN);
		end
	end

	if any(s.prePostExtrapolationDurations > 0)
		EEG = c_EEG_grow(EEG, 'timespan', [EEG.xmin EEG.xmax] + [-1 1] .* s.prePostExtrapolationDurations, 'padWith', NaN);
	end

	% blend pre-filter results
	preIndices = EEG.times < s.artifactTimespan(1) * 1e3;
	postIndices = EEG.times > s.artifactTimespan(2) * 1e3;
	blendIndices = ~preIndices & ~postIndices;
	if true
		% sigmoidish 
		fn = @(x, k)  1 - 1./(1+(1./x - 1).^-k);
		pre_weights = fn(linspace(1, 0, sum(blendIndices)), 2);
	else
		% linear
		pre_weights = linspace(1, 0, sum(blendIndices));
	end
	post_weights = 1 - pre_weights;
	tmpEEG = EEG_pre;
	tmpEEG.data(:,postIndices, :) = EEG_post.data(:, postIndices, :);
	tmpEEG.data(:, blendIndices, :) = tmpEEG.data(:, blendIndices, :) .* pre_weights + EEG_post.data(:, blendIndices, :) .* post_weights;
	
	if s.doDebug
		
		hf = figure;
		has = gobjects(0);
		
		numSubplots = 9;
		
		for iSP=1:numSubplots
			has(end+1) = c_subplot(numSubplots+1, 1, iSP);
		end
		
		colors = [...
					0 0 0;
					0.8 0 0;
					0 0 0.8;
					0 0.6 0;
					0 0.4 0];
		
		if false
			plotIndices = EEG.times >= (s.artifactTimespan(2) - s.piecewiseTimeToExtend*0.5)*1e3 & EEG.times <= (s.artifactTimespan(1) + s.piecewiseTimeToExtend*0.5)*1e3;
		else
			plotIndices = true(size(EEG.times));
		end
		assert(sum(plotIndices)>1);
		
		% plot original
		axes(has(1));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 1];
		ylabel('Original');
		
		% plot pre-stim extrapolation
		axes(has(2));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		fitIndices = plotIndices & EEG.times >= (s.artifactTimespan(1) - prePostFitDurations(1))*1e3 & EEG.times <= s.artifactTimespan(1)*1e3;
		hp = plot(EEG.times(fitIndices), EEG.data(iCh, fitIndices, iTr));
		hp.Color = [colors(1,:) 1];
		extrapIndices = plotIndices & EEG.times >= s.artifactTimespan(1)*1e3 & EEG.times <= (s.artifactTimespan(1) + s.piecewiseTimeToExtend)*1e3;
		hp = plot(EEG.times(extrapIndices), EEG_pre.data(iCh, extrapIndices, iTr));
		hp.Color = [colors(2,:), 0.5];
		if s.prePostExtrapolationDurations(1) > 0
			extrapIndices = EEG.times < s.EEG.xmin*1e3 & EEG.times > (s.EEG.xmin - s.prePostExtrapolationDurations(1))*1e3;
			hp = plot(EEG.times(extrapIndices), EEG_pre.data(iCh, extrapIndices, iTr));
			hp.Color = [colors(2,:)*0.5, 0.5];
		end
		ylabel('pre-stim extrap');
		
		% leave a place for pre-stim extrap filtered
		
		% plot post-stim extrapolation
		axes(has(4));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		fitIndices = plotIndices & EEG.times <= (s.artifactTimespan(2) + prePostFitDurations(1))*1e3 & EEG.times >= s.artifactTimespan(2)*1e3;
		hp = plot(EEG.times(fitIndices), EEG.data(iCh, fitIndices, iTr));
		hp.Color = [colors(1,:) 1];
		extrapIndices = plotIndices & EEG.times <= s.artifactTimespan(2)*1e3 & EEG.times >= (s.artifactTimespan(2) - s.piecewiseTimeToExtend)*1e3;
		hp = plot(EEG.times(extrapIndices), EEG_post.data(iCh, extrapIndices, iTr));
		hp.Color = [colors(3,:), 0.5];
		if s.prePostExtrapolationDurations(2) > 0
			extrapIndices = EEG.times > s.EEG.xmax*1e3 & EEG.times < (s.EEG.xmax + s.prePostExtrapolationDurations(2))*1e3;
			hp = plot(EEG.times(extrapIndices), EEG_post.data(iCh, extrapIndices, iTr));
			hp.Color = [colors(3,:)*0.5, 0.5];
		end
		ylabel('post-stim extrap');
		
		% leave a place for post-stim extrap filtered
		
		% plot pre-filter blended
		axes(has(6));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		thisIndices = plotIndices & EEG.times <= s.artifactTimespan(1)*1e3;
		hp = plot(EEG.times(thisIndices), EEG.data(iCh, thisIndices, iTr));
		hp.Color = [colors(1,:) 1];
		thisIndices = plotIndices & EEG.times >= s.artifactTimespan(2)*1e3;
		hp = plot(EEG.times(thisIndices), EEG.data(iCh, thisIndices, iTr));
		hp.Color = [colors(1,:) 1];
		hp = plot(EEG.times(blendIndices), tmpEEG.data(iCh, blendIndices, iTr));
		hp.Color = [colors(4,:) 0.5];
		ylabel('pre-filt blended');
	end
	
	% apply filter
	if any(s.prePostExtrapolationDurations > 0)
		% extra wrappers to remove NaN times before filtering, and replace them after
		% not very efficient...
		EEG_pre = pop_select(EEG_pre, 'time', [EEG_pre.xmin, s.EEG.xmax]);
		EEG_pre = applyFilter(EEG_pre, s);
		EEG_pre = c_EEG_grow(EEG_pre, 'timespan', [EEG_post.xmin, EEG_post.xmax], 'padWith', NaN);

		EEG_post = pop_select(EEG_post, 'time', [s.EEG.xmin, EEG_post.xmax]);
		EEG_post = applyFilter(EEG_post, s);
		EEG_post = c_EEG_grow(EEG_post, 'timespan', [EEG_pre.xmin, EEG_pre.xmax], 'padWith', NaN);
	else
		EEG_pre = applyFilter(EEG_pre, s);
		EEG_post = applyFilter(EEG_post, s);
	end
	
	if s.doDebug
		% plot pre-stim filtered
		axes(has(3));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		thisIndices = plotIndices & EEG.times <= (s.artifactTimespan(1) + s.piecewiseTimeToExtend)*1e3;
		hp = plot(EEG.times(thisIndices), EEG_pre.data(iCh, thisIndices, iTr));
		hp.Color = [colors(2,:) 0.5];
		ylabel('pre-stim filt');
		
		% plot post-stim filtered
		axes(has(5));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		thisIndices = plotIndices & EEG.times >= (s.artifactTimespan(2) - s.piecewiseTimeToExtend)*1e3;
		hp = plot(EEG.times(thisIndices), EEG_post.data(iCh, thisIndices, iTr));
		hp.Color = [colors(3,:) 0.5];
		ylabel('post-stim filt');
	end
	
	% blend post-filter results
	tmpEEG2 = EEG_pre;
	EEG_pre = [];
	tmpEEG2.data(:,postIndices, :) = EEG_post.data(:, postIndices, :);
	tmpEEG2.data(:, blendIndices, :) = tmpEEG2.data(:, blendIndices, :) .* pre_weights + EEG_post.data(:, blendIndices, :) .* post_weights;
	EEG_post = [];
	
	if s.doDebug
		% plot blended filtered
		axes(has(7));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		hp = plot(EEG.times(plotIndices), tmpEEG2.data(iCh, plotIndices, iTr));
		hp.Color = [colors(4,:) 0.5];
		ylabel('blended filt');
	end
	
	tmpEEG.data = tmpEEG.data - tmpEEG2.data;
	tmpEEG2 = [];
	
	if s.doDebug
		% plot residual
		axes(has(8));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		hp = plot(EEG.times(plotIndices), tmpEEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(4,:) 0.5];
		ylabel('residual');
		
		% prepare to plot final
		axes(has(9));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
	end
	
	EEG.data = EEG.data - tmpEEG.data;
	
	if s.doDebug
		% plot final
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(5,:) 1];
		
		c_plot_setEqualAxes(has);
		xlim(extrema(EEG.times(plotIndices)));
		xlabel('Time (ms)');
		ylabel('final');
		
		for iSP = 1:numSubplots
			if iSP < numSubplots
				has(iSP).XTickLabel = {};
			end
			has(iSP).YTickLabel = {};
			set(has(iSP).Children, 'LineWidth', 1.5);
			has(iSP).Position([2 4]) = [1/(numSubplots+1)*(numSubplots + 1 - iSP) 1/(numSubplots+1)*0.9];
		end
		
		hf.Position = [2 50 440 1300];
	end

	if any(s.prePostExtrapolationDurations > 0)
		% cut data back down to original timespan
		EEG = pop_select(EEG, 'time', [s.EEG.xmin, s.EEG.xmax]);
	end
	
	
else

	if any(s.prePostExtrapolationDurations > 0)
		error('Not implemented');
	end
	
	if s.doDebug
		
		hf = figure;
		has = gobjects(0);
		
		numSubplots = 5;
		
		for iSP=1:numSubplots
			has(end+1) = c_subplot(numSubplots+1, 1, iSP);
		end
		
		colors = [...
					0 0 0;
					0 0.6 0;
					0 0.4 0];
	end
	
	interpArgs = c_cellToStruct(s.interpolationArgs);
	if ~isfield(interpArgs, 'method')
		interpArgs.method = 'ARExtrapolation';
	end
	
	if isequal(interpArgs.method, 'ARExtrapolation')
		if ~c_isFieldAndNonEmpty(interpArgs, 'prePostFitDurations')
			% use a longer pre duration since signals are expected to be more stationary there
			interpArgs.prePostFitDurations = [100 30]*1e-3;
		end
	end
	interpArgs = c_structToCell(interpArgs);

	tmpEEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
		'timespanToReplace', s.artifactTimespan,...
		interpArgs{:});

	tmpEEG2 = applyFilter(tmpEEG, s);
	
	if s.doDebug
		if false
			plotIndices = EEG.times >= (s.artifactTimespan(2) - 0.25)*1e3 & EEG.times <= (s.artifactTimespan(1) + 0.25)*1e3;
		else
			plotIndices = true(size(EEG.times));
		end
		assert(sum(plotIndices)>1);
		
		% plot original
		axes(has(1));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 1];
		
		% plot interp
		axes(has(2));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		thisIndices = plotIndices & EEG.times <= s.artifactTimespan(1)*1e3;
		hp = plot(EEG.times(thisIndices), EEG.data(iCh, thisIndices, iTr));
		hp.Color = [colors(1,:) 1];
		thisIndices = plotIndices & EEG.times >= s.artifactTimespan(2)*1e3;
		hp = plot(EEG.times(thisIndices), EEG.data(iCh, thisIndices, iTr));
		hp.Color = [colors(1,:) 1];
		blendIndices = plotIndices & EEG.times >= s.artifactTimespan(1)*1e3 & EEG.times <= s.artifactTimespan(2)*1e3;
		hp = plot(EEG.times(blendIndices), tmpEEG.data(iCh, blendIndices, iTr));
		hp.Color = [colors(2,:) 0.5];
		
		% plot interp filtered
		axes(has(3));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		hp = plot(EEG.times(plotIndices), tmpEEG2.data(iCh, plotIndices, iTr));
		hp.Color = [colors(2,:) 0.5];
	end
	
	tmpEEG.data = tmpEEG.data - tmpEEG2.data;
	tmpEEG2 = [];
	
	if s.doDebug
		% plot residual
		axes(has(4));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
		hp = plot(EEG.times(plotIndices), tmpEEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(2,:) 0.5];
		
		% prepare to plot final
		axes(has(5));
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(1,:) 0.1];
		hold on;
	end
	
	EEG.data = EEG.data - tmpEEG.data;
	
	if s.doDebug
		hp = plot(EEG.times(plotIndices), EEG.data(iCh, plotIndices, iTr));
		hp.Color = [colors(3,:) 1];
		
		c_plot_setEqualAxes(has);
		xlim(extrema(EEG.times(plotIndices)));
		xlabel('Time (ms)');
		
		for iSP = 1:numSubplots
			if iSP < numSubplots
				has(iSP).XTickLabel = {};
			end
			has(iSP).YTickLabel = {};
			set(has(iSP).Children, 'LineWidth', 1.5);
			has(iSP).Position([2 4]) = [1/(numSubplots+1)*(numSubplots + 1 - iSP) 1/(numSubplots+1)*0.9];
		end
		
		hf.Position = [2 50 440 1300*6/10];
	end
end

end

%%
function EEG = applyFilter(EEG, s)
	switch(s.filterMethod)
		case 'eegfiltnew'
            if strcmp(s.filtType,'bandstop')
                EEG = pop_eegfiltnew(EEG, 'locutoff',s.lowCutoff,'hicutoff',s.highCutoff,'revfilt',1);
            else
    			EEG = pop_eegfiltnew(EEG, s.lowCutoff, s.highCutoff);
            end
		case 'butterworth'
			EEG = c_EEG_filter_butterworth(EEG, [ ...
				c_if(isempty(s.lowCutoff), 0, s.lowCutoff), ...
				c_if(isempty(s.highCutoff), 0, s.highCutoff)],...
                'order',c_if(isempty(s.filtOrder), 4, s.filtOrder),...
                'type',c_if(isempty(s.filtType), 'auto', s.filtType));
		otherwise
			error('Not implemented');
	end
end

%%
function isSpan = c_isSpan(x)
    isSpan = isvector(x) && isnumeric(x) && length(x)==2 && x(1) <= x(2);
end

%%
function s = c_cellToStruct(c,varargin)
% convert a cell array of named parameter arguments to a struct array

if nargin == 0 % example
	s = c_cellToStruct({'FirstArg',1,'SecondArg','example'});
	return
end

p = inputParser();
p.addParameter('RecursionLevel',1,@isscalar); % e.g. to convert an appropriately formatted cell array to a struct of structs, set RecursionLevel=2
p.parse(varargin{:});

if ~iscell(c)
	error('Input should be a cell array');
end

s = struct();

if isempty(c)
	return;
end

if mod(length(c),2)~=0
	error('Cell should consist of pairs of (name, value) elements');
end

for i=1:length(c)/2
	j = i*2-1;
	if p.Results.RecursionLevel < 2
		if false
			s.(c{j}) = c{j+1};
		else
			s = c_setField(s, c{j}, c{j+1});  % support field.subfield nesting
		end
	else
		s.(c{j}) = c_cellToStruct(c{j+1},'RecursionLevel',p.Results.RecursionLevel-1);
	end
end

end

%%
function res = c_isFieldAndNonEmpty(struct,field)
	res = c_isField(struct,field) && ~c_isEmptyOrEmptyStruct(c_getField(struct,field));
end

%%
function isFieldAndTrue = c_isFieldAndTrue(s,field)
	isFieldAndTrue = c_isField(s,field) && c_use(c_getField(s,field),@(val) islogical(val) && val);
end

%%
function res = c_isField(struct,field)
 % c_isField - like isfield(), but allows fields to be specified as nested strings 
 %
 % Example:
 % a = struct('b', struct('c', struct('d', 'e')));
 % assert(c_isField(a, 'b.c.d'))
 
	if isstruct(struct)
		isf = @isfield;
	elseif istable(struct)
		isf = @(t,str) ismember(str,t.Properties.VariableNames);
	elseif isobject(struct)
		isf = @isprop;
	else
		if isempty(struct)
			res = false;
			return;
		else
			error('first input must be a struct or object');
		end
	end
	assert(ischar(field))
	i = find(field=='.',1,'first');
	if isempty(i)
		res = isf(struct,field);
	elseif ~isf(struct,field(1:i-1))
		res = false;
	elseif ~isstruct(struct) && ~istable(struct) && ~isobject(struct)
		res = false;
	elseif isempty(struct.(field(1:i-1)))
		res = false;
	else
		% recursive call
		res = c_isField(struct.(field(1:i-1)),field(i+1:end));
	end
end

%%
function c = c_structToCell(s)
% convert a struct array to a cell array of named parameter arguments

if nargin == 0 % example
	c = c_structToCell(struct('FirstArg',1,'SecondArg','example'));
	return
end

if ~isstruct(s)
	error('Input should be a struct array');
end

c = {};
fieldNames = fieldnames(s);
for i=1:length(fieldNames)
	c = [c, fieldNames{i}, {s.(fieldNames{i})}];
end

end

%%
function EEG = c_EEG_ReplaceEpochTimeSegment(EEG,varargin)
	if nargin == 0
		testfn();
		return;
	end

	p = inputParser();
	%p.addRequired('EEG',@isstruct);
	p.addParameter('timespanToReplace',[],@(x) isnumeric(x) && length(x)==2); % in s
	p.addParameter('prePostFitDurations', [], @(x) isvector(x) && (isempty(x) || length(x)==2)); % in s
	p.addParameter('eventType','',@ischar);
	p.addParameter('method','spline',@(x) ischar(x) || isscalar(x)); % if scalar, will be used as a constant to replace values
	p.addParameter('doDebug', false, @islogical);
	p.parse(varargin{:});
	s = p.Results;
	assert(isstruct(EEG));
	
	if isempty(s.timespanToReplace)
		error('Need to specify timespanToReplace');
	end

	if c_EEG_isEpoched(EEG)
		if ~isempty(s.eventType)
			indicesToReplace = false(1,EEG.pnts,EEG.trials);
			for iE = 1:EEG.trials
				ep = EEG.epoch(iE);
				epEvtIndices = find(ismember(ep.eventtype,{s.eventType}));
				for iEvt = epEvtIndices
					indicesToReplace(:,:,iE) = indicesToReplace(:,:,iE) |...
						(EEG.times >= s.timespanToReplace(1)*1e3 + ep.eventlatency{iEvt} & ...
						 EEG.times <= s.timespanToReplace(2)*1e3 + ep.eventlatency{iEvt});
					 % (note that ep.eventlatency is in ms not the usual units of samples)
				end
			end
			keyboard
		else
			% epoch around time 0 in epoch
			indicesToReplace = EEG.times >= s.timespanToReplace(1)*1e3 & EEG.times <= s.timespanToReplace(2)*1e3;
		end
		
		if strcmp(s.method,'delete')
			% don't actually interpolate, but instead delete data entirely, and fix time and other metadata
			assert(size(indicesToReplace,3)==1,'Deletion only supported when deleting common time across all trials');
			EEG.data(:,indicesToReplace,:) = [];
			EEG.times(indicesToReplace) = [];
			EEG.pnts = length(EEG.times);
		else
			for iE=1:EEG.trials
				if size(indicesToReplace,3)==1
					thisIndicesToReplace = indicesToReplace;
				else
					thisIndicesToReplace = indicesToReplace(:,:,iE);
				end
				tmp = interpolateWithinIndices(EEG.data(:,:,iE),thisIndicesToReplace,s.method,...
					'prePostFitDurations', s.prePostFitDurations*EEG.srate,...
					'doDebug', s.doDebug,...
					'indexOfTimeZero', sum(EEG.times<0),...
					'srate', EEG.srate);
				EEG.data(:,:,iE) = tmp;
			end
		end
	else
		if isempty(s.eventType)
			error('Event type must be specified if data is not epoched');
		end
		
		if strcmp(s.method,'delete')
			error('Delete not supported for continuous data');
			%TODO: implement
		end
		
		t = c_EEG_epoch_getOriginalEventTimes(EEG,'eventType', s.eventType, 'outputUnits', 's');
		
		indicesToReplace = false(1,length(EEG.times));
		
		for i=1:length(t)
			tstart = (t(i)+s.timespanToReplace(1))*1e3;
			tend = (t(i)+s.timespanToReplace(2))*1e3;
			indicesToReplace(EEG.times >= tstart & EEG.times <= tend) = true;
		end
		
		EEG.data = interpolateWithinIndices(EEG.data,indicesToReplace,s.method,...
			'prePostFitDurations', s.prePostFitDurations*EEG.srate,...
			'doDebug', s.doDebug);
	end
end


function data = interpolateWithinIndices(data,indices, varargin)
	p = inputParser();
	p.addOptional('method', 'spline', @(x) ischar(x) || isscalar(x));
	p.addParameter('prePostFitDurations', []);
	p.addParameter('doDebug', false, @islogical);
	p.addParameter('indexOfTimeZero', [], @isscalar); % only used for debug plotting
	p.addParameter('srate', [], @isscalar); % only used for debug plotting
	p.parse(varargin{:});
	s = p.Results;
	
	if ischar(s.method) && ismember(s.method, {'localSmoothedCubic', 'ARExtrapolation'})
		assert(~isempty(s.prePostFitDurations), 'prePostFitDurations must be specified for method ''%s''', s.method);
	else
		assert(isempty(s.prePostFitDurations), 'prePostFitDurations not used for method ''%s''', c_toString(s.method));
	end

	assert(length(indices)==size(data,2));
	assert(length(size(data))==2) % code below assumes data is [nchan x ntime]	
	assert(islogical(indices)); % code below assumes logical indexing
	
	if strcmp(s.method,'zero')
		s.method = 0; % convert string to num for below
	end

	times = 1:size(data,2); % arbitrary units
	knownTimes = times(~indices);
	unknownTimes = times(indices);

	
	switch(s.method)
		case 'ARExtrapolation'
			replaceStarts = find(diff(indices)>0)+1;
			replaceEnds = find(diff(indices)<0);
			assert(length(replaceStarts)==length(replaceEnds));
			
			s.prePostFitDurations = round(s.prePostFitDurations);
			assert(all(s.prePostFitDurations >= 0));
			
			prevReplaceEnd = 0;
			for iR = 1:length(replaceStarts)
				fitStart = replaceStarts(iR) - s.prePostFitDurations(1);
				if fitStart <= prevReplaceEnd
					warning('Reduced data available for fitting prior to timespanToReplace');
					fitStart = prevReplaceEnd + 1;
					assert(replaceStarts(iR) - fitStart > 1);
				end
				prevReplaceEnd = replaceEnds(iR);
				
				fitEnd = replaceEnds(iR) + s.prePostFitDurations(2);
				if iR < length(replaceStarts)
					nextStart = replaceStarts(iR+1);
				else
					nextStart = size(data, 2);
				end
				if fitEnd >= nextStart
					warning('Reduced data available for fitting after timespanToReplace');
					fitEnd = nextStart - 1;
					assert(nextStart - fitEnd > 1);
				end
				
				
				if s.doDebug
% 				toFit = data(:, fitStart:fitEnd);
					x = 1:(fitEnd-fitStart)+1;
	% 				xAll = x; % TODO: debug, delete
				
					relReplaceStart = s.prePostFitDurations(1) + 1;
					relReplaceStop = s.prePostFitDurations(1) + (replaceEnds(iR) - replaceStarts(iR)) + 1;
				
% 				toFit(:, relReplaceStart:relReplaceStop) = [];
% 				x_toInterp = x(relReplaceStart:relReplaceStop);
% 				x(relReplaceStart:relReplaceStop) = [];
				end
				
				for iCh = 1:size(data,1)
					if s.prePostFitDurations(1) > 0
						arModelOrder = ceil(s.prePostFitDurations(1)/3); % TODO: determine better way to determine this
						extrapSig_pre = c_extrapolateSignal(data(iCh, fitStart:replaceStarts(iR)-1), replaceEnds(iR)-replaceStarts(iR) + 1, arModelOrder);
						extrapSig_pre = extrapSig_pre(replaceStarts(iR)-fitStart+1:end);
					end
					
					if s.doDebug
						extrapX = x(relReplaceStart:relReplaceStop);
					end
					
					if s.prePostFitDurations(2) > 0
						arModelOrder = ceil(s.prePostFitDurations(2)/3); % TODO: determine better way to determine this
						extrapSig_post = c_extrapolateSignal(flip(data(iCh, replaceEnds(iR)+1:fitEnd)), replaceEnds(iR)-replaceStarts(iR) + 1, arModelOrder);
						extrapSig_post = flip(extrapSig_post(fitEnd-replaceEnds(iR)+1:end));
					end
					
					if all(s.prePostFitDurations > 0)
						% blend pre/post extrapolations
						if true
							% sigmoidish 
							fn = @(x, k)  1 - 1./(1+(1./x - 1).^-k);
							pre_weights = fn(linspace(1, 0, length(extrapSig_pre)), 2)';
						else
							% linear
							pre_weights = linspace(1, 0, length(extrapSig_pre))';
						end
						post_weights = 1 - pre_weights;
						extrapSig = extrapSig_pre .* pre_weights + extrapSig_post .* post_weights;
					elseif s.prePostFitDurations(1) > 0
						extrapSig = extrapSig_pre;
					elseif s.prePostFitDurations(2) > 0
						extrapSig = extrapSig_post;
					else
						error('At least one of pre and post fit durations must be greater than zero');
					end

					if s.doDebug && ismember(iCh, 1:5)
						if false
							figure; 
							labels = {};
							hp = plot(x, data(iCh, fitStart:fitEnd), 'lineWidth', 1.5); 
							hp.Color = [hp.Color 0.5];
							labels{end+1} = 'Original';
							hold on; 
							if s.prePostFitDurations(1) > 0
								hp = plot(extrapX, extrapSig_pre, 'lineWidth', 1.5);
								hp.Color = [hp.Color 0.5];
								labels{end+1} = 'Extrap pre';
								ylim(c_limits_multiply(extrema(extrapSig_pre), 2));
							end
							if s.prePostFitDurations(2) > 0
								hp = plot(extrapX, extrapSig_post, 'lineWidth', 1.5);
								hp.Color = [hp.Color 0.5];
								labels{end+1} = 'Extrap post';
								ylim(c_limits_multiply(extrema(extrapSig_post), 2));
							end
							if all(s.prePostFitDurations>0)
								hp = plot(extrapX, extrapSig, 'lineWidth', 1.5);
								hp.Color = [hp.Color 0.5];
								labels{end+1} = 'Extrap blended';
								ylim(c_limits_multiply(extrema([extrapSig_pre; extrapSig_post]), 2));
							end
							legend(labels, 'location', 'eastoutside');
						else 
							times_real = ((1:size(data,2)) - s.indexOfTimeZero)/s.srate*1e3;  % in ms
							plotStart = replaceStarts(iR) - 1.5*s.prePostFitDurations(1);
							plotEnd = replaceEnds(iR) + 1.5*s.prePostFitDurations(2);

							hf = figure;
	% 						ht = c_GUI_Tiler('numCols', 1);
							if true
								colors = [...
									0 0 0;
									0.8 0 0;
									0 0 0.8;
									0 0.6 0];
							else
								colors = c_getColors(4);
							end
							numSubplots = 5;
							has = gobjects(0);

	% 						sht = c_GUI_Tiler('parent', ht.add(), 'SideTitle', 'Original');
	% 						has(end+1) = sht.addAxes();
							has(end+1) = c_subplot(numSubplots, 1, 1);
							hp = plot(times_real(plotStart:plotEnd), data(iCh, plotStart:plotEnd), 'lineWidth', 1.5); 
							hp.Color = [colors(1,:) 1];
	% 						ylabel('Original');

	% 						sht = c_GUI_Tiler('parent', ht.add(), 'SideTitle', sprintf('Pre-stim\nextrapolation'));
	% 						has(end+1) = sht.addAxes();
							has(end+1) = c_subplot(numSubplots, 1, 2);
							% plot pre-fit span and post-fit span in lighter color
							hp = plot(times_real(plotStart:plotEnd), data(iCh, plotStart:plotEnd), 'lineWidth', 1.5); 
							hp.Color = [colors(1,:) 0.1];
							hold on;
							% plot timespan used for fitting in darker color
							hp = plot(times_real(fitStart:replaceStarts(iR)-1), data(iCh, fitStart:replaceStarts(iR)-1), 'lineWidth', 1.5);
							hp.Color = [colors(1,:) 1];
							% plot extrapolated span in different color
							hp = plot(times_real(replaceStarts(iR):replaceEnds(iR)), extrapSig_pre, 'lineWidth', 1.5);
							hp.Color = [colors(2,:) 0.5];
							ylim(c_limits_multiply(extrema(extrapSig_pre), 2));
	% 						ylabel('Pre-stim extrapolation');

	% 						sht = c_GUI_Tiler('parent', ht.add(), 'SideTitle', sprintf('Post-stim\nextrapolation'));
	% 						has(end+1) = sht.addAxes();
							has(end+1) = c_subplot(numSubplots, 1, 3);
							% plot pre-fit span and post-fit span in lighter color
							hp = plot(times_real(plotStart:plotEnd), data(iCh, plotStart:plotEnd), 'lineWidth', 1.5); 
							hp.Color = [colors(1,:) 0.1];
							hold on;
							% plot timespan used for fitting in darker color
							hp = plot(times_real(replaceEnds(iR)+1:fitEnd), data(iCh, replaceEnds(iR)+1:fitEnd), 'lineWidth', 1.5);
							hp.Color = [colors(1,:) 1];
							% plot extrapolated span in different color
							hp = plot(times_real(replaceStarts(iR):replaceEnds(iR)), extrapSig_post, 'lineWidth', 1.5);
							hp.Color = [colors(3,:) 0.5];
							ylim(c_limits_multiply(extrema(extrapSig_post), 2));
	% 						ylabel('Post-stim extrapolation');

	% 						sht = c_GUI_Tiler('parent', ht.add(), 'SideTitle', sprintf('Blending\nweights'));
	% 						has(end+1) = sht.addAxes();
							has(end+1) = c_subplot(numSubplots, 1, 4);
							y = ones(1, size(data,2));
							y(replaceStarts(iR):replaceEnds(iR)) = 0;
							hp = plot(times_real(plotStart:plotEnd), y(plotStart:plotEnd), 'lineWidth', 2);
							hp.Color = [colors(1,:), 1];
							hold on;
							hp = plot(times_real(replaceStarts(iR):replaceEnds(iR)), pre_weights, 'lineWidth', 2);
							hp.Color = [colors(2,:) 0.5];
							hp = plot(times_real(replaceStarts(iR):replaceEnds(iR)), post_weights, 'lineWidth', 2);
							hp.Color = [colors(3,:) 0.5];
							ylim([-0.1 1.1]);
	% 						ylabel('Blending weights');

	% 						sht = c_GUI_Tiler('parent', ht.add('relHeight', 1.25), 'SideTitle', sprintf('Interpolated\nresult'));
	% 						has(end+1) = sht.addAxes();
							has(end+1) = c_subplot(numSubplots, 1, 5);
							hp = plot(times_real(plotStart:plotEnd), data(iCh, plotStart:plotEnd), 'lineWidth', 1.5); 
							hp.Color = [colors(1,:) 0.1];
							hold on;
							hp = plot(times_real(plotStart:replaceStarts(iR)-1), data(iCh, plotStart:replaceStarts(iR)-1), 'lineWidth', 1.5); 
							hp.Color = [colors(1,:) 1];
							hp = plot(times_real(replaceEnds(iR)+1:plotEnd), data(iCh, replaceEnds(iR)+1:plotEnd), 'lineWidth', 1.5); 
							hp.Color = [colors(1,:) 1];
							hp = plot(times_real(replaceStarts(iR):replaceEnds(iR)), extrapSig, 'lineWidth', 1.5);
							hp.Color = [colors(4,:) 0.5];
	% 						ylabel('Interpolated result');
							xlabel('Time (ms)');

							for iA = 1:length(has)-1
								set(has(iA), 'XTick', []);
							end

							c_plot_setEqualAxes(has, 'axesToSet', 'x');
							c_plot_setEqualAxes([has(1:end-2), has(end)], 'axesToSet', 'y');

							ylim(c_limits_multiply(extrema([extrapSig_pre; extrapSig_post]), 2.5));
							xlim(extrema(times_real(plotStart:plotEnd)));
						end
					end

					data(iCh, replaceStarts(iR):replaceEnds(iR)) = extrapSig;
				end
			end
			
		case 'localSmoothedCubic'
			% (inspired by tesa_interpdata)
			% Especially with high-sampling rate data, the default behavior of spline interpolation to 
			%  use single or pairs of endpoint samples results in a spline that doesn't reduce
			%  low-frequency discontinuities in the way we want for TMS-EEG. So instead take a local amount
			%  of data at beginning and end of timespan to replace to fit a spline better to the low-frequency 
			%  component of the signal
			
			% find continuous segments to replace
			replaceStarts = find(diff(indices)>0)+1;
			replaceEnds = find(diff(indices)<0);
			assert(length(replaceStarts)==length(replaceEnds));
			
			s.prePostFitDurations = round(s.prePostFitDurations);
			assert(all(s.prePostFitDurations > 1));
			
			prevReplaceEnd = 0;
			for iR = 1:length(replaceStarts)
				fitStart = replaceStarts(iR) - s.prePostFitDurations(1);
				if fitStart <= prevReplaceEnd
					warning('Reduced data available for fitting prior to timespanToReplace');
					fitStart = prevReplaceEnd + 1;
					assert(replaceStarts(iR) - fitStart > 1);
				end
				prevReplaceEnd = replaceEnds(iR);
				
				fitEnd = replaceEnds(iR) + s.prePostFitDurations(2);
				if iR < length(replaceStarts)
					nextStart = replaceStarts(iR+1);
				else
					nextStart = size(data, 2);
				end
				if fitEnd >= nextStart
					warning('Reduced data available for fitting after timespanToReplace');
					fitEnd = nextStart - 1;
					assert(nextStart - fitEnd > 1);
				end
				
				toFit = data(:, fitStart:fitEnd);
				x = 1:(fitEnd-fitStart)+1;
				if s.doDebug
					xAll = x;
				end
				
				relReplaceStart = s.prePostFitDurations(1) + 1;
				relReplaceStop = s.prePostFitDurations(1) + (replaceEnds(iR) - replaceStarts(iR)) + 1;
				
				toFit(:, relReplaceStart:relReplaceStop) = [];
				x_toInterp = x(relReplaceStart:relReplaceStop);
				x(relReplaceStart:relReplaceStop) = [];
				
				for iCh = 1:size(data,1)
					y_toFit = toFit(iCh,:);
					if false
						% center and scale
						mu = mean(y_toFit);
						y_toFit = y_toFit - mu;
						sd = std(y_toFit);
						y_toFit = y_toFit/sd;
						p = polyfit(x, y_toFit, 3);
						data(iCh,replaceStarts(iR):replaceEnds(iR)) = polyval(p, x_toInterp)*sd + mu;
					else
						[p,~,mu] = polyfit(x, y_toFit, 3);
						if s.doDebug
							figure; 
							hp = plot(xAll, data(iCh, fitStart:fitEnd), 'lineWidth', 1.5);
							hp.Color = [hp.Color 0.5];
							hold on;
% 							hp = plot(x, toFit(iCh,:), 'linewidth', 1.5); 
% 							hp.Color = [hp.Color 0.5];
							hp = plot(xAll, polyval(p, xAll, [], mu), 'lineWidth', 1.5); 
							hp.Color = [hp.Color 0.5];
							ha = gca;
 							ha.ColorOrderIndex = ha.ColorOrderIndex+1;
							hp = plot([xAll(relReplaceStart-1), x_toInterp, xAll(relReplaceStop+1)],...
								[data(iCh, replaceStarts(iR)-1), polyval(p, x_toInterp, [], mu), data(iCh, replaceEnds(iR)+1)],...
								'linewidth', 1.5);
							hp.Color = [hp.Color 0.5];
% 							legend('Original', 'Data to fit', 'Cubic fit', 'Replacement timespan', 'location', 'eastoutside')
							legend('Original', 'Cubic fit', 'Replacement timespan', 'location', 'eastoutside')
						end
						data(iCh, replaceStarts(iR):replaceEnds(iR)) = polyval(p, x_toInterp, [], mu);
						if s.doDebug
							ylim(c_limits_multiply(extrema(data(iCh, fitStart:fitEnd)), 1.5));
							keyboard
						end
					end
					
				end
			end
		otherwise
			if ischar(s.method)
				if s.doDebug
						tmp = interp1(knownTimes,data(:,~indices).',unknownTimes,s.method).';
						replaceStart = find(indices>0, 1, 'first');
						replaceStop = find(diff(indices)<0, 1, 'first');
						
						nearbyStart = replaceStart - 500;
						nearbyEnd = replaceStop + 500;
						
						iCh = 1;
						
						xAll = 1:nearbyEnd-nearbyStart+1;
						xReplace = (replaceStart:replaceStop)-nearbyStart+1;
						
						figure; 
						hp = plot(xAll, data(iCh, nearbyStart:nearbyEnd), 'lineWidth', 1.5);
						hp.Color = [hp.Color 0.5];
						hold on;
						hp = plot(xReplace, tmp(iCh, 1:replaceStop-replaceStart+1), 'lineWidth', 1.5); 
						hp.Color = [hp.Color 0.5];
						ylim(c_limits_multiply(extrema(tmp(iCh, 1:replaceStop-replaceStart+1)), 4));
						legend('Original', sprintf('%s interpolation', s.method), 'location', 'eastoutside')
						keyboard
						
						data(:, indices) = tmp;
				else
						data(:,indices) = interp1(knownTimes,data(:,~indices).',unknownTimes,s.method).';
				end
			elseif isscalar(s.method)
				% method is actually a constant scalar to use to replace all unknown values
				data(:,indices) = s.method;
			else
				error('Invalid method');
			end
	end
end

function [extrapSig,modelCoeff] = c_extrapolateSignal(sig,numAddedPts,modelOrder,modelCoeff)
	if isvector(sig) && size(sig,2) > 1
		sig = sig';
	end
	if nargin < 3
		modelOrder = floor(size(sig,1)-1);
	end
	if nargin < 4
		modelCoeff = arburg(sig,modelOrder);
	end
	
	extrapSig = zeros(size(sig,1) + numAddedPts,size(sig,2));
	
	extrapSig(1:size(sig,1),:) = sig; % do not extrapolate what we already know
	
	if any(isnan(modelCoeff))
		if all(sig==0)
			% special case: if input is all zeros, arburg will return [1 NaN NaN ...] for coefficients
			assert(modelCoeff(1)==1);
			assert(all(isnan(modelCoeff(2:end))));
			extrapSig((size(sig,1)+1):end,:) = 0;
			return;
		elseif c_allEqual(sig)
			% special case: if all constant values, arburg will return [1 NaN NaN ...]
			assert(modelCoeff(1)==1);
			assert(all(isnan(modelCoeff(2:end))));
			extrapSig((size(sig,1)+1):end,:) = extrapSig(1);
			return;
		else
			error('NaN coefficients in fitted model');
			% perhaps can happen if model order is too high, maybe add support by pruning tailing NaNs from coefficients
		end
	end
	
	[~,zf] = filter(-[0 modelCoeff(2:end)], 1, sig,[],1);
	extrapSig((size(sig,1)+1):end,:) = filter([0 0], -modelCoeff, zeros(size(extrapSig,1)-size(sig,1),size(sig,2)), zf,1);
end


function testfn()	
	close all;
	tmin = -0.2;
	tmax = 0.4;
	srate = 1000;
	N = (tmax-tmin)*srate;
	t = linspace(tmin,tmax,N);
	x(1,:) = sin(10*pi*t);
	x(2,:) = cos(25*pi*t);
	x(3,:) = cumsum(randn(1,length(t)));
	x(3,:) = x(3,:) / max(abs(x(3,:)));
	
	timespanToReplace = [-0.01 0.05];
	indicesToReplace = t>=timespanToReplace(1) & t<timespanToReplace(2);
	
	hf = figure;
	c_subplot(2,1,1);
	plot(t,x);
	title('Original signals');
	
	ha = c_subplot(2,1,2);
	
	debugArgs = {'doDebug', true, 'indexOfTimeZero', sum(t<0), 'srate', srate};
	
% 	y = interpolateWithinIndices(x,indicesToReplace,'localSmoothedCubic', ...
% 		'prePostFitDurations', [100 100]*1e-3*srate);
	y = interpolateWithinIndices(x,indicesToReplace,'ARExtrapolation', ...
		'prePostFitDurations', [100 100]*1e-3*srate,...
		debugArgs{:});
% 	y = interpolateWithinIndices(x,indicesToReplace,'makima');
	plot(ha, t,y);
	title(ha, sprintf('Signals with %s interpolated', c_toString(timespanToReplace)));
end

%%
function isEpoched = c_EEG_isEpoched(EEG)
	if false
		isEpoched = EEG.trials > 1;
	else
		% handle special case where data was epoched to a single epoch (e.g. by trial rejection)
		% (infer that if there is an event exactly at time 0 and there is data available for negative times,
		%  then this was probably a single epoch)
		isEpoched = EEG.trials > 1 || (EEG.xmin < 0 && all(abs(EEG.times(round(mod([EEG.event.latency], EEG.pnts))))<1e-10)); 
	end
end


%%
function out = c_if(condition, valIfTrue, valIfFalse)
% c_if - evaluates condition to determine which of two values to return
% Useful for conditional statements inside anonymous functions or other shorthand conditions
%
% Example:
%	c_if(rand(1) > 0.5,'Is large','Is small')

if condition
	out = valIfTrue;
else
	out = valIfFalse;
end
end

%%
function EEG = c_EEG_filter_butterworth(EEG, varargin)
% adapted from TESA tesa_filtbutter function

p = inputParser();
p.addRequired('cutoffFreqs',@(x) isnumeric(x) && isvector(x) && length(x)==2); % in Hz
p.addParameter('order',4,@(x) isscalar(x) && mod(x,2)==0);
p.addParameter('type','auto',@ischar);
p.parse(varargin{:});
s = p.Results;
assert(isstruct(EEG));

switch(s.type)
	case 'auto'
		if all(s.cutoffFreqs>0)
			s.type = 'bandpass';
		elseif s.cutoffFreqs(1) > 0
			s.type = 'high';
			s.cutoffFreqs = s.cutoffFreqs(1);
		elseif s.cutoffFreqs(2) > 0
			s.type = 'low';
			s.cutoffFreqs = s.cutoffFreqs(2);
		else
			% both cutoffs are 0 or < 0
			warning('Cutoff frequencies invalid: %s. Not filtering.',c_toString(s.cutoffFreqs));
			return;
		end
	case 'bandpass'
		% do nothing
	case 'bandstop'
		s.type = 'stop';
	case 'lowpass'
		s.type = 'low';
		s.cutoffFreqs = s.cutoffFreqs(2);
	case 'highpass'
		s.type = 'high';
		s.cutoffFreqs = s.cutoffFreqs(1);
	otherwise
		error('Invalid type: %s',s.type);
end

if length(s.cutoffFreqs) > 1 && s.cutoffFreqs(1) > s.cutoffFreqs(2)
	error('Cutoff frequencies not in ascending order. Swapped?')
end

[z,p] = butter(s.order/2, s.cutoffFreqs./(EEG.srate/2), s.type);

% temporarily move time dimension to first dimension for filtfilt
EEG.data = permute(EEG.data,[2 1 3]);

didConvertFromSingle = false;
if isa(EEG.data,'single')
	c_say('Temporarily converting EEG.data from single to double');
	EEG.data = double(EEG.data);
	c_sayDone();
	didConvertFromSingle = true;
end

% apply filter
EEG.data = filtfilt(z,p,EEG.data);

if didConvertFromSingle
	c_say('Converting EEG.data back from double to single');
	EEG.data = single(EEG.data);
	c_sayDone();
end

% undo rearrangement of dimensions
EEG.data = ipermute(EEG.data,[2 1 3]);

end

%%
function c_say(varargin)
% c_say - Wrapper around fprintf with other added features
% The collection of c_say* functions can be used to print formatted strings  
% with indication of level of nesting, e.g. to show that a sequence of printed lines all nest within
% a particular function call.
%
% Syntax matches that of fprintf(), except that newlines are automatically added at end of string
%
% Example:
% 	c_say('Outer print')
% 	c_saySingle('Indented print');
% 	c_say('Further indenting');
% 	c_saySingle('Indented print with formatting: %.3f',pi)
% 	c_sayDone('End of inner');
% 	c_sayDone('End of outer');
%
% See also: c_sayDone, c_saySingle, c_saySingleMultiline, c_sayResetLevel

	global sayNestLevel;
	if isempty(sayNestLevel)
		sayNestLevel = 0;
	end
	
	global saySilenceLevel;
	if ~isempty(saySilenceLevel) && sayNestLevel >= saySilenceLevel
		% don't print anything
		sayNestLevel = sayNestLevel + 1;
		return
	end
	
	global sayDateFormat;
	if isempty(sayDateFormat)
		sayDateFormat = 'HH:MM:ss';
	end
	
	if verLessThan('matlab','8.4')
		fprintf('%s ',datestr(now,13));
	else
		fprintf('%s ',datestr(datetime,sayDateFormat));
	end

	for i=1:sayNestLevel
		if mod(i,2)==0
			fprintf(' ');
		else
			fprintf('|');
		end
	end
	sayNestLevel = sayNestLevel + 1;
	fprintf('v ');
	fprintf(varargin{:});
	fprintf('\n');
end

%%
function c_sayDone(varargin)
% c_sayDone - Wrapper around fprintf with other added features
%
% See also: c_say

global sayNestLevel;
if isempty(sayNestLevel)
	sayNestLevel = 0;
end

global saySilenceLevel;
if ~isempty(saySilenceLevel) && sayNestLevel-1 >= saySilenceLevel
	% don't print anything
	sayNestLevel = max(sayNestLevel-1,0);
	return
end


global sayDateFormat;
if isempty(sayDateFormat)
	sayDateFormat = 'HH:MM:ss';
end

if verLessThan('matlab','8.4')
	fprintf('%s ',datestr(now,13));
else
	fprintf('%s ',datestr(datetime,sayDateFormat));
end

sayNestLevel = max(sayNestLevel-1,0);
for i=1:sayNestLevel
	if mod(i,2)==0
		fprintf(' ');
	else
		fprintf('|');
	end
end
fprintf('^ ');
if ~isempty(varargin)
	fprintf(varargin{:});
end
fprintf('\n');

end

%%
function h = c_subplot(varargin)
% c_subplot - wrapper around subplot() with added features
% Beyond normal subplot features and alternate input formats, adds a right-click context menu
%  option to "pop out" specific subplots into their own figure.
%
% Input can be:
%	c_subplot(numRows, numColumns, row, column)
% or 
%	c_subplot(figNum, numFigs)
% or 
%	c_subplot(numRows, numColumns, index)
% or
%	c_subplot('position',[left bottom width height])

if nargin == 0, testfn(); return; end

% p = inputParser();
% p.addOptional('num1',[],@isscalar);
% p.addOptional('num2',[],@isscalar);
% p.addOptional('num3',[],@isscalar);
% p.addOptional('num4',[],@isscalar);
% p.addParameter('position',[],@(x) isvector(x) && length(x)==4);
% p.addParameter('numRows',[],@isscalar);
% p.addParameter('numCols',[],@isscalar);
% p.addParameter('rowNum',[],@isscalar);
% p.addParameter('colNum',[],@isscalar);
% p.addParameter('plotNum',[],@isscalar);
% p.addParameter('numPlots',[],@isscalar);
% p.addParameter('Number',[],@isscalar); % legacy syntax, equivalent to 'plotNum'
% p.parse(varargin{:});
% s = p.Results;
% 
% if ~isempty(s.num1)
% 	% if using non-named parameter syntax, no conflicting named parameters should be specified
% 	assert(all(ismember({'numCols','numRows','colNum','rowNum','plotNum','Number'},p.UsingDefaults)));
% end
% if ~isempty(s.num4)
% 	s.numRows = s.num1;
% 	s.numCols = s.num2;
% 	s.rowNum = s.num3;
% 	s.colNum = s.num4;
% elseif ~isempty(s.num3)
% 	s.numRows = s.num1;
% 	s.numCols = s.num2;
% 	s.plotNum = s.num3;
% elseif ~isempty(s.num2)
% 	s.plotNum = s.num1;
% 	s.numPlots = s.num2;
% elseif ~isempty(s.num1)
% 	error('Invalid syntax: only one un-named parameter specified');
% end
% if ~isempty(s.Number)
% 	% handle legacy syntax ('Number' instead of 'plotNum')
% 	assert(isempty(s.plotNum));
% 	s.plotNum = s.Number;
% end
%TODO: finish converting to use inputParser
	
extraArgs = {};
if isscalar(varargin{1}) || (ischar(varargin{1}) && strcmpi(varargin{1},'position'))
	numRecognizedArgs = 1;
else
	error('invalid input');
end
for i=2:nargin
	if isnumeric(varargin{i}) || (ischar(varargin{i}) && strcmpi(varargin{i},'Number'))
		numRecognizedArgs = numRecognizedArgs + 1;
	else
		% start of "extra" args reached
		break;
	end
end
if numRecognizedArgs < nargin
	extraArgs = varargin(numRecognizedArgs+1:end);
	assert(mod(length(extraArgs),2)==0); % make sure there are even number of extra args (assuming they are all (name,value) pairs
end
subVarargin = varargin(1:numRecognizedArgs);

if numRecognizedArgs >= 1
	arg1 = subVarargin{1};
end
if numRecognizedArgs >= 2
	arg2 = subVarargin{2};
end
if numRecognizedArgs >= 3
	arg3 = subVarargin{3};
end
if numRecognizedArgs >= 4
	arg4 = subVarargin{4};
end




if numRecognizedArgs >= 2 && ischar(arg1)
	% (name, value) arguments
	
	p = inputParser;
	p.addParameter('Position',[0 0 1 1],@isvector);
	p.addParameter('Number',0,@isscalar);
% 		p.addParameter('FrameWidth',0.05,@isscalar);
% 		p.addParameter('FrameHeight',0.05,@isscalar);
% 		p.addParameter('LeftInset',0.05,@isscalar);
% 		p.addParameter('BottomInset',0.05,@isscalar);

	p.parse(subVarargin{:});
	
	l = p.Results.Position(1);
	b = p.Results.Position(2);
	w = p.Results.Position(3);
	h = p.Results.Position(4);
% 	fw = p.Results.FrameWidth;
% 	fh = p.Results.FrameHeight;
% 	li = p.Results.LeftInset;
% 	bi = p.Results.BottomInset;
% 	
% 	al = (l+fw/2)*(1-li)+li;
% 	ab = (b+fh/2)*(1-bi)+bi;
% 	aw = (w-fw)*(1-li);
% 	ah = (h-fh)*(1-bi);
	
	if 1
		argsToSubplot = {};
		h = axes('OuterPosition',[l b w h],'ActivePositionProperty','outerposition',extraArgs{:});
	else
		argsToSubplot = {'Position',[l b w h]};
	end
	index = p.Results.Number;
	
elseif numRecognizedArgs == 2
		figNum = arg1;
		numFigs = arg2;

		if numFigs == 3
			numCols = numFigs;
			numRows = 1;
		elseif numFigs == 6
			numCols = 3;
			numRows = 2;
		elseif numFigs == 8
			numCols = 4;
			numRows = 2;
		else
			numCols = ceil(sqrt(numFigs));
			numRows = ceil(numFigs/numCols);
		end

		index = figNum;
		
		argsToSubplot = {numRows, numCols, index};
		
elseif numRecognizedArgs == 3
	index = arg3; % third argument is actually index (i.e. same args as normal subplot() )
	argsToSubplot = {arg1, arg2, index};
else
	index = (arg3-1)*arg2 + arg4;
	argsToSubplot = {arg1, arg2, index};
end

if ~isempty(argsToSubplot)
	h = subplot(argsToSubplot{:},extraArgs{:});
end
figHandle = get(h,'parent');
c_figure_addInteractive(figHandle,h,index);

ud = get(h,'UserData');
if iscell(ud)
	index = find(cellfun(@isstruct,ud),1,'first');
	if isempty(index)
		index = length(ud)+1;
		ud{index} = struct();
	end
	ud{index}.SubplotIndex = index;
else
	ud.SubplotIndex = index;
end
set(h,'UserData',ud);

end

function c_figure_addInteractive(figHandle,subplotHandle,index)

	% Define a context menu; it is not attached to anything
	hcmenu = get(figHandle,'uicontextmenu');
	if isempty(hcmenu)
		hcmenu = uicontextmenu;
	end
	% Define callbacks for context menu items that change linestyle
	% Define the context menu items and install their callbacks
	item1 = uimenu(hcmenu, 'Label', ['Pop out ' num2str(index)], 'Callback', {@c_copySubplotToNewFig,subplotHandle,num2str(index)});
	set(figHandle,'uicontextmenu',hcmenu);
end


function c_copySubplotToNewFig(cb,eventdata,handle,label)
	parentHandle = get(handle,'parent');
	if isprop(parentHandle,'Name')
		parentName = get(parentHandle,'Name');
	elseif isprop(parentHandle,'Title')
		parentName = get(parentHandle,'Title');
	else
		parentName = '';
	end
	figName = [parentName ' subplot ' label];
	hh = copyobj(handle,figure('name',figName));
	parentFields = get(parentHandle);
	if isfield(parentFields,'Colormap')
		colormap(hh,parentFields.Colormap);
	end
	%resize the axis to fill the figure
	set(hh, 'Position', get(0, 'DefaultAxesPosition'));
end

%%
function c_plot_setEqualAxes(varargin)
% c_plot_setEqualAxes - set axis limits equal across multiple plots
% By default, this included color scales, and also links the axes so that changes to one plot
%  affect all linked plots
%
% Example:
% 	figure;
% 	h1 = c_subplot(1,2);
% 	plot(rand(10,2));
% 	h2 = c_subplot(2,2);
% 	plot(rand(20,2)*2);
% 	c_plot_setEqualAxes([h1 h2]);

p = inputParser();
p.addOptional('axisHandles',[],@(x) all(ishandle(x(:))));
p.addParameter('xlim',[nan, nan],@isvector);
p.addParameter('ylim',[nan, nan],@isvector);
p.addParameter('zlim',[nan, nan],@isvector);
p.addParameter('clim',[nan, nan],@isvector);
p.addParameter('axesToSet','xyzc',@(x) ischar(x) || iscellstr(x)); % if cell, one char per axisHandle in each element
p.addParameter('doForceSymmetric',false,@islogical);
p.addParameter('doForceEqualAspect',false,@islogical);
p.addParameter('doLink',true,@islogical);
p.parse(varargin{:});
s = p.Results;

if isempty(s.axisHandles)
	% no handles given, assume we should grab all axes from current figure
	s.axisHandles = gcf;
else
	s.axisHandles = s.axisHandles(:); % reshape any higher dimensions into vector
end

axisHandles = [];
for i=1:length(s.axisHandles)
	if isgraphics(s.axisHandles(i),'Figure')
		% figure handle given, assume we should grab all axes from the figure
		childHandles = findobj(s.axisHandles(i),'Type','axes');
		axisHandles = [axisHandles; childHandles];
	elseif isgraphics(s.axisHandles(i),'axes')
		axisHandles = [axisHandles;s.axisHandles(i)];
	else
		error('invalid handle');
	end
end
s.axisHandles = axisHandles;

if s.doForceEqualAspect
	set(s.axisHandles,'DataAspectRatio',[1 1 1]);
end

for j=1:length(s.axesToSet)
	if length(s.doForceSymmetric)==1 % one doForceSymmetric value for all axes
		doForceSymmetric = s.doForceSymmetric;
	else % one doForceSymmetric value for each axis
		assert(length(s.doForceSymmetric)==length(s.axesToSet)); 
		doForceSymmetric = s.doForceSymmetric(j);
	end
	if ~iscell(s.axesToSet)
		limfield = [lower(s.axesToSet(j)) 'lim'];
		lim = s.(limfield);
		setEqualAxis(s.axisHandles,s.axesToSet(j),lim,doForceSymmetric,s.doLink);
	else
		assert(length(s.axisHandles)==length(s.axesToSet{j}));
		lim = s.([lower(s.axesToSet{j}(1)) 'lim']);
		setEqualAxis(s.axisHandles,s.axesToSet{j},lim,doForceSymmetric,s.doLink);
	end
end

end


function setEqualAxis(axisHandles,axesToSet,lim,doSymmetry, doLink)
	if length(axesToSet)==1
		assert(ischar(axesToSet));
		axesToSet = repmat(axesToSet,1,length(axisHandles));
	end
	assert(length(axesToSet)==length(axisHandles));
	fieldsOfInterest = cell(1,length(axisHandles));
	for i=1:length(axisHandles)
		fieldsOfInterest{i} = [upper(axesToSet(i)) 'Lim'];
	end
	i=1;
	currentLim = get(axisHandles(i),fieldsOfInterest{i});
	if isdatetime(currentLim)
		% handle special case of datetime axis values
		assert(all(arrayfun(@(i) isdatetime(get(axisHandles(i),fieldsOfInterest{i})), 1:length(axisHandles))));
		if all(isnan(lim))
			i=1;
			[minVal,maxVal] = c_mat_deal(get(axisHandles(i),fieldsOfInterest{i}));
			for i=2:length(axisHandles)
				newLim = get(axisHandles(i),fieldsOfInterest{i});
				minVal = min(minVal,newLim(1));
				maxVal = max(maxVal,newLim(2));
			end
			lim = [minVal, maxVal];
		elseif any(isnan(lim))
			error('Partially specified datetime limits not currently supported');
		end
	else
		isnat_ = @(t) isa(t, 'datetim') && isnat(t);
		if isnan(lim(1)) || isnat_(lim(1)) % need to autodetect min
			minVal = inf;
			for i=1:length(axisHandles)
				newVal = paren(get(axisHandles(i),fieldsOfInterest{i}),1);
				minVal = min(newVal,minVal);
			end
			lim(1) = minVal;
		end
		if isnan(lim(2)) || isnat_(lim(2)) % need to autodetect max
			maxVal = -inf;
			for i=1:length(axisHandles)
				newVal = paren(get(axisHandles(i),fieldsOfInterest{i}),2);
				maxVal = max(newVal,maxVal);
			end
			lim(2) = maxVal;
		end
	end
	
	if doSymmetry
		lim = [-1 1]*max(abs(lim));
	end
	
	% set limits
	for i=1:length(axisHandles)
		set(axisHandles(i),fieldsOfInterest{i},lim);
	end
	
	if doLink
		if c_allEqual(axesToSet)
			fieldOfInterest = fieldsOfInterest{1};
			hlink = linkprop(axisHandles,fieldOfInterest);
			ud = get(axisHandles(1),'UserData');
			if iscell(ud) && ~isempty(ud)
				ud{end+1} = ['Link' fieldOfInterest 'Handle'];
				ud{end+1} = hlink;
			else
				ud.(['Link' fieldOfInterest 'Handle']) = hlink;
			end
			set(axisHandles(1),'UserData',ud);
		else
			% linking not currently supported for linking between different axes
		end
	end
end

%%
function out = paren(x, varargin)
	out = x(varargin{:});
end

%%
function allequal = c_allEqual(varargin)
% c_allEqual - similar to isequal() but allows arbitrary number of inputs instead of pairwise
%  (and uses isequaln instead of isequal to allow for NaNs)
%
% Examples:
%	c_allEqual(true,true,false)
%	c_allEqual(true,true,true)
%	c_allEqual([true,true,false])

	
	if length(varargin)==1
		allequal = true;
		for i=2:length(varargin{1})
			allequal = allequal && isequaln(varargin{1}(1),varargin{1}(i));
			if ~allequal
				break;
			end
		end
		return;
	end
	
	assert(length(varargin)>1);
	
	allequal = true;
	for i=2:nargin
		allequal = allequal && isequaln(varargin{1},varargin{i});
		if ~allequal
			break;
		end
	end
end

%%
function [extremeVals imin imax] = extrema(varargin)
% extrema - calculate min and max values
%
% Syntax:
%   [minmax, imin, imax] = extrema(vals)
%   [minmax, imin, imax] = extrema(vals,[],dim)
%
% Emulates syntax of min() and max(), see those functions for details.
%
% Inputs:
%	vals - array of values
%   [] - second output not used, for compatibility with min() and max() syntax
%   dim - dimension of vals along which to operate
%
% Examples:
%   extremeVals = extrema([1 2 3])
%   extremeVals = extrema(rand(2,10),[],2)
	
	if nargin > 2
		assert(isempty(varargin{2}));
	end

	if nargout >= 3
		[minval, imin] = min(varargin{:});
	else
		minval = min(varargin{:});
	end
	if nargout >= 4
		[maxval, imax] = max(varargin{:});
	else
		maxval = max(varargin{:});
	end
	
	if nargin==3 && varargin{3} ~= 1
		extremeVals = [minval,maxval];
	else
		extremeVals = [minval.', maxval.'];
	end
end