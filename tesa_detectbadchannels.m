% tesa_detectbadchannels() - Detect and optionally replace bad EEG channels
%
% Usage:
%   >> [EEG, misc] = tesa_detectbadchannels(EEG);
%   >> [EEG, misc] = tesa_detectbadchannels(EEG, 'key', value, ...);
%
% Description:
%   This function detects bad EEG channels using several possible methods,
%   including a PREP-style robust deviation metric and data-driven Wiener
%   noise estimates (DDWiener; adapted from the TESA toolbox).
%   Identified bad channels can be interpolated, removed, replaced with NaNs,
%   or left unchanged.
%
%   The function is designed for use within EEGLAB/TESA preprocessing
%   pipelines and supports optional exclusion of artifact-contaminated
%   time windows during detection.
%
%   The function is adapted from Chris Cline's AARATEP toolbox.
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%
% Optional inputs (key-value pairs):
%
%   'detectionMethod'   - Method used to identify bad channels.
%                         String or cell array of strings. Supported options:
%                           'PREP_deviation'  (default)
%                           'TESA_DDWiener'
%                           'TESA_DDWiener_PerTrial'
%                           'TESA_DDWiener_IgnoreArtifactTime'
%                           'TESA_DDWiener_PerTrial_IgnoreArtifactTime'
%                           'TESA_DDWiener_PerTrial_BaselineOnly'
%                           'fromASR'
%
%                         If a cell array is supplied, a channel is marked as
%                         bad if it is detected by ANY listed method.
%
%   'replaceMethod'     - How to handle detected bad channels:
%                           'interpolate' (default) - spherical spline interpolation
%                           'remove' or 'delete'   - remove channels from data
%                           'NaN'                   - replace channel data with NaNs
%                           'none'                  - detect only, do not modify EEG
%
%   'artifactTimespan'  - [start end] time window (in ms) to ignore during
%                         bad channel detection (e.g., TMS pulse artifact).
%                         Behavior depends on detection method:
%                           - PREP_deviation: values set to zero
%                           - DDWiener *_IgnoreArtifactTime: values set to zero
%                           - *_BaselineOnly: uses data prior to artifact onset
%
%   'threshold'         - Numeric threshold for bad channel detection.
%                         Default values are method-specific:
%                           PREP_deviation: 9 (robust SD units)
%                           DDWiener methods: 20 (MAD units above median)
%
% Outputs:
%   EEG                 - Updated EEGLAB EEG structure. Additional fields added:
%                           EEG.chansAll - Original EEG.chanlocs saved for later interpolation
%                           EEG.badChans - Cell array of channels detected as bad
%
%   misc                - Structure containing diagnostic information:
%                           misc.badChannelIndices
%                               Indices of channels detected as bad
%                               (referring to original channel order)
%
%                           misc.channelScores
%                               Per-channel detection scores
%                               (robust deviation or Wiener noise estimate)
%
%                           misc.scoreThreshold
%                               Threshold used for classification
%
% Notes:
%   • When multiple detection methods are supplied, channel replacement
%     is performed only once after all methods are evaluated.
%   • If all channels are marked as bad, interpolation is aborted.
%   • The 'fromASR' option assumes that ASR has already been run and that
%     channel rejection information is stored in EEG.etc.clean_channel_mask.
%   • Fields EEG.chansAll and EEG.badChans are automatically added if not present.
%
% Example:
%   % Detect and interpolate bad channels using PREP-style deviation
%   [EEG, misc] = tesa_detectbadchannels(EEG, ...
%       'detectionMethod', 'PREP_deviation', ...
%       'replaceMethod', 'interpolate');
%
%   % Combine PREP deviation and DDWiener per-trial detection
%   [EEG, misc] = tesa_detectbadchannels(EEG, ...
%       'detectionMethod', {'PREP_deviation','TESA_DDWiener_PerTrial'}, ...
%       'artifactTimespan', [-2 10]);
%
% See also:
%   pop_clean_rawdata, pop_interp, eeg_interp

% This script was adapted by Nigel Rogasch for the TESA toolbox. Original
% code is available from:
% https://github.com/chriscline/AARATEPPipeline

% Acknowledgement:
%   Parts of this function and its documentation were developed with
%   assistance from ChatGPT (OpenAI), used for code review, refactoring,
%   and documentation support. All scientific and methodological decisions
%   were made by the authors.

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

function [EEG, misc] = tesa_detectbadchannels(EEG,varargin)
%
% output struct `misc` includes the following:
% - misc.badChannelIndices: numeric indices of bad channels 
%   (if removing bad channels, these indices index into the channels prior to rejection)
% - misc.channelScores: per-ch scores, where greater magnitudes indicate more likely to be bad channels 
%		if detectionMethod==PREP_deviation, then these scores are robustChanDeviation scores

if nargin < 1
	error('Not enough input arguments.');
end

p = inputParser();
% p.addRequired('EEG', @isstruct);
p.addParameter('detectionMethod', 'PREP_deviation', @(x) ischar(x) || iscellstr(x));
p.addParameter('replaceMethod', 'interpolate', @ischar);
p.addParameter('artifactTimespan', [], @c_isSpan);
p.addParameter('threshold', [], @isscalar); % default value for this is method-specific
% p.addParameter('doPlot', false, @islogical);
p.parse(varargin{:});
s = p.Results;
s.EEG = EEG;

detectionMethod = s.detectionMethod;

% --- Save current channel structure for later interpolation ---
if ~isfield(EEG, 'chansAll') || isempty(EEG.chansAll)
    EEG.chansAll = EEG.chanlocs;
end

if iscellstr(s.detectionMethod)
	% if multiple methods specified, reject a channel that is identified as bad by any of the listed methods

	% if s.doPlot
	% 	misc.hf = figure;
	% 	ht = c_GUI_Tiler();
	% end

	badChannelIndices = false(EEG.nbchan, 1);
	misc.channelScores = nan(EEG.nbchan, length(s.detectionMethod));
	misc.scoreThreshold = repmat({}, 1, length(s.detectionMethod));

	EEG_tmp = EEG;

	for iDM = 1:length(s.detectionMethod)
		% kwargs = c_cellToStruct(varargin(2:end));
        kwargs = c_cellToStruct(varargin(1:end));
		kwargs.detectionMethod = s.detectionMethod{iDM};
		%kwargs.replaceMethod = 'none';
		kwargs = c_structToCell(kwargs);
		[EEG_tmp, misc_iDM] = tesa_detectbadchannels(EEG_tmp, kwargs{:});

		badChannelIndices(misc_iDM.badChannelIndices) = true;
		misc.channelScores(:, iDM) = misc_iDM.channelScores;
		misc.scoreThreshold{iDM} = misc_iDM.scoreThreshold;

		% if s.doPlot
		% 	hp = ht.add();
		% 	copyobj(misc_iDM.hf.Children(2:end), hp)
		% 	title(gca, s.detectionMethod{iDM})
		% 	close(misc_iDM.hf);
		% end
	end

	badChannels = find(badChannelIndices);

	s.detectionMethod = '__multiple__';
end

% TODO: also do artifact timespan interpolation for TESA_DDWiener* methods or issue warning if specifying a timespan that is ignored
if ismember(s.detectionMethod, '__multiple__')
	% handled within each method, do nothing here
elseif ismember(s.detectionMethod, {'PREP_deviation', 'TESA_DDWiener_PerTrial_IgnoreArtifactTime', 'TESA_DDWiener_PerTrial_BaselineOnly'})
	if ismember('artifactTimespan', p.UsingDefaults)
		warning('Artifact timespan unspecified. Will not ignore artifact timespan for bad channel detection');
	end
	if ~isempty(s.artifactTimespan)
		if ismember(s.detectionMethod, {'PREP_deviation'})
			% tmpEEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
			% 	'timespanToReplace', s.artifactTimespan,...
			% 	'method', 'zero');
            tmpEEG = tesa_removedata( EEG, s.artifactTimespan );
		elseif ismember(s.detectionMethod, {'TESA_DDWiener_IgnoreArtifactTime', 'TESA_DDWiener_PerTrial_IgnoreArtifactTime'})
			% for backwards compatibility, don't do artifact ignoring for TESA_DDWiener and TESA_DDWiener_PerTrial
			% methods unless explicitly suffixed with '_IgnoreArtifactTime'
			% tmpEEG = c_EEG_ReplaceEpochTimeSegment(EEG,...
			% 	'timespanToReplace', s.artifactTimespan,...
			% 	'method', 'delete');
            tmpEEG = tesa_removedata( EEG, s.artifactTimespan );
		elseif ismember(s.detectionMethod, {'TESA_DDWiener_PerTrial_BaselineOnly'})
			tmpEEG = pop_select(EEG, 'notime', [s.artifactTimespan(1), EEG.xmax]); 
		else
			error('not implemented')
		end
	else
		tmpEEG = EEG;
	end
else
	if ~isempty(s.artifactTimespan)
		warning('Artifact timespan argument is ignored for method %s', s.detectionMethod)
	end
	tmpEEG = EEG;
end


switch(s.detectionMethod)
	case '__multiple__'
		% handled above, do nothing here
	case 'fromASR'
		% if ASR was previously run on EEG, bad channels are recorded in EEG struct
		% (requires that chanlocs prior to ASR were saved as EEG.urchanlocs)
		c_say('Determining bad channels from ASR results');
		
		assert(isempty(s.threshold), 'Threshold not used for channels previously marked for rejection by ASR');
		
		if ~c_isFieldAndNonEmpty(EEG.etc, 'clean_channel_mask')
			goodChannels = 1:EEG.nbchan;
			badChannels = [];
			assert(length(EEG.chanlocs)==length(EEG.urchanlocs));
		else
			goodChannels = find(EEG.etc.clean_channel_mask);
			badChannels = find(~EEG.etc.clean_channel_mask);
		end
		
		if ~isempty(badChannels)
			if ~ismember(s.replaceMethod, {'interpolate', 'NaN'})
				error('Detection without rejection not implemented for fromASR method, since "bad" indices are already removed by ASR')
			end
			
			c_say('Inserting NaNs for bad channels prior to interpolation');
			nbchan = length(goodChannels)+length(badChannels);
			assert(length(EEG.urchanlocs)==nbchan);
			tmpData = nan(nbchan, EEG.pnts, EEG.trials);
			tmpData(goodChannels, :, :) = EEG.data;
			EEG.data = tmpData;
			EEG.chanlocs = EEG.urchanlocs;
			EEG.nbchan = nbchan;
			c_sayDone();
			% will be interpolated below
		end
		
		% if s.doPlot
		% 	keyboard % TODO
		% end
		
	case 'PREP_deviation'
		c_say('Detecting bad channels based on PREP deviation scores');
		% Note: only catches very obvious bad channels that may degrade 1st stage ICA performance
		% (method adapted from part of PREP pipeline)
		bad = struct();

		if isempty(s.threshold)
			s.threshold = 9;
			c_saySingle('Using default rejection threshold of %s', c_toString(s.threshold));
		else
			c_saySingle('Using specified rejection threshold of %s', c_toString(s.threshold));
		end

		tmpData = reshape(tmpEEG.data, [EEG.nbchan, EEG.pnts*EEG.trials]);

		% detect nan channels
		bad.nan = any(isnan(tmpData),2);

		% detect constant channels
		bad.constant = mad(tmpData, 1, 2) < 10e-10;

		% unusually high or low amplitude
		chanDeviation = 0.7413 * iqr(tmpData, 2);
		chanDeviationSD = 0.7413 * iqr(chanDeviation);
		chanDeviationMedian = nanmedian(chanDeviation);

		robustChanDeviation = (chanDeviation - chanDeviationMedian) / chanDeviationSD;
		bad.deviation = abs(robustChanDeviation) > s.threshold;
		
		misc.channelScores = robustChanDeviation;
		misc.scoreThreshold = [-1 1]*s.threshold;

		badChannels = find(bad.nan | bad.constant | bad.deviation);
		
		% if s.doPlot
		% 	misc.hf = figure; 
		% 	topoplot(misc.channelScores, EEG.chanlocs, 'electrodes', 'on',...
		% 		'emarker2',{badChannels,'x','k',10,2});
		% 	hc = colorbar; 
		% 	caxis([-1 1]*s.threshold); 
		% 	ylabel(hc,'Channel deviation');
		% end
		
	case {'TESA_DDWiener', 'TESA_DDWiener_PerTrial',...
			'TESA_DDWiener_IgnoreArtifactTime', 'TESA_DDWiener_PerTrial_IgnoreArtifactTime',...
			'TESA_DDWiener_PerTrial_BaselineOnly'}
		
		c_say('Detecting bad channels based on DDWiener noise estimates');
		
		if isempty(s.threshold)
			s.threshold = 20;
			c_saySingle('Using default rejection threshold of %s', c_toString(s.threshold));
		else
			c_saySingle('Using specified rejection threshold of %s', c_toString(s.threshold));
		end
		
		switch(s.detectionMethod)
			case {'TESA_DDWiener', 'TESA_DDWiener_IgnoreArtifactTime'}
				tmp = mean(tmpEEG.data,3);
				[~, sigmas] = DDWiener(tmp);
			case {'TESA_DDWiener_PerTrial', 'TESA_DDWiener_PerTrial_IgnoreArtifactTime', 'TESA_DDWiener_PerTrial_BaselineOnly'}
				[~, sigmas] = DDWiener(reshape(tmpEEG.data, [tmpEEG.nbchan, tmpEEG.pnts*tmpEEG.trials]));
			otherwise
				error('Not implemented');
		end
		
		misc.channelScores = sigmas; 
		misc.scoreThreshold = median(sigmas) + s.threshold*mad(sigmas,1);
		
		badChannels = sigmas > misc.scoreThreshold;
		
		if true
			% also reject channels with sigma=0 (indicating channel was const throughout)
			badChannels = badChannels | sigmas==0;
		end
		
		badChannels = find(badChannels);
		
		% if s.doPlot
		% 	misc.hf = figure; 
		% 	topoplot(misc.channelScores, EEG.chanlocs,...
		% 		'electrodes', 'on',...
		% 		'emarker2',{badChannels,'x','k',10,2});
		% 	hc = colorbar; 
		% 	caxis([0 misc.scoreThreshold]); 
		% 	ylabel(hc,'Channel noise');
		% end
		
	otherwise
		if iscellstr(s.detectionMethod)
			% processed above, do nothing here
		else
			error('Not implemented');
		end
end

assert(isempty(badChannels) || isnumeric(badChannels)); % should not be logical indices

misc.badChannelIndices = badChannels;

if isempty(badChannels)
	c_sayDone('No bad channels detected');
else
	c_sayDone('%d bad channel%s detected: %s', ...
		length(badChannels),...
		c_strIfNumIsPlural(length(badChannels)),...
		c_toString({EEG.chanlocs(badChannels).labels}));
end

switch(s.replaceMethod)
	case 'interpolate'
		if length(badChannels) == EEG.nbchan
			error('All channels marked as bad, cannot interpolate');
		end
		c_say('Interpolating %d bad channel%s', length(badChannels), c_strIfNumIsPlural(length(badChannels)));
		EEG = eeg_interp(EEG, badChannels);
		c_sayDone();
	case 'NaN'
		if ismember(s.detectionMethod, 'fromASR')
			% bad channels already replaced with NaNs above
		else
			EEG.data(badChannels,: ,:) = NaN;
		end
	case {'remove', 'delete'}
		c_say('Removing %d bad channel%s', length(badChannels), c_strIfNumIsPlural(length(badChannels)));
		fn = @() pop_select(EEG, 'nochannel', badChannels);
		if true
			[~, EEG] = evalc('fn()');
		else
			EEG = fn();
		end
	case 'none'
		% do not change input data
		c_saySingle('Not interpolating, replacing, or removing bad channels.');
	otherwise
		error('Not implemented');		
end

% Determine the index for the new entry
if isfield(EEG, 'badChans') && ~isempty(EEG.badChansDetect)
    newIdx = length(EEG.badChansDetect) + 1;
else
    newIdx = 1;
end

% Save information in the structured array
EEG.badChansDetect(newIdx).indices   = badChannels;                          % numeric indices of bad channels
EEG.badChansDetect(newIdx).labels    = {EEG.chanlocs(badChannels).labels};   % corresponding channel labels
EEG.badChansDetect(newIdx).method    = detectionMethod;                     % method used for detection
EEG.badChansDetect(newIdx).score     = misc.channelScores(badChannels);       % detection scores for each bad channel
EEG.badChansDetect(newIdx).threshold = misc.scoreThreshold;                   % threshold used

%% --- Update running list of all bad channels using labels ---
badLabels = {EEG.chanlocs(badChannels).labels};
if isfield(EEG, 'badChans') && ~isempty(EEG.badChans)
    EEG.badChans = unique([EEG.badChans, badLabels]);
else
    EEG.badChans = badLabels;
end

end

%%
function [y_solved, sigmas] = DDWiener(data)  
% This function computes the data-driven Wiener estimate (DDWiener),
% providing the estimated signals and the estimated noise-amplitudes
%
% .........................................................................
% From TESA toolbox
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

% Compute the sample covariance matrix
C = data*data';

gamma = mean(diag(C));

% Compute the DDWiener estimates in each channel
chanN = size(data,1);
for i=1:chanN
    idiff = setdiff(1:chanN,i);
    y_solved(i,:) = C(i,idiff)*((C(idiff,idiff)+gamma*eye(chanN-1))\data(idiff,:));
end

% Compute the noise estimates in all channels 
sigmas = sqrt(diag((data-y_solved)*(data-y_solved)'))/sqrt(size(data,2));

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
function res = c_isFieldAndNonEmpty(struct,field)
	res = c_isField(struct,field) && ~c_isEmptyOrEmptyStruct(c_getField(struct,field));
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
function varargout = c_saySingle(varargin)
% c_saySingle - Wrapper around fprintf with other added features
%
% See also: c_say

	global sayNestLevel;
	
	if nargout > 0
		varargout{1} = ''; %TODO: possibly change in the future to return meaningful strings
	end

	if isempty(sayNestLevel)
		sayNestLevel = 0;
	end
	
	global saySilenceLevel
	if ~isempty(saySilenceLevel) && sayNestLevel >= saySilenceLevel
		% don't print anything
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
	fprintf('> ');
	fprintf(varargin{:});
	fprintf('\n');
end

%%
function s = c_toString(c,varargin)
% c_toString Convert various data types to a string
% Handles reasonable printing of scalars, vectors, matrices, cells, etc.
%
% Examples:
%	c_toString(rand(2,3))
%	c_toString(rand(2,3),'doPreferMultiline',true)
%	c_saySingleMultiline('%s',c_toString(rand(2,3),'doPreferMultiline',true))
%   c_toString({'test',pi,[1 2 3]}')
%	c_toString({'test',pi,[1 2 3]}','precision',3)

	if nargin == 0, testfn(); s=''; return; end
	
	doPreferMultiline = false;
	doQuoteStrings = true;
	precision = [];
	indentation = 0;
	printLimit = 100;
	
	if nargin > 1 
		p = inputParser();
		p.addParameter('doPreferMultiline',doPreferMultiline,@islogical);
		p.addParameter('doQuoteStrings',doQuoteStrings,@islogical);
		p.addParameter('printLimit',printLimit,@isscalar);
		p.addParameter('precision',precision,@isscalar);
		p.addParameter('indentation',indentation,@isscalar);
		p.parse(varargin{:});
		doPreferMultiline = p.Results.doPreferMultiline;
		doQuoteStrings = p.Results.doQuoteStrings;
		printLimit = p.Results.printLimit;
		precision = p.Results.precision;
		indentation = p.Results.indentation;
	end

	num2strArgs = {};
	if ~isempty(precision)
		num2strArgs = {precision};
	end
	
	if iscell(c)
		s = c_cellToString(c,varargin{:},'indentation',indentation+1);
	elseif isempty(c) && isnumeric(c)
			s = '[]';
	elseif isempty(c) && ischar(c)
		if doQuoteStrings
			s = '''''';
		else
			s = '';
		end
	elseif isscalar(c) && isnumeric(c)
		s = num2str(c,num2strArgs{:});
	elseif isnumeric(c)
		if isvector(c) && length(c) > 2 && all(diff(c)==1) && c_isinteger(c)
			s = ['[' num2str(c(1)) ':' num2str(c(end)) ']'];
			if size(c,1) > size(c,2)
				s = [s,'.'''];
			end
		else
			if numel(c) > printLimit
				%warning('Too many elements to print');
				s = sprintf('<Array of size %s>',c_toString(size(c)));
			else
				if ismatrix(c)
					if size(c,1) > size(c,2) && ~doPreferMultiline
						c = c.';
						didTranspose = true;
					else
						didTranspose = false;
					end
					s = '[';
					
					numStrs = arrayfun(@(x) num2str(x, num2strArgs{:}),c,'UniformOutput',false);
					
					if doPreferMultiline 
						% use less compact format
						% make all strings the same length to line up columns neatly
						maxStrLength = max(cellfun(@length,numStrs(:)));
						templateStr = repmat(' ',1,maxStrLength);
						for i=1:numel(numStrs)
							tmp = numStrs{i};
							numStrs{i} = templateStr;
							numStrs{i}(1:length(tmp)) = tmp;
						end
						withinLineDelim = sprintf('\t');
						betweenLineDelim = sprintf(';\n ');
					else
						withinLineDelim = ' ';
						betweenLineDelim = '; ';
					end
					for i=1:size(c,1)
						s = [s, strjoin(numStrs(i,:),withinLineDelim)];
						if i~=size(c,1)
							s = [s, betweenLineDelim];
						end
					end
					s = [s,']'];
					if didTranspose
						s = [s,'.'''];
					end
				else
					s = sprintf('<Array of size %s>',c_toString(size(c)));
					%TODO: add support for showing values of arrays of higher dimensions
				end
			end
		end
	elseif ischar(c)
		if numel(c) > length(c)
			s = sprintf('<Char array of size %s>',c_toString(size(c)));
		else
			if doQuoteStrings
				s = ['''' c ''''];
			else
				s = c;
			end
		end
	elseif islogical(c)
		s = num2str(c);
		s = strrep(s,'0','false');
		s = strrep(s,'1','true');
	elseif isstruct(c)
		if length(c) > 1
			s = sprintf('<Struct array>:\n\t');
			if length(c) < printLimit
				for i=1:length(c)
					s = [s, sprintf('<Element %d>:\n\t',i)];
					s = [s, indentLines(c_toString(c(i),varargin{:}))];
					if i~=length(c)
						s = [s, sprintf('\n')];
					end
				end
			else
				s = [s, '<too long to print>'];
			end
		else
			fields = fieldnames(c);
			longestFieldLength = max(cellfun(@length,fields));
			if length(c) > 0
				s = sprintf('<Struct>:\n');
			else
				s = sprintf('<Empty struct>:\n');
			end
			for i=1:length(fields)
				if i~=length(fields)
					s = [s,' |'];
				% else
				% 	s = [s,' |'];
				end
				s = [s, repmat('_',1,longestFieldLength - length(fields{i}) + 1) ' '];
				s = [s, sprintf('%s',fields{i})];
				if length(c)>0
					s = [s, ': ', indentLines(c_toString(c.(fields{i}),varargin{:}))];
				end
				if i~=length(fields)
					s = [s, sprintf('\n')];
				end
			end
		end
	elseif isdatetime(c)
		s = strtrim(evalc('disp(c)'));
	elseif iscategorical(c)
		if length(c)==1
			s = c_toString(char(c),varargin{:});
		else
			s = c_toString(cellstr(c),varargin{:});
		end
	elseif isobject(c)
		s = sprintf('[Object]:');
		tmp = indentLines([sprintf('\n') evalc('disp(c)')]);
		if length(tmp(:)) > printLimit*10
			tmp = '<too long to print>';
		end
		s = [s, tmp];
	elseif isa(c,'function_handle')
		s = sprintf('[function_handle]: %s',char(c));
	else
		error('unsupported type');
	end
end

function s = indentLines(s)
	s = strrep(s,sprintf('\n'),sprintf('\n\t'));
end

function testfn()
	c_toString({'test',[1 2 3],'a1',{'test inner', 5}});
	c_toString(struct('parent1','child1','parent2',struct('child2',1,'child3','subchild')));
end

%%
function str = c_strIfNumIsPlural(varargin)
p = inputParser();
p.addRequired('num',@isscalar);
p.addOptional('strIfPlural','s',@ischar);
p.addParameter('elseStr','',@ischar);
p.parse(varargin{:});
s = p.Results;

if abs(s.num)~=1
	str = s.strIfPlural;
else
	str = s.elseStr;
end
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
function isEmpty = c_isEmptyOrEmptyStruct(struct)
	isEmpty = isempty(struct) || ((isstruct(struct) || isobject(struct)) && isempty(fieldnames(struct)));
end

%%
function varargout = c_getField(struct,fieldName)
% c_getField - like getfield(), but allows nested field names
%
% Example:
%	a_struct = struct('outer',struct('inner',1))
%	c_getField(a_struct,'outer.inner')

	assert(ischar(fieldName))
	i = find(fieldName=='.',1,'first');
	if isempty(i)
		[varargout{1:nargout}] = struct.(fieldName);
	else
		% recursive call
		if length(struct)==1
			[varargout{1:nargout}] = c_getField(struct.(fieldName(1:i-1)),fieldName(i+1:end));
		else
			varargout = cell(1,length(struct));
			for iS = 1:length(struct)
				varargout{iS} = c_getField(struct(iS).(fieldName(1:i-1)),fieldName(i+1:end));
			end
		end
	end
end

%%
function s = c_setField(s,fieldName,value)

if nargout == 0 && ~isa(s,'handle')
	warning('Must assign output of %s to store changes to struct',mfilename);
end

if length(s) > 1
	for i=1:numel(s)
		s(i) = c_setField(s(i),fieldName,value);
	end
	return;
end
assert(ischar(fieldName));
i = find(fieldName=='.',1,'first');
if isempty(i)
	s.(fieldName) = value;
else
	% recursive call
	if ~isfield(s, (fieldName(1:i-1)))
		s.(fieldName(1:i-1)) = struct();
	end
	s.(fieldName(1:i-1)) = c_setField(s.(fieldName(1:i-1)),fieldName(i+1:end),value);
end

end

%%
function out = c_isinteger(x)
	out = isinteger(x) || (isnumeric(x) && all(ceil(x(:)) == floor(x(:))));
end

%%
function s = c_cellToString(c,varargin)
	if nargin == 0
		% test
		c_cellToString({'test',[1 2 3],'a1',{'test inner', 5}})
		return
	end
	
	doPreferMultiline = false;
	doQuoteStrings = true;
	precision = [];
	indentation = 0;
	
	if nargin > 1
		p = inputParser();
		p.addParameter('doPreferMultiline',doPreferMultiline,@islogical);
		p.addParameter('doQuoteStrings',doQuoteStrings,@islogical);
		p.addParameter('printLimit',[],@isscalar);
		p.addParameter('precision',precision,@isscalar);
		p.addParameter('indentation',indentation,@isscalar);
		p.parse(varargin{:});
		doPreferMultiline = p.Results.doPreferMultiline;
		doQuoteStrings = p.Results.doQuoteStrings;
		precision = p.Results.precision;
		indentation = p.Results.indentation;
	end

	assert(iscell(c));
	
	if isempty(c)
		s = '{}';
		return;
	end
	
	num2strArgs = {};
	if ~isempty(precision)
		num2strArgs = {precision};
	end
	
	s = '{';
	if doPreferMultiline
		s = [s sprintf('\t')];
	end
	assert(ndims(c)==2);
	for i=1:size(c,1)
		for j=1:size(c,2)
			if iscell(c{i,j})
				s = [s c_cellToString(c{i,j},varargin{:},'indentation',indentation+1) ','];
			elseif ischar(c{i,j})
				if doQuoteStrings
					s = [s '''' c{i,j} '''' ','];
				else
					s = [s c{i,j} ','];
				end
			else
				s = [s c_toString(c{i,j},varargin{:}) ','];
			end
		end
		s = s(1:end-1); % remove comma
		if i ~= size(c,1)
			s = [s ';'];
			if doPreferMultiline
				s = [s,sprintf('\n ')];
				for ii = 1:indentation
					s = [s,sprintf('\t')];
				end
			end
		end
	end
	s = [s '}'];
end