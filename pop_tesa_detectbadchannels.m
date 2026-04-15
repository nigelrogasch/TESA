function [EEG, com] = pop_tesa_detectbadchannels(EEG, varargin)
% pop_tesa_detectbadchannels() - GUI wrapper for tesa_detectbadchannels
%
% Description:
%   Detects bad EEG channels using various TESA detection methods and optionally
%   replaces, interpolates, or removes them. This function wraps the core
%   tesa_detectbadchannels() function in an EEGLAB-style GUI, allowing users
%   to select detection parameters interactively or call the function
%   programmatically.
%
% Inputs:
%   EEG              - EEGLAB EEG structure
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
%   EEG  - Updated EEGLAB EEG structure. Additional fields added:
%          EEG.chansAll - Original EEG.chanlocs saved for later interpolation
%          EEG.badChans - Cell array of channels detected as bad
% 
%   com  - Command string that reproduces the call (for EEG.history)
%
% Usage:
%   % GUI mode
%   EEG = pop_tesa_detectbadchannels(EEG);
%
%   % Command-line mode
%   EEG = pop_tesa_detectbadchannels(EEG, 'detectionMethod', 'PREP_deviation', ...
%                                  'replaceMethod', 'interpolate', ...
%                                  'artifactTimespan', [-2 10], ...
%                                  'threshold', 9);
%
%   % Combine PREP deviation and DDWiener per-trial detection
%   EEG = pop_tesa_detectbadchannels(EEG, ...
%       'detectionMethod', {'PREP_deviation','TESA_DDWiener_PerTrial'}, ...
%       'artifactTimespan', [-2 10]);
%
% Notes:
%   - If multiple detection methods are specified as a cell array, a channel
%     will be marked as bad if identified by any method.
%   - Fields EEG.chansAll and EEG.badChans are automatically added if not present.
%   - This wrapper supports EEGLAB history tracking and works both interactively
%     and from the command line.
%
% See also:
%   tesa_detectbadchannels

com = '';

if nargin < 1
    help pop_tesa_detectbadchannels;
    return;
end

% -------------------------------------------------------------------------
% GUI call
% -------------------------------------------------------------------------
if nargin == 1

    detectionMethods = { ...
        'PREP_deviation', ...
        'TESA_DDWiener', ...
        'TESA_DDWiener_PerTrial', ...
        'TESA_DDWiener_IgnoreArtifactTime', ...
        'TESA_DDWiener_PerTrial_IgnoreArtifactTime', ...
        'TESA_DDWiener_PerTrial_BaselineOnly', ...
        'fromASR' };

    replaceMethods = {'interpolate','remove','NaN','none'};

    uilist = { ...
        {'style','text','string','Detection method','fontweight','bold'} ...
        {'style','popupmenu','string', detectionMethods,'tag','detectionMethod'} ...
        {} ...
        {'style','text','string','Replace method','fontweight','bold'} ...
        {'style','popupmenu','string', replaceMethods,'tag','replaceMethod'} ...
        {} ...
        {'style','text','string','Artifact timespan [start end] (s)'} ...
        {'style','edit','string',''} ...
        {} ...
        {'style','text','string','Detection threshold (optional)'} ...
        {'style','edit','string',''} ...
    };

    geom = { [1] [1] [1] [1] [1] [1] [1] [1] [1] [1] };

    [res, ~, ~, outstruct] = inputgui( ...
        'uilist', uilist, ...
        'geometry', geom, ...
        'title', 'TESA: Detect bad channels');

    if isempty(res)
        return;
    end

    args = {};
    args = [args {'detectionMethod', detectionMethods{outstruct.detectionMethod}}];
    args = [args {'replaceMethod', replaceMethods{outstruct.replaceMethod}}];

    % artifact timespan
    if ~isempty(outstruct.edit1)
        artifactTimespan = str2num(outstruct.edit1); %#ok<ST2NM>
        if numel(artifactTimespan) == 2
            args = [args {'artifactTimespan', artifactTimespan}];
        end
    end

    % threshold
    if ~isempty(outstruct.edit2)
        threshold = str2double(outstruct.edit2);
        if ~isnan(threshold)
            args = [args {'threshold', threshold}];
        end
    end

    % Call underlying function
    [EEG, ~] = tesa_detectbadchannels(EEG, args{:});

    % History
    com = sprintf('EEG = pop_tesa_detectbadchannels(EEG, %s);', ...
        vararg2str(args));
    EEG = eeg_hist(EEG, com);

    return;
end

% -------------------------------------------------------------------------
% Command-line call
% -------------------------------------------------------------------------
[EEG, ~] = tesa_detectbadchannels(EEG, varargin{:});

com = sprintf('EEG = pop_tesa_detectbadchannels(EEG, %s);', ...
    vararg2str(varargin));
EEG = eeg_hist(EEG, com);

end