% pop_tesa_fitartifactmodel() - this script fits a model to the capacitive
%                           discharge artifact following TMS and then 
%                           subtracts it from the data as described in:
%
% Freche D, Naim-Feil J, Peled A, Levit- Binnun N, Moses E (2018) A 
% quantitative physical model of the TMS-induced discharge artifacts in
% EEG. PLoS Comput Biol 14(7): e1006177. 
% https://doi.org/10.1371/journal.pcbi.1006177
%
% Please cite this article if you use this function.
%
% Usage:
%   >>  EEG = pop_tesa_fitartifactmodel( EEG ); % pop up window
%   >>  EEG = pop_tesa_fitartifactmodel( EEG, 'key1',value1... );
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs (varargin):
%   'tPulseMs', int   - Time point of the TMS pulse (ms)
%                   Default: 0
%   'tPulseDurationMs', int - Duration of the pulse artifact (ms)
%                   Default: 2
%   'tSkipMs', int - Time to be skipped after tPulseDurationMs (in ms)
%                   Default: 8
%   'tSelectMs', int - Time point used to select the two traces following 
%                   tPulseMs (in ms).
%                   Default: 12
%   't_fitmodel_p_ms', [start,end] - Time ranges for model fitting, 
%                   positive traces (ms). Note that more that one trace can
%                   be definded.
%                   Default: [0,60];
%                   Example: [0,20;0,20]; provide time range for two
%                   traces.
%   't_fitmodel_n_ms',[start,end] - Time ranges for model fitting, 
%                   negative traces (ms). Note that more that one trace can
%                   be definded.
%                   Default: [0,60];
%                   Example: [0,20;0,20]; provide time range for two
%                   traces.
%   't_lincomb_ms', [start,end] - Time range where the linear coefficients 
%                   are to be found (ms)
%   'epoch_range', [start,end] - Epoch over which to apply the correction (ms).
%                   Leave empty to apply over whole epoch.
%                   Default: [-500,500];
%                   Example: [];
% Outputs:
%   EEG             - EEGLAB EEG structure
% 
% Examples:
%   EEG = pop_tesa_fitartifactmodel(EEG); % default use with pop up window
%   EEG = pop_tesa_fitartifactmodel(EEG, 'tPulseDurationMs', 5 ); % Increase the TMS pulse duration that is excluded to 5 ms.

% This script was adapted by Nigel Rogasch for the TESA toolbox. Original
% code is available from:
% https://osf.io/q3vjd/

% Copyright Weizmann Institute of Science (2018)

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, com] = pop_tesa_fitartifactmodel( EEG, varargin )

com = ''; 

    if nargin < 1
        error('Not enough input arguments.');
    end
    
    %check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

    %define defaults
    options = struct('tPulseMs',0,'tPulseDurationMs',2, 'tSkipMs',8,'tSelectMs',12,'t_fitmodel_p_ms',[0,60],'t_fitmodel_n_ms',[0,60],'t_lincomb_ms',[0,100],'epoch_range',[-500,500]);
    
    % read the acceptable names
    optionNames = fieldnames(options);

    % count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('tesa_compselect needs key/value pairs')
    end

    for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
       inpName = pair{1}; % make case insensitive

       if any(strcmpi(inpName,optionNames))%looks for known options and replaces these in options
          options.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end


% Inputs (defaults)
tPulseMs = options.tPulseMs; % time point of the TMS pulse (ms)
tPulseDurationMs = options.tPulseDurationMs; % duration of the pulse artifact (ms)
tSkipMs = options.tSkipMs; % time to be skipped after tPulseDuration (in ms)
tSelectMs = options.tSelectMs; % time point used to select the two traces (indexed by I) following tPulse (in ms)
t_fitmodel_p_ms = options.t_fitmodel_p_ms; % Time ranges for model fitting, positive traces (ms)
t_fitmodel_n_ms = options.t_fitmodel_n_ms; % time ranges for model fitting, negative traces (ms)
t_lincomb_ms = options.t_lincomb_ms; % time range where the linear coefficients for V(t,i) are to be found (in ms)
epoch_range = options.epoch_range; % epoch over which to apply the correction (in ms; leave empty to apply over whole epoch)

% Generate approriate strings
if size(t_fitmodel_p_ms,1) >1
    temp = num2str(t_fitmodel_p_ms);
    for sx = 1:size(temp,1)-1
        temp(sx,size(temp,2)+1) = ';';
        temp(sx,size(temp,2)+1) = ' ';
    end
    tempStrP = reshape(temp',1,[]);
else
    tempStrP = num2str(t_fitmodel_p_ms);
end

if size(t_fitmodel_n_ms,1) >1
    temp = num2str(t_fitmodel_n_ms);
    for sx = 1:size(temp,1)-1
        temp(sx,size(temp,2)+1) = ';';
        temp(sx,size(temp,2)+1) = ' ';
    end
    tempStrN = reshape(temp',1,[]);
else
    tempStrN = num2str(t_fitmodel_n_ms);
end

% pop up window
% -------------
if nargin < 2
    
geometry = {1 ... % Title
            [1 0.5] ... % Time point TMS pulse
            [1 0.5] ... % Duration pulse artifact
            [1 0.5] ... % Duration skip
            [1 0.5] ... % Time point for +/- traces
            [1 0.5] ... % Time range for + traces
            [1 0.5] ... % Time range for - traces
            [1 0.5] ... % Time range for linear combination
            [1 0.5] ... % Epoch over which to apply correction     
        };
    
uilist = {{'style', 'text', 'string', 'Fit model to decay artifact and subtract from data','fontweight','bold'} ...
    {'style', 'text', 'string', 'Time point of TMS pulse (ms)'} ...
    {'style', 'edit', 'string', num2str(tPulseMs)} ...
    {'style', 'text', 'string', 'Duration of the pulse artifact (ms)'} ...
    {'style', 'edit', 'string', num2str(tPulseDurationMs)} ...
    {'style', 'text', 'string', 'Time range to be skipped after pulse artifact (ms)'} ...
    {'style', 'edit', 'string', num2str(tSkipMs)} ...
    {'style', 'text', 'string', 'Time point used to select the two traces following TMS pulse (ms)'} ...
    {'style', 'edit', 'string', num2str(tSelectMs)} ...
    {'style', 'text', 'string', 'Time ranges for model fitting, positive traces (ms; to fit multiple traces use: t1 t2; t3 t4)'} ...
    {'style', 'edit', 'string', tempStrP} ...
    {'style', 'text', 'string', 'Time ranges for model fitting, negative traces (ms; to fit multiple traces use: t1 t2; t3 t4)'} ...
    {'style', 'edit', 'string', tempStrN} ...
    {'style', 'text', 'string', 'Time range where the linear coefficients are to be found (ms)'} ...
    {'style', 'edit', 'string', num2str(t_lincomb_ms)} ...
    {'style', 'text', 'string', 'Epoch over which to apply the correction (ms; leave empty to apply over whole epoch)'} ...
    {'style', 'edit', 'string', num2str(epoch_range)} ...
    };

    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Fit and subtract decay artifact -- pop_tesa_fitartifactmodel()', 'helpcom', 'pophelp(''pop_tesa_fitartifactmodel'')');
    if isempty(result), return; end
    
    % Time point TMS pulse
    tPulseMs = str2num(result{1});
    
    % Duration pulse artifact
    tPulseDurationMs = str2num(result{2});
    
    % Duration skip
    tSkipMs = str2num(result{3});
    
    % Time point for +/- traces
    tSelectMs = str2num(result{4});
    
    % Time range for + traces
    t_fitmodel_p_ms = str2num(result{5});
    
    % Time range for - traces
    t_fitmodel_n_ms = str2num(result{6});

    % Time range for linear combination
    t_lincomb_ms = str2num(result{7});

    % Epoch over which to apply correction   
    if strcmp(result{8},'')
        epoch_range = [];
    else
        epoch_range = str2num(result{8}); 
    end
    
end

% Epoch data for interval required for correction
if ~isempty(epoch_range)
    EEG = pop_epoch( EEG, {  }, [epoch_range(1)./1000 epoch_range(2)./1000], 'newname', 'Merged datasets epochs', 'epochinfo', 'yes');
end

% Populate some defaults from the data (in samples)
samplerate = EEG.srate; % in Hz
[~,tPulse] = min(abs(EEG.times-tPulseMs)); % sample point index marking TMS pulse application in EEGdata
tPulseDuration = round(tPulseDurationMs.*(samplerate./1000)); % duration of the pulse artifact (samples)
tSkip   = round(tSkipMs.*(samplerate./1000)); % number of samples to be skipped after tPulseDuration (in samples)
tSelect = round(tSelectMs.*(samplerate./1000)); % sample point used to select the two traces (indexed by I) following tPulse (in samples)

% t_fitmodel_p{i}     time range for trace i in I
t_fitmodel_p = []; % cell array of time ranges for model fitting.
for px = 1:size(t_fitmodel_p_ms,1)
    t_fitmodel_p{px} = round(t_fitmodel_p_ms(px,1).*(samplerate./1000)):round(t_fitmodel_p_ms(px,2).*(samplerate./1000));
end

% t_fitmodel_n{i}     time range for trace i in I
t_fitmodel_n = []; % cell array of time ranges for model fitting.
for nx = 1:size(t_fitmodel_n_ms,1)
    t_fitmodel_n{nx} = round(t_fitmodel_n_ms(nx,1).*(samplerate./1000)):round(t_fitmodel_n_ms(nx,2).*(samplerate./1000));
end

t_lincomb = round(t_lincomb_ms(1).*(samplerate./1000)):round(t_lincomb_ms(2).*(samplerate./1000));

% Epoching settings (in samples)
tPrePulse  = EEG.xmin*samplerate:-1;
tSkipped   = 1:tSkip-1;

% Calculate the time periods for each section of data
t1 = tPrePulse;
t2 = 0:tPulseDuration;
t3 = tPulseDuration + tSkipped;
tInt = [t1,t2,t3];
tPostPulse = 0:(length(EEG.times) - length(tInt))-1;
t4 = tPulseDuration + tSkip + tPostPulse;


for iT = 1:size(EEG.data,3)
    
    fprintf('Fitting %d of %d\n', iT, size(EEG.data,3));

    electrodeLst = 1:size(EEG.data,1);

    EEGdata = double( EEG.data(electrodeLst,:,iT) );

    %%% *** uncomment the following line *** to subtract common mode signal (subtract average of all electrodes)
    % EEGdata = EEGdata - mean(EEGdata,1);

    [ V, u, PrepulseOfs ] = tesa_fitartifactmodel(EEGdata, tPulse, tPulseDuration, tSkip, tSelect, t_fitmodel_p, t_fitmodel_n, t_lincomb, samplerate);

    A = V(tPostPulse,1:size(u,1)); % reconstructed artifacts for time range tPostPulse

    y1 = [];
    y2 = [];
    y3 = [];
    y4 = [];
    mOut = [];
    
    for iE = 1:size(EEGdata,1)
        % pre-pulse data
        y1(iE,:) = EEGdata(iE,tPulse + t1);
        y2(iE,:) = EEGdata(iE,tPulse + t2);
        y3(iE,:) = EEGdata(iE,tPulse + t3);
        % subtract reconstructed artifact
        y4(iE,:) = EEGdata(iE,tPulse + t4)' - A*u(:,iE);
        mOut(iE,:) = A*u(:,iE);
    end
    
    dataOut(:,:,iT) = [y1,y2,y3,y4];
    modelOut(:,:,iT) = mOut;
    modelTime = tPulse+t4;

end  

% Replace data with corrected data
EEG.data = dataOut;
EEG.modelOutput = modelOut;
EEG.modelTime = modelTime;

%Run script from input
if nargin <2
    com = sprintf('%s = pop_tesa_fitartifactmodel( %s, ''tPulseMs'', %s, ''tPulseDurationMs'', %s, ''tSkipMs'', %s, ''tSelectMs'', %s, ''t_fitmodel_p_ms'', [%s], ''t_fitmodel_n_ms'', [%s], ''t_lincomb_ms'', [%s], ''epoch_range'', [%s] );', inputname(1), inputname(1),result{1}, result{2}, result{3}, result{4}, result{5}, result{6}, result{7}, result{8});
elseif nargin > 1
    com = sprintf('%s = pop_tesa_fitartifactmodel( %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end