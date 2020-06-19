% tesa_sspsir()-      Uses SSP-SIR method to suppress TMS-evoked muscle artifacts [1] OR
%                     control data [2] from the provided EEGLAB dataset
%
%                     [1] Mutanen, T. P., Kukkonen, M., Nieminen, J. O., Stenroos, M., Sarvas, J.,
%                     & Ilmoniemi, R. J. (2016). Recovering TMS-evoked EEG responses masked by
%                     muscle artifacts. Neuroimage, 139, 157-166.
%
%                     [2] Biabani, M, Fornito, A, Mutanen, T. P., Morrow, J, & Rogasch, N. C.(2019).
%                     Characterizing and minimizing the contribution of sensory inputs to TMS-evoked
%                     potentials.Brain Stimulation, 12(6):1537-1552.
%                     
%
%
% Usage:
%  >>  [EEG] = tesa_sspsir(EEG); % run SSP-SIR using default values
%  >>  [EEG] = tesa_sspsir(EEG, 'key1',value1... ); % run SSP-SIR using customised 
%      inputs
%
% Optional input pairs (varargin)
%
%
% 'leadfieldIn',      - Lead field matrix : channel x vertices - If not provided, a simple spherical                     
%                      lead field for the EEGLAB structure will be constructed.
%
% 'leadfieldChansFile'- Leadfield channel matrix (channelName x 1). If provided, sorts the leadfield
%                     channel order to match EEGLAB data (if required).
%
%
% 'artScale',         - Choose the method to estimate the artefacts in the PCA :
%                     'automatic' (default): Uses a sliding window to scale the signal relative to 
%                     the amplitude
%                     'manual': Only uses the data within a window provided by the user in timeRange
%                     'manualConstant': Works like 'manual' but projects out the artifact dimensions 
%                     estimated from timeRange across the whole time domain. Useful for the validation 
%                     of SSPSIR.
%                     'control': The only option to remove control data. Uses the  data from the 
%                     control data as provided in EEG_control.
%                     NOTE: 
%                     Setting artScale to 'automatic' only cleans the detected time window of
%                     muscle artifacts.
%                     Setting the artScale to 'manual' only cleans the specified time window of muscle 
%                     artifacts. Setting the artScale to 'manualConstant' estimates muscle artifacts 
%                     from the given time window (e.g. [0 50]), but cleans the dataset indentically                     
%                     across the whole time domain. Option 'manualConstant' is particularly useful for 
%                     controlling how much the neuronal signals of interest might be attenuated by 
%                     SSPSIR.
%
% 'timeRange',        - Required for artScale :'manual' and 'manualConstant'. Vector with start and 
%                     end times in ms of window containing TMS-evoked muscle response.
%                     - Optional for artScale :'control'. Vector with start and end times in ms of 
%                     window containing control responses.
%                     Note that multiple time windows can be defined.
%                     Example: [5,50] - one window
%                     Example: [5,50;100,120] - multiple windows
%                     Multiple windows are useful e.g. when comparing pre- and post datasets and 
%                     attempting to clean the datasets identically.
%                     In such a case the datasets should be concatenated prior to TESA_SSPSIR and 
%                     the time windows should correspond the artifactual window of both of the 
%                     compared datasets.
%
%
% 'PC',               - The number of artifact PCs to be removed. If not provided plots the PC 
%                     dimentions to choose from (default).
%                     [N]: N = The number of principal components to remove
%                     {'data', [N]} : Uses the provided EEG data to find the number of components
%                     which explain N% of the variations
%                     Example:  {'data', [90]} - removes the components which explain more than 90%
%                     of the variance in the EEG.data
%
%
% 'M',                - The truncation dimension used in the SIR step. If not provided estimates it 
%                     from the EEG data as> rank(data) - PC
%
%
% 'EEG_control'       - EEG data (*.set) from control condition. Required for artScale = 'control'.
%                     [fileName,filePath] Example: ['/Users/myPC/Desktop/controlResponse.set/']
%                     NOTE: EEG_control is expected to be provided in an evoked form.
%
% 
% Outputs:
% EEG                 - EEGLAB EEG structure ( Output of SSPSIR applied on single trials )
%                     NOTE: 
%                     EEG.meanTrials  is the output of SSPSIR applied to the average of all trials
%                     saved on a new field   
%
% Examples:
%  >> [EEG] = tesa_sspsir( EEG ); % default use
%  >> [EEG] = tesa_sspsir( EEG, 'artScale', 'manual','timeRange',[5,50], 'PC',...  
%     {'data', [90]} ); Suppresses muscle artefacts by removing the data components explaining more 
%     than 90% of variance in the time winodw of 5-50ms 
%  >> [EEG] = tesa_sspsir( EEG, 'artScale', 'control','PC',  [5], 'EEG_control',...
%     ['/Users/myPC/Desktop/controlResponse.set/']); % Suppresses control data by removing the first
%     5 principal components of controlResponse.Data from EEG.data 
%
%
% Copyright (C) 2019 
% Mana Biabani Mana.Biabanimoghadam@monash.edu, Monash University,
% Tuomas Mutanen tuomas.mutanen@aalto.fi, Aalto University
% Nigel Rogasch,nigel.rogasch@adelaide.edu.au University of Adelaide,
% Jukka Sarvas Jukka.Sarvas@aalto.fi, Aalto University
%
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
% 30/09/2019


function [EEG] = tesa_sspsir(EEG,varargin)
% ------------------------------------ Check EEG inputs ----------------------------------------
% Check that EEG channels have been correctly specified
if nargin <1
    error('Not enough inputs. Please provide the EEG structure as an input.');
end

for i=1:length(EEG.chanlocs)
    if isempty(EEG.chanlocs(i).type)
        EEG.chanlocs(i).type = 'EEG';
    end
end

% ------------------------------------ set defaults --------------------------------------------
% Define defaults
options = struct('leadfieldIn',[],'leadfieldChansFile',[],'artScale',[],'timeRange',[],'PC',[],'M',[],'EEG_control',[]);

% Read the acceptable names
optionNames = fieldnames(options);

% Count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('tesa_ssp_sir needs key/value pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
    inpName = pair{1}; % make case insensitive
    
    if any(strcmpi(inpName,optionNames))%looks for known options and replaces these in options
        options.(inpName) = pair{2};
    else
        error('%s is not a recognized parameter name',inpName)
    end
end

if isempty(options.leadfieldIn)
    leadfieldIn = [];
else
    options.leadfieldIn =  struct2array(load(options.leadfieldIn)); %options.leadfieldIn;    
end

if isempty(options.leadfieldChansFile)
    leadfieldChansFile = [];
else
    leadfieldChansFile = struct2array(load(options.leadfieldChansFile));
    options.leadfieldChansFile = leadfieldChansFile;
end

if isempty(options.artScale)
    artScale = 'automatic';
else
    artScale = options.artScale;
end

if isempty(options.timeRange)
    timeRange = [];
else
    timeRange = options.timeRange;
end

if isempty(options.PC)
    PC = [];
else
    PC = options.PC;
end

if isempty(options.M)
    M = [];
else
    M = options.M;
end

if isempty(options.EEG_control)
    EEG_control = [];
else
    EEG_control = pop_loadset(options.EEG_control);
end

% ----------------------------------- Check leadfield inputs -----------------------------------

if ~isempty(options.leadfieldIn) && isempty(options.leadfieldChansFile)
    warning('leadfieldChansFile is not provided. Please make sure that the EEG recordings and the individualized head model have the same order of channels!')
end

% If lead-field is not specified, construct a simple spherical lead field for the EEGLAB structure.
if isempty(options.leadfieldIn)
    % Build the spherical lead field, using the theoretical electrode-locations of the data.
    fprintf('No leadfield provided. Calculating based on spherical model\n');
    leadfieldIn = construct_spherical_lead_field(EEG);
elseif ~isempty(options.leadfieldIn)
    leadfieldIn = options.leadfieldIn;
    fprintf('Using customised leadfield\n');
    size(leadfieldIn,1);
    if size(EEG.data,1) ~= size(leadfieldIn,1) % If not replacing channels
        error('Number of channels in data (%d) does not match number of channels in the leadfield matrix (%d).', size(EEG.data,1), size(leadfieldIn,1));
    end
end

% Sort the leadfield channel order to match EEGLAB data (if required)
if ~isempty(options.leadfieldChansFile)
    leadfieldChans = options.leadfieldChansFile;
    eeglabChans = {EEG.allchan.labels}';    
    if size(leadfieldIn,1) ~= length(leadfieldChans) % Check the size of the lead field matrix matches the number of channels.
        error('Number of channels in lead field matrix (%d) does not match number of channels in the lead field channel file (%d)\n', size(leadfieldIn,1), length(leadfieldChans));
    end
    
    for i = 1:size(eeglabChans,1)
        [~,chanIndex(i,1)] = ismember(lower(eeglabChans{i}),lower(leadfieldChans));
    end
    leadfieldIn = leadfieldIn(chanIndex,:);
    fprintf('Lead field matrix channel dimension sorted to match EEGLAB data.\n');
end

% Check that EEG channels have been correctly specified and trying to correct if necessary.

    all_chan_labels = {'Fp1','EEG';'Fpz','EEG';'Fp2','EEG';'AF9','EEG';'AF7','EEG';'AF5','EEG';'AF3','EEG';'AF1','EEG';'AFz','EEG';'AF2','EEG';'AF4','EEG';'AF6','EEG';'AF8','EEG';'AF10','EEG';'F9','EEG';'F7','EEG';'F5','EEG';'F3','EEG';'F1','EEG';'Fz','EEG';'F2','EEG';'F4','EEG';'F6','EEG';'F8','EEG';'F10','EEG';'FT9','EEG';'FT7','EEG';'FC5','EEG';'FC3','EEG';'FC1','EEG';'FCz','EEG';'FC2','EEG';'FC4','EEG';'FC6','EEG';'FT8','EEG';'FT10','EEG';'T9','EEG';'T7','EEG';'C5','EEG';'C3','EEG';'C1','EEG';'Cz','EEG';'C2','EEG';'C4','EEG';'C6','EEG';'T8','EEG';'T10','EEG';'TP9','EEG';'TP7','EEG';'CP5','EEG';'CP3','EEG';'CP1','EEG';'CPz','EEG';'CP2','EEG';'CP4','EEG';'CP6','EEG';'TP8','EEG';'TP10','EEG';'P9','EEG';'P7','EEG';'P5','EEG';'P3','EEG';'P1','EEG';'Pz','EEG';'P2','EEG';'P4','EEG';'P6','EEG';'P8','EEG';'P10','EEG';'PO9','EEG';'PO7','EEG';'PO5','EEG';'PO3','EEG';'PO1','EEG';'POz','EEG';'PO2','EEG';'PO4','EEG';'PO6','EEG';'PO8','EEG';'PO10','EEG';'O1','EEG';'Oz','EEG';'O2','EEG';'I1','EEG';'O9','EEG';'Iz','EEG';'I2','EEG';'O10','EEG';'AFp9h','EEG';'AFp7h','EEG';'AFp5h','EEG';'AFp3h','EEG';'AFp1h','EEG';'AFp2h','EEG';'AFp4h','EEG';'AFp6h','EEG';'AFp8h','EEG';'AFp10h','EEG';'AFF9h','EEG';'AFF7h','EEG';'AFF5h','EEG';'AFF3h','EEG';'AFF1h','EEG';'AFF2h','EEG';'AFF4h','EEG';'AFF6h','EEG';'AFF8h','EEG';'AFF10h','EEG';'FFT9h','EEG';'FFT7h','EEG';'FFC5h','EEG';'FFC3h','EEG';'FFC1h','EEG';'FFC2h','EEG';'FFC4h','EEG';'FFC6h','EEG';'FFT8h','EEG';'FFT10h','EEG';'FTT9h','EEG';'FTT7h','EEG';'FCC5h','EEG';'FCC3h','EEG';'FCC1h','EEG';'FCC2h','EEG';'FCC4h','EEG';'FCC6h','EEG';'FTT8h','EEG';'FTT10h','EEG';'TTP9h','EEG';'TTP7h','EEG';'CCP5h','EEG';'CCP3h','EEG';'CCP1h','EEG';'CCP2h','EEG';'CCP4h','EEG';'CCP6h','EEG';'TTP8h','EEG';'TTP10h','EEG';'TPP9h','EEG';'TPP7h','EEG';'CPP5h','EEG';'CPP3h','EEG';'CPP1h','EEG';'CPP2h','EEG';'CPP4h','EEG';'CPP6h','EEG';'TPP8h','EEG';'TPP10h','EEG';'PPO9h','EEG';'PPO7h','EEG';'PPO5h','EEG';'PPO3h','EEG';'PPO1h','EEG';'PPO2h','EEG';'PPO4h','EEG';'PPO6h','EEG';'PPO8h','EEG';'PPO10h','EEG';'POO9h','EEG';'POO7h','EEG';'POO5h','EEG';'POO3h','EEG';'POO1h','EEG';'POO2h','EEG';'POO4h','EEG';'POO6h','EEG';'POO8h','EEG';'POO10h','EEG';'OI1h','EEG';'OI2h','EEG';'Fp1h','EEG';'Fp2h','EEG';'AF9h','EEG';'AF7h','EEG';'AF5h','EEG';'AF3h','EEG';'AF1h','EEG';'AF2h','EEG';'AF4h','EEG';'AF6h','EEG';'AF8h','EEG';'AF10h','EEG';'F9h','EEG';'F7h','EEG';'F5h','EEG';'F3h','EEG';'F1h','EEG';'F2h','EEG';'F4h','EEG';'F6h','EEG';'F8h','EEG';'F10h','EEG';'FT9h','EEG';'FT7h','EEG';'FC5h','EEG';'FC3h','EEG';'FC1h','EEG';'FC2h','EEG';'FC4h','EEG';'FC6h','EEG';'FT8h','EEG';'FT10h','EEG';'T9h','EEG';'T7h','EEG';'C5h','EEG';'C3h','EEG';'C1h','EEG';'C2h','EEG';'C4h','EEG';'C6h','EEG';'T8h','EEG';'T10h','EEG';'TP9h','EEG';'TP7h','EEG';'CP5h','EEG';'CP3h','EEG';'CP1h','EEG';'CP2h','EEG';'CP4h','EEG';'CP6h','EEG';'TP8h','EEG';'TP10h','EEG';'P9h','EEG';'P7h','EEG';'P5h','EEG';'P3h','EEG';'P1h','EEG';'P2h','EEG';'P4h','EEG';'P6h','EEG';'P8h','EEG';'P10h','EEG';'PO9h','EEG';'PO7h','EEG';'PO5h','EEG';'PO3h','EEG';'PO1h','EEG';'PO2h','EEG';'PO4h','EEG';'PO6h','EEG';'PO8h','EEG';'PO10h','EEG';'O1h','EEG';'O2h','EEG';'I1h','EEG';'I2h','EEG';'AFp9','EEG';'AFp7','EEG';'AFp5','EEG';'AFp3','EEG';'AFp1','EEG';'AFpz','EEG';'AFp2','EEG';'AFp4','EEG';'AFp6','EEG';'AFp8','EEG';'AFp10','EEG';'AFF9','EEG';'AFF7','EEG';'AFF5','EEG';'AFF3','EEG';'AFF1','EEG';'AFFz','EEG';'AFF2','EEG';'AFF4','EEG';'AFF6','EEG';'AFF8','EEG';'AFF10','EEG';'FFT9','EEG';'FFT7','EEG';'FFC5','EEG';'FFC3','EEG';'FFC1','EEG';'FFCz','EEG';'FFC2','EEG';'FFC4','EEG';'FFC6','EEG';'FFT8','EEG';'FFT10','EEG';'FTT9','EEG';'FTT7','EEG';'FCC5','EEG';'FCC3','EEG';'FCC1','EEG';'FCCz','EEG';'FCC2','EEG';'FCC4','EEG';'FCC6','EEG';'FTT8','EEG';'FTT10','EEG';'TTP9','EEG';'TTP7','EEG';'CCP5','EEG';'CCP3','EEG';'CCP1','EEG';'CCPz','EEG';'CCP2','EEG';'CCP4','EEG';'CCP6','EEG';'TTP8','EEG';'TTP10','EEG';'TPP9','EEG';'TPP7','EEG';'CPP5','EEG';'CPP3','EEG';'CPP1','EEG';'CPPz','EEG';'CPP2','EEG';'CPP4','EEG';'CPP6','EEG';'TPP8','EEG';'TPP10','EEG';'PPO9','EEG';'PPO7','EEG';'PPO5','EEG';'PPO3','EEG';'PPO1','EEG';'PPOz','EEG';'PPO2','EEG';'PPO4','EEG';'PPO6','EEG';'PPO8','EEG';'PPO10','EEG';'POO9','EEG';'POO7','EEG';'POO5','EEG';'POO3','EEG';'POO1','EEG';'POOz','EEG';'POO2','EEG';'POO4','EEG';'POO6','EEG';'POO8','EEG';'POO10','EEG';'OI1','EEG';'OIz','EEG';'OI2','EEG';'T3','EEG';'T5','EEG';'T4','EEG';'T6','EEG';'M1','EEG';'M2','EEG';'A1','EEG';'A2','EEG'};

    channel_type_warning_made = 0;
    wrong_type_channels = {};
    

    for i=1:length(EEG.chanlocs)
        if ~strcmp(EEG.chanlocs(i).type,'EEG')
            if ~channel_type_warning_made
                warning('All channel types are not specified! Attempting to detect and classify EEG channels based on labels.' )
                channel_type_warning_made = 1;
            end
            if isempty(find(ismember(all_chan_labels, EEG.chanlocs(i).labels))) && isempty(find(ismember(lower(all_chan_labels), EEG.chanlocs(i).labels))) && isempty(find(ismember(upper(all_chan_labels), EEG.chanlocs(i).labels))) && isempty(find(ismember(upper(all_chan_labels), upper(EEG.chanlocs(i).labels))))
                    warning(['The data contain non-EEG channels! Rejecting the ',EEG.chanlocs(i).labels,' channel from further analysis.']);
                    wrong_type_channels{end + 1} = EEG.chanlocs(i).labels;
            else    
                    EEG.chanlocs(i).type = 'EEG';
           end
        end
    end


if ~isempty(wrong_type_channels)
 
    if (strcmp(options.artScale,'control'))
        EEG_control =  pop_select( EEG_control,'nochannel',wrong_type_channels);  
    end
 
    EEG = pop_select( EEG,'nochannel',wrong_type_channels);

end
clear channel_type_warning_made;
clear wrong_type_channels;

%--------------------------------------------------------------------------

% Check that data in average reference:
if ~strcmp(EEG.ref,'averef')
    warning('The data is not in average reference. Note that tesa_sspsir returns the data in the average reference.')
    EEG = pop_reref( EEG, []);
end
    
% Calculate the average of trials
% (if the average has been already computed this does not change it)
EEG_evo = mean(EEG.data,3);

% ----------------------------------- Check leadfield inputs -----------------------------------

% When artScale = 'control' check if control data�is provided to be removed
if strcmp(options.artScale,'control') && isempty(EEG_control)
    uiwait( msgbox('Please provide EEG_control data as an input.'));
    [fileName,filePath]= uigetfile;
    EEG_control = pop_loadset([filePath,fileName]);
end

% Re-reference the data and the lead field to the channel average:
[L_ave] = ref_ave(leadfieldIn);
[EEG_evo_ave] = ref_ave(EEG_evo);

% Check if EEG channels have been correctly specified in control condition
if (strcmp(options.artScale,'control'))
    for i=1:length(EEG_control.chanlocs)
        if isempty(EEG_control.chanlocs(i).type)
            EEG_control.chanlocs(i).type = 'EEG';
        end
    end
    % Re-reference the control data to the channel average:
    EEG_control_evo = mean(EEG_control.data,3);
    [EEG_control_evo] = ref_ave(EEG_control_evo);
end

% -------------------------------------- apply SSP-SIR -----------------------------------------

if (strcmp(options.artScale,'control'))
    
    if ~isempty(options.timeRange)
        
        warning('Cleaning only the time window based on timeRange. NOTE: Not likely to work correctly if the timeaxis of the control and main data are not identical!')
        [data_correct, artifact_topographies, ~,filt_ker] = SSP_SIR_control(EEG_evo_ave, L_ave, EEG.times,EEG_control_evo, PC, M, timeRange);
    else
        [data_correct, artifact_topographies, ~] = SSP_SIR_control(EEG_evo_ave, L_ave, EEG.times,EEG_control_evo, PC, M);
        filt_ker = [];
    end
    
else
    
    [data_correct, artifact_topographies, ~, filt_ker] = SSP_SIR(EEG_evo_ave, L_ave, EEG.srate, EEG.times, artScale, timeRange, PC, M);
    
    if  isempty(filt_ker)
        filt_ker = -1*sigmf(EEG.times,[0.05*(timeRange(2)-timeRange(1)) timeRange(2)])+ sigmf(EEG.times,[0.05*(timeRange(2)-timeRange(1)) timeRange(1)]);
    end
end

EEG.meanTrials = data_correct;
% -------------------------------------- Return single trials ----------------------------------

% Clean each trial separately with the obtained artifact topograhies.

%EEG = pop_reref( EEG, []);

[EEG_trials] = SSP_SIR_trials(EEG, L_ave, artifact_topographies, filt_ker, M);
EEG.data =  EEG_trials.data;
end

function [data_correct, artifact_topographies, data_clean, filt_ker] = SSP_SIR(data, L, Fs, time, artScale, timeRange, PC,M)

if nargin < 4
    artScale = 'automatic';
end

% High-pass filtering the data from 100 Hz:
[b,a] = butter(2,100/(Fs/2),'high');
data_high = (filtfilt(b,a,double(data')))';

if strcmp(artScale,'automatic')
    
    % Estimating the relative-amplitude kernel of the artifact:
    tmp = data_high.^2;
    % x_scal = 73; %estimating the RMS value in 50-ms sliding windows
    x_scal = round(Fs./1000.*50); %estimating the RMS value in 50-ms sliding windows
    x = int8(linspace(-x_scal/2, x_scal/2,x_scal));
    
    for i=1:size(tmp,1)
        filt_ker(i,:) = conv(tmp(i,:),ones(1,x_scal),'same')/x_scal;
    end
    
    filt_ker = sum(filt_ker,1)/size(tmp,1);
    filt_ker = filt_ker/max(filt_ker);
    filt_ker = sqrt(filt_ker);
    filt_ker = repmat(filt_ker,[size(data_high,1),1]);
    
    % Estimating the artifact-suppression matrix P:
    [U,singular_spectum,~] = svd(filt_ker.*data_high,'econ');
    
elseif strcmp(artScale,'manual') || strcmp(artScale,'manualConstant')
    
    if isempty(timeRange)
        error('Please include a time range for the TMS-evoked muscle response');
    end
    
    % RUN PCA JUST ON DATA CONTAINING THE TMS-EVOKED MUSCLE RESPONSE
    data_short = [];
    for tidx = 1:size(timeRange,1)
        [~,t1] = min(abs(time-timeRange(tidx,1)));
        [~,t2] = min(abs(time-timeRange(tidx,2)));
        data_short = [data_short, data_high(:,t1:t2)];
    end
    
    % Estimating the artifact-suppression matrix P:
    [U,singular_spectum,~] = svd(data_short,'econ');
    
    % Calculating the smooth filtering function for later combination of
    % the projected and non-projected data.
    for tidx = 1:size(timeRange,1)
        
        smoothLength =  10; % The step function transition from 10% to 90% takes roughly 10 ms
        smooth_weighting_function(tidx,:) = dsigmf(time,[4/(smoothLength) timeRange(tidx,1) 4/(smoothLength) timeRange(tidx,2)]);
        
    end
    % Combining different moments of artifacts:
    smooth_weighting_function = sum(smooth_weighting_function,1);
    
    % Copying the same time scaling to all EEG channels:
    smooth_weighting_function = repmat(smooth_weighting_function,[size(data_high,1),1]);
end

if nargin < 6 || isempty(PC)

 % Test (-100 to 500 ms)
    fig1_1 = figure;
    fig1_2 = figure;
    for pcdx = 1:10
        figure(fig1_1);
        subplot(2,5,pcdx)
        datatmp = U(:,pcdx)'*data;
        plot(time,datatmp,'k','LineWidth',1);
        if time(1) <= -20 && time(end) >= 100
            set(gca,'xlim',[-20,100]);
        end
        title(['PC ',num2str(pcdx)]);
        
        xlabel('Time [ms]');
        
        if pcdx == 1 || pcdx ==6
            ylabel('V [\mu V]')
        end
        
        [tf, freqs, times] = timefreq(datatmp, Fs,'tlimits' , [time(1) time(end)]);
        
        figure(fig1_2);
        subplot(2,5,pcdx)
        imagesc(times, freqs, 20*log10(abs(tf)));
        set(gca,'YDir','normal')
        
        if times(1) <= -500 && times(end) > 500
            set(gca,'xlim',[-500, 500]);
        end
        
       if freqs(end) > 200
           set(gca,'ylim',[0, 200]);
           title(['PC ',num2str(pcdx)]);
        end
        
        
        xlabel('Time [ms]');
        
        if pcdx == 1 || pcdx ==6
            ylabel('Freq. [Hz]')
        end
        
    end

        fig2 = figure;
        plot(diag(singular_spectum)/singular_spectum(1,1),'b','LineWidth',2)
        title('Click line to choose the artifact dimension, then press Return.')
        xlabel('PC');
        ylabel('Signal magnitude > 100 Hz (Normalized to 1st PC)');
        datacursormode on;
        dcm_obj = datacursormode(fig2);
        % Wait while the user does this.
        pause
        c_info = getCursorInfo(dcm_obj);
        PC = c_info.Position(1);
        close(fig2);
close(fig1_1);
close(fig1_2);
    elseif iscell(PC)

        if     nargin > 5 && any(strcmp(PC{1},'data')) && size(PC,2)<2
            error('Please provide the percentage of variation that should be covered by PC in the PC cell.')

        elseif nargin > 5 && strcmp(PC{1}, 'data') && size(PC,2)==2
            PC = PC{2}/100;
            i = 1;
            d = diag(singular_spectum);
            a = sum(d(1:i)).^2/sum(d(1:end)).^2;
            while a < PC
                i = i+1;
                a = sum(d(1:i)).^2/sum(d(1:end)).^2;
            end
            PC = i;
        end
        
    end


if nargin < 7 || isempty(M)
    M = rank(data) - PC;
end

P = eye(size(data,1)) - U(:,1:PC)*(U(:,1:PC))';
artifact_topographies = U(:,1:PC);

%Suppressing the artifacts:
data_clean = P*data;

%Performing SIR for the suppressed data:
PL = P*L;
tau_proj = PL*PL';
[U,S,V] = svd(tau_proj);
S_inv = zeros(size(S));
S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
tau_inv = V*S_inv*U';
suppr_data_SIR = L*(PL)'*tau_inv*data_clean;

%Performing SIR for the original data:
tau_proj = L*L';
[U,S,V] = svd(tau_proj);
S_inv = zeros(size(S));
S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
tau_inv = V*S_inv*U';
orig_data_SIR = L*(L)'*tau_inv*data;

if strcmp(artScale,'automatic')
    data_correct = filt_ker.*suppr_data_SIR + orig_data_SIR - filt_ker.*orig_data_SIR;
    filt_ker = filt_ker(1,:);
    
    
elseif strcmp(artScale,'manual')
    %     data_correct = smooth_weighting_function.*suppr_data_SIR + (1 - smooth_weighting_function).*orig_data_SIR;
    data_correct = smooth_weighting_function.*suppr_data_SIR + orig_data_SIR - smooth_weighting_function.*orig_data_SIR;
    
    filt_ker =  smooth_weighting_function(1,:);
elseif strcmp(artScale,'manualConstant')
    data_correct = suppr_data_SIR;
    
    filt_ker =   ones(1,size( data_correct,2));
    
end

end
function [EEG_out] = SSP_SIR_trials(EEG_in, L, art_topographies, filt_ker, M )

%This function cleans each trial separately with a given artifact
%topographies

EEG_out = EEG_in;
P = eye(size(EEG_in.data,1)) - art_topographies*art_topographies';

for i = 1:size(EEG_in.data,3)
    
    data = EEG_in.data(:,:,i);
    
    %Suppressing the artifacts:
    data_clean = P*data;
    
    %Performing SIR for the suppressed data:
    PL = P*L;
    
    if isempty (M)
        M = rank(data_clean) -  1 ;
    end
    
    tau_proj = PL*PL';
    [U,S,V] = svd(tau_proj);
    S_inv = zeros(size(S));
    S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
    tau_inv = V*S_inv*U';
    suppr_data_SIR = L*(PL)'*tau_inv*data_clean;
    
    %Performing SIR for the original data:
    tau_proj = L*L';
    [U,S,V] = svd(tau_proj);
    S_inv = zeros(size(S));
    S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
    tau_inv = V*S_inv*U';
    orig_data_SIR = L*(L)'*tau_inv*data;
    
    if isempty(filt_ker)
        data_correct = suppr_data_SIR;
    else
        filt_ker_B = repmat(filt_ker,[size(suppr_data_SIR,1),1]);
        data_correct = filt_ker_B.*suppr_data_SIR + orig_data_SIR - filt_ker_B.*orig_data_SIR;
        filt_ker_B = [];
    end
    
    EEG_out.data(:,:,i) = data_correct;
    
end
end
function [data_correct, artifact_topographies, data_clean, filt_ker] = SSP_SIR_control(data, leadfieldIn, time, controlData, PC, M, timeRange)

if nargin < 7 || isempty (timeRange)
    
    % Estimating the artifact-suppression matrix P
    [U,singular_spectum,~] = svd(controlData,'econ');
    filt_ker = [];
else
    
    [~,t1] = min(abs(time-timeRange(1)));
    [~,t2] = min(abs(time-timeRange(2)));
    data_control_short = controlData(:,t1:t2);
    
    % Estimating the artifact-suppression matrix P:
    [U,singular_spectum,~] = svd(data_control_short,'econ');
    
    % Calculate scaling function for manual data
    smoothLength = 10; % Length of smooting window either side of 'includeLength'
    filt_ker = dsigmf(time,[4/(smoothLength) timeRange(1) 4/(smoothLength) timeRange(2)]);
    
end

   
if nargin < 6 || isempty(PC)

 fig1_1 = figure;
    fig1_2 = figure;
    for pcdx = 1:10
        figure(fig1_1);
        subplot(2,5,pcdx)
        datatmp = U(:,pcdx)'*data;
        plot(time,datatmp,'k','LineWidth',1);
        if time(1) <= -20 && time(end) >= 100
            set(gca,'xlim',[-20,100]);
        end
        title(['PC ',num2str(pcdx)]);
        
        xlabel('Time [ms]');
        
        if pcdx == 1 || pcdx ==6
            ylabel('V [\mu V]')
        end
        
        Fs = 1000/(time(2) - time(1))
        [tf, freqs, times] = timefreq(datatmp, Fs,'tlimits' , [time(1) time(end)]);
        
        figure(fig1_2);
        subplot(2,5,pcdx)
        imagesc(times, freqs, 20*log10(abs(tf)));
        set(gca,'YDir','normal')
        
        if times(1) <= -500 && times(end) > 500
            set(gca,'xlim',[-500, 500]);
        end
        
       if freqs(end) > 200
           set(gca,'ylim',[0, 200]);
           title(['PC ',num2str(pcdx)]);
        end
        
        xlabel('Time [ms]');
        
        if pcdx == 1 || pcdx ==6
            ylabel('Freq. [Hz]')
        end
        
    end

        fig2 = figure;
        plot(diag(singular_spectum)/singular_spectum(1,1),'b','LineWidth',2)
        title('Click line to choose the artifact dimension, then press Return.')
        xlabel('PC');
        ylabel('Signal magnitude > 100 Hz (Normalized to 1st PC)');
        datacursormode on;
        dcm_obj = datacursormode(fig2);
        % Wait while the user does this.
        pause
        c_info = getCursorInfo(dcm_obj);
        PC = c_info.Position(1);
        close(fig2);
        close(fig1_1);
        close(fig1_2);


elseif iscell(PC)
    if nargin > 4 && any(strcmp(PC{1},'data')) && size(PC,2)<2
        error('Please provide the percentage of variation that should be covered by PC in the PC cell.')
    elseif nargin > 4 && strcmp(PC{1}, 'data') && size(PC,2)==2
        PC = PC{2}/100;
        i = 1;
        d = diag(singular_spectum);
        a = sum(d(1:i)).^2/sum(d(1:end)).^2;
        while a < PC
            i = i+1;
            a = sum(d(1:i)).^2/sum(d(1:end)).^2;
        end
        PC = i;
    end
    
end

% Default M
if nargin < 6 || isempty(M)
    M = rank(data) - PC;
end
P = eye(size(data,1)) - U(:,1:PC)*(U(:,1:PC))';
artifact_topographies = U(:,1:PC);

%Suppressing the artifacts:
data_clean = P*data;

%Performing SIR for the suppressed data:
PL = P*leadfieldIn;
tau_proj = PL*PL';
[U,S,V] = svd(tau_proj);
S_inv = zeros(size(S));
S_inv(1:M,1:M) = diag(1./diag(S(1:M,1:M)));
tau_inv = V*S_inv*U';
suppr_data_SIR = leadfieldIn*(PL)'*tau_inv*data_clean;
data_correct = suppr_data_SIR;
end

function [LFM_sphere] = construct_spherical_lead_field (EEG,r0,r1,r2,r3,sig1,sig2,sig3,DipN)
% This functions creates a spherical lead-field that can be used in SIR
% to correct the data.
%
% Input variable:
%
% EEG = The studied EEGLAB dataset. The function expects that EEG contains
% the chanlocs on a spherical surface.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university
% .........................................................................

% r0 = distance of the current dipoles from the origin in [meters]
% (default, if input r0=[], r0=76*1e-3, as in the original article)
if nargin < 2 || isempty(r0)
    r0 = 76*1e-3;
end

% r1 = Radius of the inner skull surface in [meters] (default, if input
% r1=[], r1=81*1e-3, as in the original article)
if nargin < 3 || isempty(r1)
    r1 = 81*1e-3;
end

% r2 = Radius of the outer skull surface in [meters] (default, if input
% r2=[], r2=85*1e-3, as in the original article)
if nargin < 4 || isempty(r2)
    r2 = 85*1e-3;
end

% r3 = Radius of the scalp surface in [meters] (default, if input
% r3=[], r3=88*1e-3, as in the original article)
if nargin < 5 || isempty(r3)
    r3 = 88*1e-3;
end

% sig1 = Conductivity of th brain in [1/0mega*1/m] (default, if input
% sig1=[], sig1=0.33, as in the original article)
if nargin < 6 || isempty(sig1)
    sig1 = 0.33;
end

% sig2 = Conductivity of the skull in [1/0mega*1/m] (default, if input
% sig2=[], sig2=0.33/50, as in the original article)
if nargin < 7 || isempty(sig2)
    sig2 = 0.33/50;
end

% sig3 = Conductivity of th scalp in [1/0mega*1/m] (default, if input
% sig3=[], sig3=0.33, as in the original article)
if nargin < 8 || isempty(sig3)
    sig3 = 0.33;
end

% DipN = The number of the brain noise source dipoles (default, in input
% DipN=[], DipN=5000, as in the original article)
if nargin < 9 || isempty(DipN)
    DipN = 5000;
end


% Reading the electrode locations from the chanlocs field in the given EEG
% file:

k = 1;
for i=1:length(EEG.chanlocs)
    if strcmp(EEG.chanlocs(i).type,'EEG')
        elec_coords(k,1) = EEG.chanlocs(i).X;
        elec_coords(k,2) = EEG.chanlocs(i).Y;
        elec_coords(k,3) = EEG.chanlocs(i).Z;
        k = k+1;
    end
end

% Projecting the coordinates on a unit sphere:
elec_coords = elec_coords./repmat(sqrt(sum(elec_coords.^2,2)),[1,3]);

%%%% Computing the spherical lead field using function leadfield1

% THE SPHERICAL HEAD MODEL
rad=[r1,r2,r3];
sig=[sig1,sig2,sig3];

% the radius of the cortex hemisphere on which the dipoles lie
R=r3*elec_coords; % electrode positions,

P=randn(3,DipN);
P=r0*P./(ones(3,1)*sqrt(sum(P.^2)));
Pns=[P(1,:);P(2,:);(P(3,:))]; % dipole positions on the upper hemisphere
% with radius r0,
%Moments of the unit dipoles
Q = Pns;
Qns = Q./repmat(sqrt(sum(Q.^2,1)),[3,1]);

%Computing the leadfield accoring to [1]
LFM_sphere = leadfield1(R',Pns,Qns,rad,sig,30);

%[1] Mutanen, Tuomas P., et al. "Recovering TMS-evoked EEG responses masked by muscle artifacts." Neuroimage 139 (2016): 157-166.

end
function V=leadfield1(X,Y,Q,rad,sig,nmax)
% In a layered sphere the program returns the potential V, at  points
% given in X and lying on the outer surface of the sphere, due
% to  current dipoles  located at points in Y with moments in Q;
% points in Y must lie in the innermost  region of  the sphere.
%
% INPUT: X = (3,m) matrix of field points on the outer spherical surface,
%        Y = (3,n) matrix of  the locations of the dipoles,
%        Q = (3,n) matrix of the of moments of the dipoles,
%        rad = the vector of the radii of the shperical regions,
%            rad(1) < rad(2) < ... < rad(M), M=length(r),
%        sig = the vector of the conductivities sig(j) in the j:th
%              region for j=1:M,
%        nmax = the number of terms in the multipole expansion of the
%               potential; if norm(Y) is not very close to r(1), nmax=30 is
%               usually O.K.
%
% OUTPUT: V = (m,n) matrix containing the potentials: V(j,k) = potential
%              at X(:,k) due to dipole (Y(:,j),Q(:,j)), k=1:m, j = 1:n,
% .........................................................................
% 15 August 2013 : Jukka Sarvas, BECS, Aalto university
% .........................................................................

m=size(X,2);
n=size(Y,2);
Xhat=X./(ones(3,1)*sqrt(sum(X.^2)));
normY=sqrt(sum(Y.^2));
Yhat=Y./(ones(3,1)*normY);

alpha=ones(m,1)*sum(Q.*Yhat);
beta=Xhat'*Q;
gamma=Xhat'*Yhat;
[L,dL]=legen(nmax,gamma(:));
P=L(2:nmax+1,:);
dP=dL(2:nmax+1,:);
gam=potgam(rad,sig,nmax);

a=0;
b=0;
for k=1:nmax
    c=gam(k)*(ones(m,1)*(normY.^(k-1)));
    a=a+k*c.*reshape(P(k,:)',m,n);
    b=b+c.*reshape(dP(k,:)',m,n);
end
V=1/(4*pi*sig(end))*(alpha.*a+(beta-alpha.*gamma).*b);

end
function [L,dL]=legen(N,x)
% The function returns the values of the Legendre polynomial P_n(x) of
% order n at x in L(n+1,:) for n=0:N .
% If nargout=2, also the derivatives P'_n(x) are returned in dL(n+1,:).
% Note the the legendre functions P_n_1(x) are obtained from P'_n(x) as
%            P_n_1(x)=-sqrt(1-x^2)*P'_n(x)

L=ones(N+1,length(x));
x=x(:).';
if N>0
    L(2,:)=x;
end
if N<2
    L=L(1:N+1,:);
else
    for n=1:N-1
        L(n+2,:)=(2*n+1)/(n+1)*x.*L(n+1,:)-n/(n+1)*L(n,:);
    end
end

if nargout==2
    dL=zeros(N+1,length(x));
    dL(2,:)=ones(1,length(x));
    for n=1:N-1
        dL(n+2,:)=dL(n,:)+(2*n+1)*L(n+1,:);
    end
end
end
function gam=potgam(r,sig,nmax)
% The function returns the translation coefficients gam(n),
% INPUT: r = vector of (outer) radii of the layers,
%        sig = vector of the conductivities of the layers,
%        nmax = the maximal n index of the multipole expansion,
% OUTPUT: gam = the column vector of the translation coefficients,

M=length(r);
gam=zeros(nmax,1);
s=sig;
for n=1:nmax
    C=eye(2);
    for j=1:M-1
        c=1/((2*n+1)*s(j+1));
        Cj=c*[n*s(j)+(n+1)*s(j+1),(n+1)*(s(j+1)-s(j))/r(j)^(2*n+1);...
            n*r(j)^(2*n+1)*(s(j+1)-s(j)),(n+1)*s(j)+n*s(j+1)];
        C=Cj*C;
    end
    a=(2*n+1)/(n+1)*r(M)^n;
    b=-C(2,1)+n/(n+1)*r(M)^(2*n+1)*C(1,1);
    gam(n)=a/b;
end
end
function [data_ave] = ref_ave(data)
% This function re-references the data to the average reference.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university
% .........................................................................

data_ave = data - repmat(mean(data,1),[size(data,1),1]);

end

