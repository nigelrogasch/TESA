% tesa_sound() - performs noise suppression on the data using the SOUND algorithm.
%                   Note that data are rereferenced to average during SOUND
%                   For further details on SOUND, or if you use this
%                   algorithm, please cite:
%                   
%                   Mutanen, TP, et al. Automatic and robust noise suppression 
%                   in EEG and MEG: The SOUND algorithm. Neuroimage. 2018 
%                   Feb 1;166:135-151.   
%
% Usage:
%   >>  EEG = tesa_sound( EEG ); % run SOUND using default values and a spherical leadfield
%   >>  EEG = tesa_sound( EEG, 'key1',value1... );% run SOUND using customised inputs
%
% Inputs:
%   EEG             - (required) EEGLAB EEG structure. EEG.data can include single trials, or the average of across trials. 
%                         
% 
% Optional input pairs (varargin):
%
%   'lambdaValue', double       - (optional) double providing the lambda value for SOUND. Default is 0.1.
%                               Example: 0.1
%
%   'iter', int                 - (optional) integer providing the number of iterations for SOUND. Default is 5.
%                               Example: 5
%
%   'leadfieldInVar', matrix    - (optional) n x m matrix with individualised leadfield where n = channels and m = dipoles. 
%                               IMPORTANT - channel number and order must match the channels in the analyzed EEGLAB structure OR 
%                               The channel order should be sorted to match EEGLAB using the 'leadfieldChansFile' input below.
%                               Note: The data will be re-referenced to an average output.
%                               Example: leadfield [e.g. where leadfield = 62 x 15006 matrix (channel x dipole) stored in MATLAB workspace]
%
%   'leadfieldInFile', str      - (optional) file path and name of individualised leadfield matrix where n = channels and m = dipoles.
%                               This command can be used as an alternative to 'leadfieldInVar' which will load the required leadfield matrix in to the workspace.
%                               Note that the file must only contain the leadfield matrix and nothing else.
%                               IMPORTANT - channel number and order must match channels in EEG structure OR the channel
%                               order should be sorted to match EEGLAB using the 'leadfieldChansFile' input below.
%                               Example: 'C:\data\myLeadfield.mat' which contains leadfield variable [e.g. where leadfield = 62 x 15006 matrix (channel x dipole)]
%
%   'leadfieldChansFile', str   - (optional) file path and name of cell array listing the channel order of the individual leadfield matrix. 
%                               If called, this command will load the indicated cell array and then sort the leadfield matrix channel order to match the EEGLAB data channel order.
%                               This is important if the leadfield has been generated in another program which stores the EEG data in a different order.
%                               Note that the channel names must match those in the EEGLAB data.
%                               Example: 'C:\data\myLeadfieldChans.mat' which contains cell array with channel order.
%
%   'multipleConds', cell array - (optional) cell string array indicating the variables of the additional datasets that need to be cleaned simultaneously with input variable EEG.  
%                               Use this input if you would like to apply SOUND identically to multiple conditions (e.g. pre and post an intervention). 
%                               Conditions must first be epoched. Epochs from different conditions identified in the 'multipleConds' input will be separated, 
%                               averaged across epochs within a condition, and then concatenated in time before performing SOUND.
%                               Example: {'EEG2', 'EEG3'}; 
%                               For instance, EEG1 = pop_tesa_sound( EEG1, 'multipleConds', {'EEG2', 'EEG3'} ); would clean simultaneously identically the
%                               datasets EEG1, EEG2, EEG3. Note, only the cleaned version of EEG1 is returned as an output and the datasets
%                               defined within the 'multipleConds' variable are directly saved to the workspace with new names EEG2_after_SOUND and EEG3_after_SOUND.
%                               The names of the additional datasets will be also saved into an extra field in the output variable EEG, i.e., EEG.multipleCondsClean 
%
%   'replaceChans',             - (optional) path and file name for EEGLAB EEG structure containing all required output channels stored in the EEG.chanlocs field.
%                               This command uses SOUND to replace missing channels (e.g. ones that have been removed earlier in the processing pipeline).
%                               If called, this command will replace any channels missing from the current data set relative to the data set called by 'replaceChans'.
%                               Example: 'C:\data\fullChannelEegData.set';
%                               NOTE: This works also with 'multipleConds' option with the restriction that the same channels must have been removed in all
%                               datasets prior to SOUND.
%                               NOTE: When using your own leadfield, it should have channels matching the 'replaceChans'-input dataset. 
%                                
% 
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = tesa_sound( EEG ); % run SOUND using default values and a spherical leadfield
%   EEG = tesa_sound( EEG, 'lambdaValue', 0.2, 'iter', 10 ); % run SOUND using customised input values 
%   EEG = tesa_sound( EEG, 'leadfieldInVar', leadfield, 'leadfieldChansFile', 'C:\data\myLeadfieldChans.mat' ); % run SOUND using default values and an individual leadfield matrix (stored in variable called leadfield). Sort the leadfield channel order (stored in 'C:\data\myLeadfieldChans.mat') to match the EEGLAB data channel order.
%   EEG = tesa_sound( EEG, 'leadfieldInFile', 'C:\data\myLeadfield.mat' ); % run SOUND using default values and an individual leadfield matrix (stored in variable called leadfield loaded from 'C:\data\myLeadfield.mat')
%   EEG1 = pop_tesa_sound( EEG1, 'multipleConds', {'EEG2', 'EEG3'} ); % Clean simultaneously identically the  datasets EEG1, EEG2, EEG3. 
%   EEG = tesa_sound( EEG, 'replaceChans', 'C:\data\fullChannelEegData.set' ); % run SOUND using default values and replacing missing channels in the current data set (all channels defined in the loaded file). 
% 
% See also:
%   tesa_sspsir 

% Copyright (C) 2020 Tuomas Mutanen, Aalto University, Tuomas.Mutanen@gmail.com   
% Nigel Rogasch, Monash University, nigel.rogasch@monash.edu
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

function [EEG, ERP_cleaned, LFM_sphere_mean] = tesa_sound(EEG, varargin)

% Check that there is enough inputs
if nargin <1
    error('Not enough inputs. Please provide the EEG structure as an input');
end

% Define defaults
options = struct('lambdaValue',[],'iter',[],'leadfieldInVar',[],'leadfieldInFile',[],'multipleConds',[],'leadfieldChansFile',[],'replaceChans',[]);

% Read the acceptable names
optionNames = fieldnames(options);

% Count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('tesa_sound needs key/value pairs')
end

for pair = reshape(varargin,2,[]) % pair is {propName;propValue}
   inpName = pair{1}; % make case insensitive

   if any(strcmpi(inpName,optionNames))%looks for known options and replaces these in options
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

%Set defaults
if isempty(options.lambdaValue)
    lambdaValue = 0.1;
else
    lambdaValue = options.lambdaValue;
end
if isempty(options.iter)
    iter = 5;
else
    iter = options.iter;
end
if isempty(options.leadfieldInVar)
    leadfieldIn = [];
else
    leadfieldIn = options.leadfieldInVar;
end
if isempty(options.leadfieldInFile)
    leadfieldInFile = [];
else
    leadfieldInFile = options.leadfieldInFile;
end
if isempty(options.leadfieldChansFile)
    leadfieldChansFile = [];
else
    leadfieldChansFile = options.leadfieldChansFile;
end
if isempty(options.multipleConds)
    multipleConds = [];
else
    multipleConds = options.multipleConds;
end
if isempty(options.replaceChans)
    replaceChans = [];
else
    replaceChans = options.replaceChans;
end

% Check iter inputs
if mod(iter,1) ~= 0
    error('Input ''iter'' needs to be an integer');
end

% Check that data in average reference:
if ~strcmp(EEG.ref,'averef')
    warning('The data is not in average reference. Note that tesa_SOUND returns the data in the average reference.')
end

% Check if channels are being replaced and load data if so
if ~isempty(replaceChans)
    EEGallChans = pop_loadset(replaceChans);
    allChan = {EEGallChans.chanlocs.labels};
    EEG.allChan = allChan;
    currentChan = {EEG.chanlocs.labels};
    EEG.goodC = ismember(allChan,currentChan);
    EEG.badC = EEG.goodC == 0;
    EEG.badCName = allChan(EEG.badC);
    fprintf('%d channel(s) are being replaced by SOUND\n',sum(EEG.badC));
end

% Check leadfield inputs from Var and File
if ~isempty(leadfieldIn) && ~isempty(leadfieldInFile)
    warning('A leadfield has been specified from both the workspace and a file. Using the leadfield from file.');
end

% Load the leadfield (if required)
if ~isempty(leadfieldInFile)
    leadfieldInLoad = load(leadfieldInFile);
    fieldsLoad = fieldnames(leadfieldInLoad);
    if length(fieldsLoad) > 1
        error('Lead field file must only contain a single matrix with channels x dipoles');
    end
    leadfieldIn = leadfieldInLoad.(fieldsLoad{1});
    
end

% Check leadfield inputs
if ~isempty(leadfieldIn)
    fprintf('Using customised leadfield\n');
    if ~isempty(replaceChans)
        if size(EEGallChans.data,1) ~= size(leadfieldIn) % If replacing channels
            error('Number of channels in data (all channels) (%d) does not match number of channels in the leadfield matrix (%d).', size(EEGallChans.data,1), size(leadfieldIn,1));
        end  
    else
        if size(EEG.data,1) ~= size(leadfieldIn) % If not replacing channels
            error('Number of channels in data (%d) does not match number of channels in the leadfield matrix (%d).', size(EEG.data,1), size(leadfieldIn,1));
        end
    end
else
    fprintf('No leadfield provided. Calculating based on spherical model\n');
end

% Sort the leadfield channel order to match EEGLAB data (if required)
if ~isempty(leadfieldChansFile)
    leadfieldChansIn = load(leadfieldChansFile);
    fieldsIn = fieldnames(leadfieldChansIn);
    if length(fieldsIn) > 1
        error('Lead field channel file must only contain a single cell array with electrode names');
    end
    leadfieldChans = leadfieldChansIn.(fieldsIn{1});
    if ~isempty(replaceChans) % if replacing channels
        eeglabChans = {EEGallChans.chanlocs.labels};
    else % if not replacing channels
        eeglabChans = {EEG.chanlocs.labels};
    end
    
    if size(leadfieldIn,1) ~= length(leadfieldChans) % Check the size of the lead field matrix matches the number of channels.
        error('Number of channels in lead field matrix (%d) does not match number of channels in the lead field channel file (%d)\n', size(leadfieldIn,1), length(leadfieldChans));
    end
    
    for i = 1:size(eeglabChans,2)
        [~,chanIndex(i)] = ismember(lower(eeglabChans{i}),lower(leadfieldChans));
    end
    leadfieldIn = leadfieldIn(chanIndex,:);
    fprintf('Lead field matrix channel dimension sorted to match EEGLAB data.\n');
end

    % Check that EEG channels have been correctly specified and trying to correct if necessary.

    all_chan_labels = {'Fp1','EEG';'Fpz','EEG';'Fp2','EEG';'AF9','EEG';'AF7','EEG';'AF5','EEG';'AF3','EEG';'AF1','EEG';'AFz','EEG';'AF2','EEG';'AF4','EEG';'AF6','EEG';'AF8','EEG';'AF10','EEG';'F9','EEG';'F7','EEG';'F5','EEG';'F3','EEG';'F1','EEG';'Fz','EEG';'F2','EEG';'F4','EEG';'F6','EEG';'F8','EEG';'F10','EEG';'FT9','EEG';'FT7','EEG';'FC5','EEG';'FC3','EEG';'FC1','EEG';'FCz','EEG';'FC2','EEG';'FC4','EEG';'FC6','EEG';'FT8','EEG';'FT10','EEG';'T9','EEG';'T7','EEG';'C5','EEG';'C3','EEG';'C1','EEG';'Cz','EEG';'C2','EEG';'C4','EEG';'C6','EEG';'T8','EEG';'T10','EEG';'TP9','EEG';'TP7','EEG';'CP5','EEG';'CP3','EEG';'CP1','EEG';'CPz','EEG';'CP2','EEG';'CP4','EEG';'CP6','EEG';'TP8','EEG';'TP10','EEG';'P9','EEG';'P7','EEG';'P5','EEG';'P3','EEG';'P1','EEG';'Pz','EEG';'P2','EEG';'P4','EEG';'P6','EEG';'P8','EEG';'P10','EEG';'PO9','EEG';'PO7','EEG';'PO5','EEG';'PO3','EEG';'PO1','EEG';'POz','EEG';'PO2','EEG';'PO4','EEG';'PO6','EEG';'PO8','EEG';'PO10','EEG';'O1','EEG';'Oz','EEG';'O2','EEG';'I1','EEG';'O9','EEG';'Iz','EEG';'I2','EEG';'O10','EEG';'AFp9h','EEG';'AFp7h','EEG';'AFp5h','EEG';'AFp3h','EEG';'AFp1h','EEG';'AFp2h','EEG';'AFp4h','EEG';'AFp6h','EEG';'AFp8h','EEG';'AFp10h','EEG';'AFF9h','EEG';'AFF7h','EEG';'AFF5h','EEG';'AFF3h','EEG';'AFF1h','EEG';'AFF2h','EEG';'AFF4h','EEG';'AFF6h','EEG';'AFF8h','EEG';'AFF10h','EEG';'FFT9h','EEG';'FFT7h','EEG';'FFC5h','EEG';'FFC3h','EEG';'FFC1h','EEG';'FFC2h','EEG';'FFC4h','EEG';'FFC6h','EEG';'FFT8h','EEG';'FFT10h','EEG';'FTT9h','EEG';'FTT7h','EEG';'FCC5h','EEG';'FCC3h','EEG';'FCC1h','EEG';'FCC2h','EEG';'FCC4h','EEG';'FCC6h','EEG';'FTT8h','EEG';'FTT10h','EEG';'TTP9h','EEG';'TTP7h','EEG';'CCP5h','EEG';'CCP3h','EEG';'CCP1h','EEG';'CCP2h','EEG';'CCP4h','EEG';'CCP6h','EEG';'TTP8h','EEG';'TTP10h','EEG';'TPP9h','EEG';'TPP7h','EEG';'CPP5h','EEG';'CPP3h','EEG';'CPP1h','EEG';'CPP2h','EEG';'CPP4h','EEG';'CPP6h','EEG';'TPP8h','EEG';'TPP10h','EEG';'PPO9h','EEG';'PPO7h','EEG';'PPO5h','EEG';'PPO3h','EEG';'PPO1h','EEG';'PPO2h','EEG';'PPO4h','EEG';'PPO6h','EEG';'PPO8h','EEG';'PPO10h','EEG';'POO9h','EEG';'POO7h','EEG';'POO5h','EEG';'POO3h','EEG';'POO1h','EEG';'POO2h','EEG';'POO4h','EEG';'POO6h','EEG';'POO8h','EEG';'POO10h','EEG';'OI1h','EEG';'OI2h','EEG';'Fp1h','EEG';'Fp2h','EEG';'AF9h','EEG';'AF7h','EEG';'AF5h','EEG';'AF3h','EEG';'AF1h','EEG';'AF2h','EEG';'AF4h','EEG';'AF6h','EEG';'AF8h','EEG';'AF10h','EEG';'F9h','EEG';'F7h','EEG';'F5h','EEG';'F3h','EEG';'F1h','EEG';'F2h','EEG';'F4h','EEG';'F6h','EEG';'F8h','EEG';'F10h','EEG';'FT9h','EEG';'FT7h','EEG';'FC5h','EEG';'FC3h','EEG';'FC1h','EEG';'FC2h','EEG';'FC4h','EEG';'FC6h','EEG';'FT8h','EEG';'FT10h','EEG';'T9h','EEG';'T7h','EEG';'C5h','EEG';'C3h','EEG';'C1h','EEG';'C2h','EEG';'C4h','EEG';'C6h','EEG';'T8h','EEG';'T10h','EEG';'TP9h','EEG';'TP7h','EEG';'CP5h','EEG';'CP3h','EEG';'CP1h','EEG';'CP2h','EEG';'CP4h','EEG';'CP6h','EEG';'TP8h','EEG';'TP10h','EEG';'P9h','EEG';'P7h','EEG';'P5h','EEG';'P3h','EEG';'P1h','EEG';'P2h','EEG';'P4h','EEG';'P6h','EEG';'P8h','EEG';'P10h','EEG';'PO9h','EEG';'PO7h','EEG';'PO5h','EEG';'PO3h','EEG';'PO1h','EEG';'PO2h','EEG';'PO4h','EEG';'PO6h','EEG';'PO8h','EEG';'PO10h','EEG';'O1h','EEG';'O2h','EEG';'I1h','EEG';'I2h','EEG';'AFp9','EEG';'AFp7','EEG';'AFp5','EEG';'AFp3','EEG';'AFp1','EEG';'AFpz','EEG';'AFp2','EEG';'AFp4','EEG';'AFp6','EEG';'AFp8','EEG';'AFp10','EEG';'AFF9','EEG';'AFF7','EEG';'AFF5','EEG';'AFF3','EEG';'AFF1','EEG';'AFFz','EEG';'AFF2','EEG';'AFF4','EEG';'AFF6','EEG';'AFF8','EEG';'AFF10','EEG';'FFT9','EEG';'FFT7','EEG';'FFC5','EEG';'FFC3','EEG';'FFC1','EEG';'FFCz','EEG';'FFC2','EEG';'FFC4','EEG';'FFC6','EEG';'FFT8','EEG';'FFT10','EEG';'FTT9','EEG';'FTT7','EEG';'FCC5','EEG';'FCC3','EEG';'FCC1','EEG';'FCCz','EEG';'FCC2','EEG';'FCC4','EEG';'FCC6','EEG';'FTT8','EEG';'FTT10','EEG';'TTP9','EEG';'TTP7','EEG';'CCP5','EEG';'CCP3','EEG';'CCP1','EEG';'CCPz','EEG';'CCP2','EEG';'CCP4','EEG';'CCP6','EEG';'TTP8','EEG';'TTP10','EEG';'TPP9','EEG';'TPP7','EEG';'CPP5','EEG';'CPP3','EEG';'CPP1','EEG';'CPPz','EEG';'CPP2','EEG';'CPP4','EEG';'CPP6','EEG';'TPP8','EEG';'TPP10','EEG';'PPO9','EEG';'PPO7','EEG';'PPO5','EEG';'PPO3','EEG';'PPO1','EEG';'PPOz','EEG';'PPO2','EEG';'PPO4','EEG';'PPO6','EEG';'PPO8','EEG';'PPO10','EEG';'POO9','EEG';'POO7','EEG';'POO5','EEG';'POO3','EEG';'POO1','EEG';'POOz','EEG';'POO2','EEG';'POO4','EEG';'POO6','EEG';'POO8','EEG';'POO10','EEG';'OI1','EEG';'OIz','EEG';'OI2','EEG';'T3','EEG';'T5','EEG';'T4','EEG';'T6','EEG';'M1','EEG';'M2','EEG';'A1','EEG';'A2','EEG'};

    channel_type_warning_made = 0;
    wrong_type_channels = {};
    
if ~isempty(replaceChans)
    for i=1:length(EEGallChans.chanlocs)
        if ~strcmp(EEGallChans.chanlocs(i).type,'EEG')
            if ~channel_type_warning_made
                warning('All channel types are not specified! Attempting to detect the EEG channels based on labels.' )
                channel_type_warning_made = 1;
            end
            if isempty(find(ismember(all_chan_labels, EEGallChans.chanlocs(i).labels))) && isempty(find(ismember(lower(all_chan_labels), EEGallChans.chanlocs(i).labels))) && isempty(find(ismember(upper(all_chan_labels), EEGallChans.chanlocs(i).labels))) && isempty(find(ismember(upper(all_chan_labels), upper(EEGallChans.chanlocs(i).labels))))
                    warning(['The reference data with all channels contain non-EEG channels! Rejecting the ',EEGallChans.chanlocs(i).labels,' channel from further analysis.'])
                    wrong_type_channels{end + 1} = EEGallChans.chanlocs(i).labels;    
            else    
                    EEGallChans.chanlocs(i).type = 'EEG';
           end
        end
    end
else
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
end

if ~isempty(wrong_type_channels)
EEG = pop_select( EEG,'nochannel',wrong_type_channels);
end
clear channel_type_warning_made;
clear wrong_type_channels;

% Check whether one or multiple datasets are being used
if ~isempty(multipleConds) % More than one datasets are being contanenated
    for eegx = 1:length(multipleConds)
        tempCond = evalin('base', multipleConds{eegx})
        EEG_evo{eegx} = mean(tempCond.data,3);
    end
    
% Contatenate data sets (if multiple are being used)
    data_cat = [mean(EEG.data,3)];
    
    for lpx = 1:length(EEG_evo)
        data_cat = [data_cat,EEG_evo{lpx}];
    end
    
else % One dataset is being used
    
    % Check if data have been averaged
    if size(EEG.data,3) > 1
        % Estimate the evoked responses as the mean over trials
        EEG_evo = mean(EEG.data,3);
    else
        EEG_evo = EEG.data;
    end
    
    % Contatenate data sets (if multiple are being used)
    data_cat = EEG_evo;
end


% Build the spherical lead field, using the theoretical electrode-locations
% of the data, or an use an individual lead field.

if isempty(leadfieldIn)
    
    if isempty(replaceChans)
        [LFM_sphere] = construct_spherical_lead_field(EEG);
    else
        [LFM_sphere] = construct_spherical_lead_field(EEGallChans); % Calculate based on all final channels, including missing channels
    end
else
    % When using SOUND only to clean the channel signals (like here), and not to
    % improve source estimation per se, we can replace the full lead-field matrix with a square matrix as follows.
    % This enhances the computation significantly.
    
    [U,S,~] = svd(leadfieldIn,'econ');
    LFM_sphere =  U*S; 
    
end

% Re-reference the data and the lead field to the channel with the least
% noise (based on simple Wiener estimate)
[~, sigmas] = DDWiener(data_cat);

[~,bestC] = min(sigmas);
[datatmp] = ref_best(data_cat, bestC);

% Calculate the leadfields for the final reconstruction step of SOUND. 
% The final step will transform the data into average reference.
[LFM_sphere_mean] = ref_ave(LFM_sphere);

% Check if any channels, which will be interpolated in the end, are missing in EEG
% and modify the temporary lead field accordingly to match EEG. 
if ~isempty(replaceChans)
    LFM_sphere_goodC = LFM_sphere(EEG.goodC,:);
    [LFM_sphere_tmp] = ref_best(LFM_sphere_goodC, bestC);
else
    [LFM_sphere_tmp] = ref_best(LFM_sphere, bestC);
end

% Run the SOUND algorithm to the average data to estimate the
% channel-specific noise levels:
chans = setdiff(1:size(data_cat,1),bestC);
[~,x,sigmas,~] = SOUND(datatmp(chans,:), LFM_sphere_tmp(chans,:), iter, lambdaValue);
 
% Re-reference the data and the lead field to the channel average:
ERP_cleaned = LFM_sphere_mean*x;

% REMOVED AS UNNECESSARY AND AS A POTENTIAL FAULT WITH MULTIPLE
% 
%if size(EEG.data,3) > 1 % Data are multiple single trials
    
    % Apply SOUND to the single trial data using the estimated noise levels
    fprintf('Applying SOUND to single trials\n');
    
    if ~isempty(replaceChans)
        LFM_sphere_goodC = LFM_sphere(EEG.goodC,:);
        [LFM_sphere_tmp] = ref_best(LFM_sphere_goodC, bestC);
    else
        [LFM_sphere_tmp] = ref_best(LFM_sphere, bestC);
    end
    for k = 1:size(EEG.data,3)
        [datatmp] = ref_best(EEG.data(:,:,k), bestC);
        [~, x] = correct_with_known_noise(datatmp(chans,:),...
            LFM_sphere_tmp(chans,:), lambdaValue, sigmas);
        EEGtmp(:,:,k) = LFM_sphere_mean*x;
    end
    EEG.data = EEGtmp;


% Store the leadfield matrix
EEG.soundLeadfield = LFM_sphere_mean;

% Update EEGLAB fields if replacing electrodes
if ~isempty(replaceChans)
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs = EEGallChans.chanlocs;   
end

% Update the reference in the bookkeeping:

EEG = pop_reref( EEG, []);

 % If more than one datasets are cleaned simultaneously:
if ~isempty(multipleConds)

    multipleCondsClean = cell(1,length(multipleConds));
    
    for eegx = 1:length(multipleConds)

        tempCond = evalin('base', multipleConds{eegx});
        
        % Check if also missing channels will be interpolated
        if ~isempty(replaceChans)
            EEGtmp = zeros(size(LFM_sphere_mean,1) , size(tempCond.data, 2));
        else
            EEGtmp = tempCond.data;
        end
        
        % Correct each of the extra datasets to be cleaned with EEG. 
        for k = 1:size(tempCond.data,3)
            [datatmp] = ref_best(tempCond.data(:,:,k), bestC);
            [corrected_data, x] = correct_with_known_noise(datatmp(chans,:),...
                LFM_sphere_tmp(chans,:), lambdaValue, sigmas);
            EEGtmp(:,:,k) = LFM_sphere_mean*x;            
        end
        
        % Bookkeeping for the extra datasets:
        tempCond.data = EEGtmp;
        tempCond.soundLeadfield = LFM_sphere_mean;
        
        % Update EEGLAB fields if replacing electrodes
        if ~isempty(replaceChans)
            tempCond.nbchan = size(EEG.data,1);
            tempCond.chanlocs = EEGallChans.chanlocs;
        end
        
        % Store the names of all the other condition cleaned with EEG:
        EEG.multipleCondsClean = multipleCondsClean;
        
        % Save the multiple conditions directly to the workspace:
        assignin('base', [multipleConds{eegx},'_after_SOUND'] , tempCond);
        
        fprintf(['Saving the cleaned version of ', multipleConds{eegx}, ...
            ' into a new variable. \n'])

    end
    
end

fprintf('SOUND complete\n');

end

% ##### OTHER FUNCTIONS USED IN TESA_SOUND #######

function [LFM_sphere] = construct_spherical_lead_field(EEG,r0,r1,r2,r3,sig1,sig2,sig3,DipN)
% This functions creates a spherical lead-field that can be used in SOUND
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
   sig2 = 0.33/100;
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
        if isempty(EEG.chanlocs(i).X) || isempty(EEG.chanlocs(i).Y) || isempty(EEG.chanlocs(i).Z)
            error('All the coordinates for all the included EEG channels are not specified! Unable to build a spherical lead field.')
        else
            elec_coords(k,1) = EEG.chanlocs(i).X;
            elec_coords(k,2) = EEG.chanlocs(i).Y;
            elec_coords(k,3) = EEG.chanlocs(i).Z;
            k = k+1;
        end
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

rng('default'); % Reset the random number generator to fix dipoles. This minimises variability in SOUND output.
P=randn(3,DipN);
rng('shuffle'); % Re-shuffle the random number generator
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

function [y_solved, sigmas] = DDWiener(data)  
% This function computes the data-driven Wiener estimate (DDWiener),
% providing the estimated signals and the estimated noise-amplitudes
%
% .........................................................................
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

function [data_bestC] = ref_best(data, bestC)
% This function re-references the data to the channel with the least noise, 
% indicated with bestC.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

data_bestC = data - repmat(data(bestC,:),[size(data,1),1]);

end

function [data_ave] = ref_ave(data)
% This function re-references the data to the average reference.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

data_ave = data - repmat(mean(data,1),[size(data,1),1]);

end

function [corrected_data, x, sigmas,dn] = SOUND(data, LFM, iter,lambda0)
% This function performs the SOUND algorithm for a given data.
%
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................

chanN = size(data,1);
sigmas = ones(chanN,1);
[y_solved, sigmas] = DDWiener(data);

if nargin<4
    lambda0 = 1;
end

% Number of time points
T = size(data,2);

% Going through all the channels as many times as requested
for k=1:iter
    sigmas_old = sigmas;
    
    disp(['Performing SOUND. Iteration round: ',num2str(k)])
    
    %Evaluating each channel in a random order
    for i=randperm(chanN)
        chan = setdiff(1:chanN,i);
        % Defining the whitening operator with the latest noise
        % estimates
        W = diag(1./sigmas);

        % Computing the whitened version of the lead field
        WL = (W(chan,chan))*(LFM(chan,:));
        WLLW = WL*WL';

        % Computing the MNE, the Wiener estimate in the
        % studied channel, as well as the corresponding noise estimate
        x = (WL)'*((WLLW + lambda0*trace(WLLW)/(chanN-1)*eye(chanN-1))\((W(chan,chan))*(data(chan,:))));
        y_solved = LFM*x;
        sigmas(i) = sqrt((y_solved(i,:)-data(i,:))*(y_solved(i,:)-data(i,:))')/sqrt(T);
    end
    
    % Following and storing the convergence of the algorithm
    dn(k) = max(abs(sigmas_old - sigmas)./sigmas_old);

end

% Final data correction based on the final noise-covariance estimate.

W = diag(1./sigmas);
WL = W*LFM;
WLLW = WL*WL';
x = WL'*((WLLW + lambda0*trace(WLLW)/chanN*eye(chanN))\(W*data));

corrected_data = LFM*x;

end

function [corrected_data, x] = correct_with_known_noise(data, LFM, lambda0,  sigmas)

% This function corrects a data segment when the noise distribution is already known.
% .........................................................................
% 24 September 2017: Tuomas Mutanen, NBE, Aalto university  
% .........................................................................
    chanN = length(sigmas);
    W = diag(1./sigmas);
    WL = W*LFM;
    WLLW = WL*WL';
    x = WL'*((WLLW + lambda0*trace(WLLW)/chanN*eye(chanN))\(W*data));
    corrected_data = LFM*x;

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

% When using SOUND only to clean the channel signals, and not to
% improve source estimation per se, we can replace the full lead-field matrix with a square matrix as follows.
% This enhances the computation significantly. 
[U,S,~] = svd(V,'econ');
V = U*S;

end

function [L,dL]=legen(N,x);
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

function gam=potgam(r,sig,nmax);
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
