% pop_tesa_sound() - performs noise suppression on the data using the SOUND algorithm.
%                   Note that data are rereferenced to average during SOUND
%                   For further details on SOUND, or if you use this
%                   algorithm, please cite:
%                   
%                   Mutanen, TP, et al. Automatic and robust noise suppression 
%                   in EEG and MEG: The SOUND algorithm. Neuroimage. 2018 
%                   Feb 1;166:135-151.   
%
% Usage:
%   >>  EEG = pop_tesa_sound( EEG ); % pop up window
%   >>  EEG = pop_tesa_sound( EEG ); % run SOUND using default values and a spherical leadfield
%   >>  EEG = pop_tesa_sound( EEG, 'key1',value1... );% run SOUND using customised inputs
%
% Inputs:
%   EEG             - (required) EEGLAB EEG structure. EEG.data can include single trials, or the average of across trials. 
%
% Optional input pairs (varargin):
%   'lambdaValue', double       - (optional) double providing the lambda value for SOUND. Default is 0.1.
%                               Example: 0.1
%   'iter', int                 - (optional) integer providing the number of iterations for SOUND. Default is 5.
%                               Example: 5
%   'leadfieldInVar', matrix    - (optional) n x m matrix with individualised leadfield where n = channels and m = dipoles. 
%                               IMPORTANT - channel number and order must match the channels in the analyzed EEGLAB structure OR 
%                               The channel order should be sorted to match EEGLAB using the 'leadfieldChansFile' input below.
%                               Note: The data will be re-referenced to an average output.
%                               Example: leadfield [e.g. where leadfield = 62 x 15006 matrix (channel x dipole) stored in MATLAB workspace]
%   'leadfieldInFile', str      - (optional) file path and name of individualised leadfield matrix where n = channels and m = dipoles.
%                               This command can be used as an alternative to 'leadfieldInVar' which will load the required leadfield matrix in to the workspace.
%                               Note that the file must only contain the leadfield matrix and nothing else.
%                               IMPORTANT - channel number and order must match channels in EEG structure OR the channel
%                               order should be sorted to match EEGLAB using the 'leadfieldChansFile' input below.
%                               Example: 'C:\data\myLeadfield.mat' which contains leadfield variable [e.g. where leadfield = 62 x 15006 matrix (channel x dipole)]
%   'leadfieldChansFile', str   - (optional) file path and name of cell array listing the channel order of the individual leadfield matrix. 
%                               If called, this command will load the indicated cell array and then sort the leadfield matrix channel order to match the EEGLAB data channel order.
%                               This is important if the leadfield has been generated in another program which stores the EEG data in a different order.
%                               Note that the channel names must match those in the EEGLAB data.
%                               Example: 'C:\data\myLeadfieldChans.mat' which contains cell array with channel order.
%   'multipleConds', cell array - (optional) cell string array indicating the variables of the additional datasets that need to be cleaned simultaneously with input variable EEG.  
%                               Use this input if you would like to apply SOUND identically to multiple conditions (e.g. pre and post an intervention). 
%                               Conditions must first be epoched. Epochs from different conditions identified in the 'multipleConds' input will be separated, 
%                               averaged across epochs within a condition, and then concatenated in time before performing SOUND.
%                               Example: {'EEG2', 'EEG3'}; 
%                               For instance, EEG1 = pop_tesa_sound( EEG1, 'multipleConds', {'EEG2', 'EEG3'} ); would clean simultaneously identically the
%                               datasets EEG1, EEG2, EEG3. Note, only the cleaned version of EEG1 is returned as an output and the datasets
%                               defined within the 'multipleConds' variable are directly saved to the workspace with new names EEG2_after_SOUND and EEG3_after_SOUND.
%                               The names of the additional datasets will be also saved into an extra field in the output variable EEG, i.e., EEG.multipleCondsClean                             
%   'replaceChans',             - (optional) path and file name for EEGLAB EEG structure containing all required output channels stored in the EEG.chanlocs field.
%                               This command uses SOUND to replace missing channels (e.g. ones that have been removed earlier in the processing pipeline).
%                               If called, this command will replace any channels missing from the current data set relative to the data set called by 'replaceChans'.
%                               Example: 'C:\data\fullChannelEegData.set';
%                               NOTE: This works also with 'multipleConds' option with the restriction that the same channels must have been removed in all
%                               datasets prior to SOUND.
%                               NOTE: When using your own leadfield, it should have channels matching the 'replaceChans'-input dataset. 

% 
% Outputs:
%   EEG             - EEGLAB EEG structure
%
% Examples
%   EEG = pop_tesa_sound( EEG, 'lambdaValue', 0.2, 'iter', 10 ); % run SOUND using customised input values 
%   EEG = pop_tesa_sound( EEG, 'leadfieldInVar', leadfield, 'leadfieldChansFile', 'C:\data\myLeadfieldChans.mat' ); % run SOUND using default values and an individual leadfield matrix (stored in variable called leadfield). Sort the leadfield channel order (stored in 'C:\data\myLeadfieldChans.mat') to match the EEGLAB data channel order.
%   EEG = pop_tesa_sound( EEG, 'leadfieldInFile', 'C:\data\myLeadfield.mat' ); % run SOUND using default values and an individual leadfield matrix (stored in variable called leadfield loaded from 'C:\data\myLeadfield.mat')
%   EEG1 = pop_tesa_sound( EEG1, 'multipleConds', {'EEG2', 'EEG3'} ); % Clean simultaneously identically the  datasets EEG1, EEG2, EEG3. 
%   EEG = pop_tesa_sound( EEG, 'replaceChans', 'C:\data\fullChannelEegData.set' ); % run SOUND using default values and replacing missing channels in the current data set (all channels defined in the loaded file). 
% 
% See also:
%   tesa_sound, tesa_sspsir, pop_tesa_sspsir 

% Copyright (C) 2018 Tuomas Mutanen, University of Glasgow, Tuomas.Mutanen@glasgow.ac.uk;   
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

 function [EEG, com] = pop_tesa_sound( EEG, varargin )

com = '';          

%check that EEG has been provided as input
if nargin < 1
    error('Please provide EEGLAB structure EEG as input');
end

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    geometry = {1 [1 0.5] [1 0.5] 1 1 1 [1 0.2] 1 1 1 [1 0.2] 1 1 1 [1 0.2] 1 1 [1 0.5] 1};
    
    commandload1 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
    'if filename ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''leadfield''), ''string'', [ filepath filename ]);' ...
    'end;' ...
    'clear filename filepath tagtest;' ];

    commandload2 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
    'if filename ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''order''), ''string'', [ filepath filename ]);' ...
    'end;' ...
    'clear filename filepath tagtest;' ];

    commandload3 = [ '[filename, filepath] = uigetfile(''*'', ''Select a text file'');' ...
    'if filename ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''elec''), ''string'', [ filepath filename ]);' ...
    'end;' ...
    'clear filename filepath tagtest;' ];

    uilist = {{'style', 'text', 'string', 'Perform SOUND noise suppression','fontweight','bold'} ... %1
              {'style', 'text', 'string', 'Lambda value'} ... %[1
              {'style', 'edit', 'string', '0.1', 'tag', 'chans' } ... %0.5]
              {'style', 'text', 'string', 'Iterations',}... %[1
              {'style', 'edit', 'string', '5', 'tag', 'chans' } ... %0.5]
              {}... %1
              {'style', 'text', 'string', 'Optional inputs','fontweight','bold'}...%1
              {'style', 'text', 'string', 'Load individual lead field matrix  ','fontweight','bold'} ... %1
              { 'style' 'edit'       'string' '' 'tag' 'leadfield' } ... %[1
              { 'style' 'pushbutton' 'string' '...' 'callback' commandload1 }... %0.2]
              {'style', 'text', 'string', 'Select .mat file with channel x dipole matrix','fontAngle','italic'}...%1
              {}... %1
              {'style', 'text', 'string', 'Load channel order for lead field matrix','fontweight','bold'} ... %1
              { 'style' 'edit'       'string' '' 'tag' 'order' } ... %[1
              { 'style' 'pushbutton' 'string' '...' 'callback' commandload2 }... %0.2]
              {'style', 'text', 'string', 'Select .mat file with channel order in cell array','fontAngle','italic'}...%1
              {}... %1
              {'style', 'text', 'string', 'Replace missing electrodes','fontweight','bold'} ... %1
              { 'style' 'edit'       'string' '' 'tag' 'elec' } ... %[1
              { 'style' 'pushbutton' 'string' '...' 'callback' commandload3 }... %0.2]
              {'style', 'text', 'string', 'Select .set file containing all required channels','fontAngle','italic'}...%1
              {}...%1
              {'style', 'text', 'string', 'Apply to multiple conditions','fontweight','bold'} ... %[1
              {'style', 'edit', 'string', '', 'tag', 'chans' }... %0.5]
              {'style', 'text', 'string', 'Enter the list of names of all datasets to be cleaned (e.g. EEG_2 EEG_3)','fontAngle','italic'}}; %1
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Perform SOUND -- pop_tesa_sound()', 'helpcom', 'pophelp(''pop_tesa_sound'')');
    if isempty(result), return; end;
    
    %Extract data for single pulse artifact find
    if strcmp(result{1,1},'')
        lambda = [];
        lambdaString = '[]';
    else
        lambda = str2num(result{1,1});
        lambdaString = result{1,1};
    end
    
    if strcmp(result{1,2},'')
        iter = [];
        iterString = '[]';
    else
        iter = str2num(result{1,2});
        iterString = result{1,2};
    end
    
    if strcmp(result{1,3},'')
        leadfieldInFile = [];
        leadfieldInFileString = '[]';
    else
        leadfieldInFile = result{1,3};
        leadfieldInFileString = ['''',result{1,3},''''];
    end
    
    if strcmp(result{1,4},'')
        leadfieldChansFile = [];
        leadfieldChansFileString = '[]';
    else
        leadfieldChansFile = result{1,4};
        leadfieldChansFileString = ['''',result{1,4},''''];
    end
    
    if strcmp(result{1,5},'')
        replaceChans = [];
        replaceChansString = '[]';
    else
        replaceChans = result{1,5};
        replaceChansString = ['''',result{1,5},''''];
    end
    
    if strcmp(result{1,6},'')
        multipleConds = [];
        multipleCondsString = '[]';
    else
        multipleConds = strsplit(result{1,6});
        multipleCondsString = ['''',result{1,6},''''];
        multipleCondsString = strrep(multipleCondsString,' ',''',''');
        multipleCondsString = strrep(multipleCondsString,' ','');
    end
  
end

%Run script from input
if nargin > 1 
    EEG = tesa_sound(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_sound(%s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

%Run script from GUI input
if nargin < 2
    EEG = tesa_sound(EEG,'lambdaValue',lambda,'iter',iter,'leadfieldInFile',leadfieldInFile,'leadfieldChansFile',leadfieldChansFile,'replaceChans',replaceChans,'multipleConds',multipleConds);
    
    if isempty(multipleConds)
       com = sprintf('%s = pop_tesa_sound(%s, ''lambdaValue'', %s, ''iter'', %s, ''leadfieldInFile'' ,%s, ''leadfieldChansFile'', %s, ''replaceChans'', %s, ''multipleConds'', %s);', inputname(1), inputname(1), lambdaString, iterString, leadfieldInFileString, leadfieldChansFileString, replaceChansString, multipleCondsString);
    else
       com = sprintf('%s = pop_tesa_sound(%s, ''lambdaValue'', %s, ''iter'', %s, ''leadfieldInFile'', %s, ''leadfieldChansFile'', %s, ''replaceChans'', %s, ''multipleConds'', {%s});', inputname(1), inputname(1), lambdaString, iterString, leadfieldInFileString, leadfieldChansFileString, replaceChansString, multipleCondsString);
    end
end

 end
