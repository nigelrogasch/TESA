% pop_tesa_edm()    - find artefactual components automatically by using the 
%                     EDM algorithm and FastICA as advocated in 
%                     the following paper:
%         
%                           Korhonen, Hernandez-Pavon et al (2011) Removal of 
%                           large muscle artifacts from transcranial magnetic
%                           stimulation-evoked EEG by independent component
%                           analysis. Med Biol Eng Compt, 49:397-407.
%                           
% 
% You must install the FastICA toolbox and include it in the Matlab path
% before using this function
% For downloading the FastICA toolbox go to:
% http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
%
% Usage:
%   >>  EEG = tesa_edm(EEG, chanlocs);
%   >>  EEG = tesa_edm(EEG, chanlocs, compVal, Nic);
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   chanlocs            - Channel locations 
%   compVal             - [Integer] Number of dimensions to compress data
%                           Default = 25.
%   Nic              - [Integer] Number of independent components to look for
%                           Default = rank(EEG)-5 (to make sure the
%                           algorithm converges).
% 
% Outputs:
%   EEG                 - EEGLAB EEG structure, data after removing
%                         artefactual ICs
%
% See also:
%   SAMPLE, EEGLAB 
%
% Copyright (C) 2016  Julio Cesar Hernandez Pavon, Aalto University,
% Finland, julio.hpavon@gmail.com
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
 function [EEG com] = pop_tesa_edm(EEG, chanlocs, compVal, Nic);
 
 
com = '';          

%check that data is present
if isempty(EEG.data);
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 3
   
   geometry = { 1 [1 0.5] [1 0.5] 1};
    
    uilist =  {{'style', 'text', 'string', 'The number of ICs must be smaller or equal than the remaining dimension of the data','fontweight','bold'}...
              {'style', 'text', 'string', 'Number of dimensions for compression'} ...
              {'style', 'edit', 'string', '25'}...
              {'style', 'text', 'string', 'Number of ICs to search'} ...
              {'style', 'edit', 'string', '25'}...
              {'style', 'text', 'string', 'Press right click for selecting neuronal components and any other key for selecting artefactual ones','fontweight','bold'}};  
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Compress data dimensionality-- pop_tesa_edm()', 'helpcom', 'pophelp(''tesa_edm'')');
      
    %extract compression and IC values
    compVal = str2num(result{1,1});
    Nic=str2num(result{1,2});
    
    if Nic  > compVal
    error('Please enter a number of ICs smaller that the number of dimensions');
    end

end


%Compressed data

EEG = tesa_edm(EEG, chanlocs, compVal, Nic);

com = sprintf('%s = pop_tesa_edm( %s, %s, %s );', inputname(1), inputname(1), compVal,Nic );

end
