% pop_tesa_pcasupress()    - Supresses data by top n-dimensions
%                          
%                           Hernandez-Pavon et al (2012) Uncovering neural  
%                           independent components from highly artifactual 
%                           TMS-evoked EEG data. J Neurosci Meth,
%                           209:144-57
% 
% Usage:
%   >>  EEG = pop_tesa_pcasupress( EEG );
%   >>  EEG = pop_tesa_pcasupress( EEG, timeWin );
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   timeWin             - [Vector] Time window for data supression in ms.
%                       [start, end]
% Outputs:
%   EEG                 - EEGLAB EEG structure
%
% See also:
%   SAMPLE, EEGLAB 

% Copyright (C) 2015  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu
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

 function [EEG com] = pop_tesa_pcasupress( EEG, timeWin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 0.5] [1 0.5]};
    
    uilist = {{'style', 'text', 'string', 'Time range for supression','fontweight','bold'} ...
              {'style', 'text', 'string', 'Start time (in ms)'} ...
              {'style', 'edit', 'string', ''} ...
              {'style', 'text', 'string', 'End time (in ms)'} ...
              {'style', 'edit', 'string', ''}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Supress data -- pop_tesa_pcasupress()', 'helpcom', 'pophelp(''tesa_pcasupress'')');
    
    %extract time window
   
    if strcmp(result(1,1),'') || strcmp(result(1,2),'')
        error('Please enter the time range for detrend');  
    else
        timeWin = cellfun(@str2num,(result(:,1:2)));
    end 
    
end

%rsupress data
EEG = tesa_pcasupress( EEG, timeWin );

% return the string command
com = sprintf('%s = pop_tesa_pcasupress( %s, %s );', inputname(1), inputname(1), mat2str(timeWin));

end
