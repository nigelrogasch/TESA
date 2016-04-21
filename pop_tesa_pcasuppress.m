% pop_tesa_pcasuppress()    - suppresses data by top PCA n-dimensions. Note that
%                      PCA is generated from average data, but subtracted
%                      from single trial data.
%                          
%                      Hernandez-Pavon et al (2012) Uncovering neural  
%                      independent components from highly artifactual 
%                      TMS-evoked EEG data. J Neurosci Meth,
%                      209:144-57
% 
% Usage:
%   >>  EEG = pop_tesa_pcasuppress( EEG ); % pop up window
%   >>  EEG = pop_tesa_pcasuppress( EEG, timeWin );
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   timeWin             - [Vector] - required. Time window for data suppression in ms.
%                       [start, end]
% Outputs:
%   EEG                 - EEGLAB EEG structure
%
% Example:  
%   EEG = pop_tesa_pcasuppress( EEG, [11, 50] );
% 
% See also:
%   pop_tesa_pcacompress, pop_tesa_edm 

% Copyright (C) 2016  Nigel Rogasch & Julio C. Hernandez-Pavon
% Monash University and Aalto University
% nigel.rogasch@monash.edu; julio.hpavon@gmail.com
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USAA

 function [EEG com] = pop_tesa_pcasuppress( EEG, timeWin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
   
    geometry = {1 [1 0.3]};
    
    uilist = {{'style', 'text', 'string', 'Suppress TMS-evoked artefacts','fontweight','bold'} ...
              {'style', 'text', 'string', 'Time range for suppression (ms): start,end'} ...
              {'style', 'edit', 'string', '11, 50'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Supress data -- pop_tesa_pcasuppress()', 'helpcom', 'pophelp(''pop_tesa_pcasuppress'')');
    if isempty(result), return; end;
    
    %extract time window
    timeWin = str2num(result{1,1});
    
%     if strcmp(result(1,1),'') || strcmp(result(1,2),'')
%         error('Please enter the time range for suppression');  
%     else
%         timeWin = cellfun(@str2num,(result(:,1:2)));
%     end 
%     
end

%supress data
EEG = tesa_pcasuppress( EEG, timeWin );

% return the string command
com = sprintf('%s = pop_tesa_pcasuppress( %s, %s);', inputname(1), inputname(1), mat2str(timeWin));

end
