% pop_tesa_peakoutput() - returns the results for the peak analysis in a table
%                       in the workspace. Users can also opt to have the
%                       average amplitude incorporating data points either 
%                       side of the peak instead of the absolute
%                       peak amplitude.
%
% Usage:
%   >>  output = pop_tesa_peakoutput( EEG ) %use pop-up window
%   >>  output = pop_tesa_peakoutput( EEG, varargin )
%
% Inputs:
%   EEG             - EEGLAB EEG structure
% 
% Optional input pairs
%   'averageWin',int  - integer describing a time window +/- the peak (in
%                       ms) in which an average amplitude will be taken. If left
%                       empty, the absolute amplitude at the peak will be
%                       returned.
%     
% Outputs:
%   output             - table with results from peak analysis
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

 function [output, com] = pop_tesa_peakoutput( EEG, varargin )

com = '';          

%check that data is present
if isempty(EEG.data)
    error('Data is empty');
end

% pop up window
% -------------
if nargin < 2
    
    geometry = {1 1 1 [1 0.5] 1 1};

    uilist = {{'style', 'text', 'string', 'Output peak analysis results','fontweight','bold'} ...
              {'style', 'text', 'string', 'Results from peak analyses will appear in the output variable in the workspace.'} ...
              {}...
              {'style', 'text', 'string', 'Window for averaging over peak (+/- ms) [optional]','fontweight','bold'} ...
              {'style', 'edit', 'string', ''}...
              {'style', 'text', 'string', '   Instead of returning the absolute amplitude, entering a value here will return','fontangle','italic'}...
              {'style', 'text', 'string', '   the average value of the window around the peak','fontangle','italic'}};
             
    result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'Output peak -- pop_tesa_peakoutput()', 'helpcom', 'pophelp(''tesa_peakoutput'')');
      
    %Run script
    if ~strcmp(result{1,1},'')
        averageWin = str2num(result{1,1});
    end       
   
end

%Run script from input
if nargin > 1;
    output = tesa_peakoutput(EEG,varargin{:});
    com = sprintf('%s = pop_tesa_peakoutput( %s, %s, %s );', inputname(1), inputname(1), vararg2str(varargin) );
end

%perform roi analysis and return the string command
if nargin == 1;
    if strcmp(result{1,1},'')
        output = tesa_peakoutput(EEG);
        com = sprintf('%s = pop_tesa_peakoutput( %s, ''averageWin'', ''0'' );', inputname(1), inputname(1));
    else
        output = tesa_peakoutput(EEG,'averageWin',averageWin);
        com = sprintf('%s = pop_tesa_peakoutput( %s, ''averageWin'', %s );', inputname(1), inputname(1), num2str(averageWin) );
    end
end

end
