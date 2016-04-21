% tesa_sortcomps()      - sorts ICA components by percentage variance explained by
%                           each time course. Use this to help identify
%                           muscle/decay artefacts using pop_selectcomps.
% 
% 
%                           Note that this was designed primarily for fastica which
%                           does not sort components. Infomax ICA
%                           algorithms already sort components so this will be of limited use. 
% Usage:
%   >>  EEG = tesa_sortcomps( EEG );
%   >>  [EEG, varsPerc] = tesa_sortcomps( EEG ); % Outputs the percentage of overall variance accounted for by the average time course of each component
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
% 
% Outputs:
%   EEG                 - EEGLAB EEG structure
%   varsPerc            - vector with % variance accounted for by the
%                       time course (averaged across trials) for each
%                       component. Output is sorted by variance.
%
% See also:
%   tesa_autocompselect 

% Copyright (C) 2016  Nigel Rogasch, Monash University,
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

function [EEG, varsPerc] = tesa_sortcomps( EEG )

%Check input for ICA weights
if isempty(EEG.icawinv)
	error('No ICA weights. Please run ICA first.');
end

vars = arrayfun(@(x)var(mean(eeg_getdatact(EEG, 'component', [x], 'projchan', []),3)),1:size(EEG.icawinv,2)); %extracts time course for each component
vars_norm = vars/sum(vars)*100; %calculates the % variance of each component relative to all components
[xSorted, ixsSort] = sort(vars_norm, 'descend'); %ranks components based on %variance
varsPerc = vars_norm(ixsSort); 

EEG.icawinv = EEG.icawinv(:,ixsSort); %alters inverse weights matrix based on component variance order
EEG.icaweights = EEG.icaweights(ixsSort,:); %alters weights matrix based on component variance order
if ~isempty(EEG.icaact)
    EEG.icaact = EEG.icaact(ixsSort,:,:); %alters time course matrix based on component variance order
end

fprintf('ICA weights sorted by time course variance\n');
        
end
