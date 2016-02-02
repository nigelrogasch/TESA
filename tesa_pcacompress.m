% tesa_pcacompress()    - Compresses data to n-dimensions as advocated in
%                           the following papers:
%
%                           Korhonen, hernandez-pavon et al (2011) Removal of large muscle artifacts 
%                           from transcranial magnetic stimulation-evoked EEG 
%                           by independent component analysis. Med Biol Eng
%                           Compt, 49:397-407.
%
%                           Hernandez-Pavon et al (2012) Uncovering neural  
%                           independent components from highly artifactual 
%                           TMS-evoked EEG data. J Neurosci Meth,
%                           209:144-57
% 
% Usage:
%   >>  EEG = tesa_pcacompress( EEG );
%   >>  EEG = tesa_pcacompress( EEG , compVal);
%
% Inputs:
%   EEG                 - EEGLAB EEG structure
%   compVal             - [Integer] Number of dimensions to compress data
%                           Default = 25.
% 
% Outputs:
%   EEG                 - EEGLAB EEG structure
%
% See also:
%   SAMPLE, EEGLAB 
%
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG,inMat,outMat] = tesa_pcacompress( EEG, compVal )

%Sets truncate to default value if not specified
if nargin <2;
	compVal = 25;
end
if isempty(compVal)
    compVal = 25;
end

%Reshapes 3D matrix to 2D
inMat=reshape(EEG.data,size(EEG.data,1),[],1);

%Checks that number of dimensions is larger than compression value
rankMat = rank(inMat);
if rankMat <= compVal
    error('Dimension of data (%d) is lower than the compression value (%d). Function terminated.', rankMat, compVal);
end

%Runs singular value decomposition
[U,S,V]=svd(inMat*inMat');
d=diag(S);

%Compresses data by truncating to 'compVal' dimensions 
C=U(:,1:compVal)*U(:,1:compVal)';
outMat=C*inMat;

%Reshapes 2D matrix to 3D matrix
EEG.data = reshape(outMat,size(EEG.data,1),size(EEG.data,2),size(EEG.data,3));

%Figures before and after suppression
h1=figure; subplot(1,2,1), bar(d);grid;set(gca,'fontsize',16);title('Data before suppression'); 
xlabel ('Dimensions'); ylabel ('Amplitude (A.U)');
subplot(1,2,2), bar(d(1:compVal));grid;set(gca,'fontsize',16);title('Data after suppression'); 
xlabel ('Dimensions'); ylabel ('Amplitude (A.U)');

fprintf('Data compressed to %d dimensions\n', compVal);
        
end
