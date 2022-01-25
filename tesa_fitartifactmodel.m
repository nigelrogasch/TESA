% tesa_fitartifactmodel() - this script fits a model to the capacitive
%                           discharge artifact following TMS and then 
%                           subtracts it from the data as described in:

% Freche D, Naim-Feil J, Peled A, Levit- Binnun N, Moses E (2018) A 
% quantitative physical model of the TMS-induced discharge artifacts in
% EEG. PLoS Comput Biol 14(7): e1006177. 
% https://doi.org/10.1371/journal.pcbi.1006177
%
% Please cite this article if you use this function.
%
% Usage:
% To call this script, use pop_tesa_fitartifactmodel.m or see in the script
% below for further instructions.

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
%
function [ V, u, PrepulseOfs ] = tesa_fitartifactmodel(EEGdata, tPulse, tPulseDuration, tSkip, tSelect, t_fitmodel_p, t_fitmodel_n, t_lincomb, srate)

% 
% Input:
%   EEGdata(iE, t)      EEG data matrix, where iE is the electrode index and t is the time (in samples) index
% 
%   tPulse              sample point index marking TMS pulse application in EEGdata
%   tPulseDuration      duration of the pulse artifact
%   tSkip               number of samples to be skipped after tPulseDuration samples following tPulse
%   tSelect             sample point used to select the two traces (indexed by I) following tPulse (see below)
% 
%   t_fitmodel_p        cell array of time ranges for model fitting.
%                       The model is fit to the positive traces from index set I (see below).
%   t_fitmodel_n        cell array of time ranges for model fitting.
%                       The model is fit to the negative traces from index set I (see below).
%   t_fitmodel_p{i}     time range for trace i in I
%   t_fitmodel_n{i}     time range for trace i in I
%   t_lincomb           time range where the linear coefficients for V(t,i) are to be found (see below)
% 
% Output:
%   V(t,i)              a matrix-valued function such that V(t,i) (for t from 0 to infinity)
%                       are the functions such that all artifacts can be obtained as linear combination of them
%   u(i,iE)             the coefficients u(i,iE) for the linear combination such that
%                           sum_(i in I) u(i,iE)*V(t,i)
%                       approximates the artifact (as function of t from 0 to infinity) for EEGdata(iE, :)
%                       where iE is the electrode index
%   PrepulseOfs(iE)     the EEG trace offset right before the TMS pulse


% For every EEG data trace, we subtract the pre-pulse offset
% (obtained as the mean of a few samples immediately preceding the pulse)
% to take care of drifts in the data.
PrepulseOfs = mean(EEGdata(:,tPulse - (1:round(srate/1000)) ),2);

% Select the traces for fitting:
% Suitable candidates are artifacts which show sortof extreme and opposite behavior.
% For example, we select here the N_p most positive and the N_n most negative artifact
% evaluated at the the point tSelect after the pulse application.
% tSelect should be at least as long as the pulse duration and possibly an additional
% few samples more to make sure the artifact shape is not affected by any fast transients
% such as hardware specifics such as internal filters, etc.
[~, idataLst] = sort(EEGdata(:,tPulse + tSelect)' - PrepulseOfs', 'descend');
I = idataLst([(1:length(t_fitmodel_p)) length(idataLst)+1-(1:length(t_fitmodel_n))]);
t_fitmodel = {t_fitmodel_p{:} t_fitmodel_n{:}};

% Parameters for rational function definition:
% The artifact model requires that dd-dn = 2 (namely dn = 1 and dd = 3).
% Asymptotically, dn = 0 and dd = 2 is also possible, which in this case may require
% increasing tSkip to exclude the O(t^-3) transients which then cannot be fitted anymore.
dn = 1; % numerator degree
dd = 3; % denominator degree
ratcoeff = zeros(dn+1+dd+1-1,length(I));
for iE = I
    
    t = t_fitmodel{I==iE} + tSkip;
    
    y = EEGdata(iE, tPulse + tPulseDuration + t) - PrepulseOfs(iE);
    
    % Fitting:
    % We fit a rational function V with numerator degree dn and denominator degree de to the
    % data y. We use an iterative nonlinear fitting algorithm to find the coefficients of V.
    % This requires a suitable starting point P0, which is obtained as follows.
    % Approximate the nonlinear fitting problem by a linear problem of the form:
    %       find P0 which minimizes norm(C*P0-d,2) such that A*P0 <= b
    % which can be solved by quadratic optimization.
    
    % Solve linearized constrained minimization problem:
    C = zeros(length(t),dn+1 + dd+1-1);
    for k = 0:dn
        C(:,k+1) = (t.^k);
    end
    for k = 0:dd-1
        C(:,k+1+dn+1) = -y.*( t.^k );
    end
    d = ( y.*(t.^dd) )';
    A = zeros(dd+dn,dn+1 + dd+1 -1);
    for k = 0:dn-1
        A(k+1,k+1) = -1;
    end
    for k = 0:dd-1
        A(k+1+dn,k+1+dn+1) = -1;
    end
    b = zeros(1,dd+dn);
    loptions = optimset('algorithm','interior-point','display','on');
    P0 = lsqlin(C,d,A,b, [], [], [], [], [], loptions);
    
    % Solve original nonlinear constrained minimization problem:
    fitfunction = @(P, time) arrayfun(@(t) ( t.^(0:dn)*P(1:dn+1) )./( t.^(0:dd)*[P(dn+1+(1:dd)) ; 1] ), time);
    lb = zeros(dn+1+dd,1);
    lb(1:dn+1) = -Inf(dn+1,1);
    nloptions = optimset('algorithm','trust-region-reflective','display','on');
    [Pmin, ~, residual] = lsqcurvefit(fitfunction, P0, t', y', lb, [], nloptions);
    
    ratcoeff(:,I==iE) = Pmin;
    
    % *** uncomment to plot ***
    %     hold on
    %     zoom on
    %     
    %     % rational function corresponding to P0
    %     fP0 = @(time) arrayfun(@(t) (t.^(0:dn)*P0(1+(0:dn)) )./( t.^(0:dd)*[P0(1+dn+1+(0:dd-1)) ; 1] ), time);
    %     % rational function corresponding to Pmin
    %     fP = @(time) fitfunction(Pmin, time);
    %     
    %     plot(t, y,       'color', [1 0 0  0.3], 'linewidth', 3) % data used for fitting
    %     plot(t, fP0(t),  'color', [0 1 0  0.4], 'linewidth', 2) % starting point
    %     plot(t, fP(t),   'color', [0 0 1  0.8], 'linewidth', 1) % fit (reconstructed artifact)
    %     plot(t, residual,'color', [0 0 0  0.5]) % residuals (artifact-free EEG data plus fitting error)
    %     
    %     drawnow;
end

% We conveniently express the computed artifacts for the each trace in I
% as a matrix-valued function V(t,i) with
% trace number i and time t (in samples) ranging from 0 to infinity
fP_cell = @(time,index) arrayfun(@(i) arrayfun(@(t) (t.^(0:dn)*ratcoeff(1:dn+1,i) )./( t.^(0:dd) * [ratcoeff(dn+1+(1:dd),i);1] ), time), index, 'uniformoutput', false);
V = @(time,index) cell2mat( fP_cell(tSkip + time,index)' )';

% The artifact in each EEG trace can now be expressed as linear combination
% of V(t,i) with coefficients u(i,iE) for trace iE.
% Because the EEG traces include also other contributions (such as physiological
% data of interest like EEG, muscle artifacts, etc) and noise, u can only be approximated.
% We therefore compute u in the sense of least squares for a prescribed time span t_lincomb.
A = V(t_lincomb,1:length(I));
u = zeros(length(I),size(EEGdata,1));
for iE = 1:size(EEGdata,1)
    u(:,iE) = A \ (EEGdata(iE,tPulse + tPulseDuration + tSkip + t_lincomb)-PrepulseOfs(iE))';
end

end
