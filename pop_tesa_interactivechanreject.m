% pop_tesa_interactivechanreject() - Interactively mark bad channels in EEG data.
%
% Usage:
%   >> [EEG, com] = pop_tesa_interactivechanreject(EEG)
%
% Inputs:
%   EEG         - EEGLAB EEG structure containing:
%                   EEG.data      : channels × time × trials (or channels × time)
%                   EEG.chanlocs  : channel locations (with .labels)
%
% Outputs:
%   EEG         - Updated EEG structure containing:
%                   EEG.badChans  : cell array of bad channel labels (cumulative)
%   com         - Command string for EEGLAB history
%
% Description:
%   pop_tesa_interactivechanreject() opens a GUI showing channel time-series.
%   Click to select (blue + brought to front), "Mark/Unmark as Bad" toggles
%   bad status (red). Two exit options:
%     - "Remove Bad Channels": removes them with pop_select and returns EEG.
%     - "Close Without Removing": returns EEG with EEG.badChans, but no removal.
%
% Author: Nigel Rogasch & ChatGPT (OpenAI), 2025
% -------------------------------------------------------------------------

function [EEG, com] = pop_tesa_interactivechanreject(EEG)

    com = ''; % default

    %% ---------------- Input checks ----------------
    if nargin < 1
        error('You must provide an EEG structure as input.');
    end
    if ~isstruct(EEG) || ~isfield(EEG, 'data') || ~isfield(EEG, 'chanlocs')
        error('Invalid EEG structure. Ensure EEG.data and EEG.chanlocs exist.');
    end

    % time axis: use EEG.times if present, otherwise fallback to 1:N
    if isfield(EEG, 'times') && ~isempty(EEG.times)
        time = EEG.times;
    else
        time = 1:size(EEG.data, 2);
    end

    % prepare plotting data (average across epochs if needed)
    data = EEG.data;
    if ndims(data) > 2
        data = mean(data, 3);
    end

    % channel labels
    if isfield(EEG.chanlocs, 'labels')
        chanLabels = {EEG.chanlocs.labels};
    else
        nChans = size(data, 1);
        chanLabels = arrayfun(@(x) sprintf('Ch%d', x), 1:nChans, 'UniformOutput', false);
    end
    nChans = numel(chanLabels);

    %% ---------------- Pre-existing bad channels ----------------
    if isfield(EEG, 'badChans') && ~isempty(EEG.badChans)
        allBadChans = EEG.badChans;
    else
        allBadChans = {};
    end

    % determine which of the stored bad channels are present now
    if ~isempty(allBadChans)
        [presentFlags, ~] = ismember(allBadChans, chanLabels);
        existingBadChans = allBadChans(presentFlags);
        missingBadChans  = allBadChans(~presentFlags);
        if ~isempty(missingBadChans)
            fprintf('\nThe following channels in EEG.badChans are not present in EEG.chanlocs (likely already removed):\n');
            fprintf('   %s\n', strjoin(missingBadChans, ', '));
            fprintf('These channels will remain in EEG.badChans for tracking.\n\n');
        end
    else
        existingBadChans = {};
        missingBadChans = {};
    end

    % session bad list (those we will display & edit in this GUI)
    sessionBadChans = existingBadChans;

    %% ---------------- Create figure ----------------
    f = figure('Name', 'TESA - Interactive Channel Rejection', ...
               'Color', 'w', 'Position', [200 150 1000 600], ...
               'CloseRequestFcn', @close_noaction); % handle close via buttons

    ax = axes('Parent', f, 'Position', [0.30 0.08 0.65 0.88]);
    hold(ax, 'on');

    % plot channels using EEG.times on x-axis
    h = gobjects(1, nChans);
    for i = 1:nChans
        h(i) = plot(ax, time, data(i,:), 'Color', [0.6 0.6 0.6], 'Tag', chanLabels{i});
    end

    % highlight pre-existing bad channels (sessionBadChans)
    if ~isempty(sessionBadChans)
        badIdxPresent = find(ismember(chanLabels, sessionBadChans));
        for ii = 1:numel(badIdxPresent)
            set(h(badIdxPresent(ii)), 'Color', 'r', 'LineWidth', 1.5);
        end
    else
        badIdxPresent = [];
    end

    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'Amplitude');
    title(ax, 'Click a line to select a channel');
    box(ax, 'on');
    hold(ax, 'off');

    %% ---------------- UI elements ----------------
    selectedLabel = uicontrol('Style', 'text', 'String', 'Selected: None', ...
                              'Position', [20 540 230 25], 'HorizontalAlignment', 'left');

    listBox = uicontrol('Style', 'listbox', 'Position', [20 140 230 380], ...
                        'String', sessionBadChans, 'Max', 2, 'Min', 0);

    uicontrol('Style', 'pushbutton', 'String', 'Mark/Unmark as Bad', ...
              'Position', [20 90 230 35], 'Callback', @toggleBad, 'FontSize', 10);

    uicontrol('Style', 'pushbutton', 'String', 'Remove Bad Channels', ...
              'Position', [20 50 230 35], 'BackgroundColor', [1 0.8 0.8], ...
              'Callback', @closeAndRemove, 'FontSize', 10);

    uicontrol('Style', 'pushbutton', 'String', 'Close Without Removing', ...
              'Position', [20 10 230 35], 'BackgroundColor', [0.8 1 0.8], ...
              'Callback', @closeOnly, 'FontSize', 10);

    % capture clicks on axes to select lines
    set(f, 'WindowButtonDownFcn', @selectLine);

    % state variables used by nested functions
    selectedIdx = [];                     % index of currently selected channel
    sessionBadIdx = find(ismember(chanLabels, sessionBadChans)); % indices of bad in current dataset
    userAction = '';                       % 'remove' or 'close' or '' if cancelled

    %% ---------------- Wait for user action ----------------
    uiwait(f);  % will resume when one of closeOnly / closeAndRemove / close_noaction calls uiresume

    %% ---------------- After GUI closed: prepare outputs ----------------
    % ensure EEG.badChans contains full cumulative list (existing + missing previously removed)
    EEG.badChans = unique([allBadChans, sessionBadChans], 'stable');

    % EEG.badIdx: indices of bad channels currently present in chanLabels
    [~, badIdx] = ismember(sessionBadChans, chanLabels);
    badIdx = badIdx(badIdx > 0);

    % handle user actions and set command string
    switch userAction
        case 'remove'
            if ~isempty(badIdx)
                % Remove channels by index via pop_select
                com = 'EEG = pop_tesa_interactivechanreject(EEG);';
                % actually run removal
                EEG = pop_select(EEG, 'nochannel', badIdx);
            else
                fprintf('No bad channels selected.\n');
                com = 'EEG = pop_tesa_interactivechanreject(EEG);';
            end
        case 'close'
            % do not remove; return EEG with EEG.badChans
            fprintf('Bad channels stored but not removed.\n');
            com = 'EEG = pop_tesa_interactivechanreject(EEG);';
        otherwise
            % user closed window using window manager or cancelled
            com = 'EEG = pop_tesa_interactivechanreject(EEG);';
            fprintf('TESA interactive rejection cancelled or closed without action.\n');
    end

    % ensure figure deleted
    if isvalid(f)
        delete(f);
    end

    %% ---------------- Nested callbacks ----------------

    function selectLine(~, ~)
        % Find nearest line to the click (in data units)
        cp = get(ax, 'CurrentPoint');
        xC = cp(1,1); yC = cp(1,2);

        minDist = inf; nearest = [];
        for idx = 1:nChans
            xd = get(h(idx), 'XData');
            yd = get(h(idx), 'YData');
            % find nearest x index
            [~, xi] = min(abs(xd - xC));
            d = abs(yd(xi) - yC);
            if d < minDist
                minDist = d;
                nearest = idx;
            end
        end

        if ~isempty(nearest)
            selectedIdx = nearest;
            selectedChan = chanLabels{selectedIdx};
            set(selectedLabel, 'String', ['Selected: ' selectedChan]);

            % Reset non-bad lines to grey and linewidth to 1
            for ii = 1:nChans
                set(h(ii), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
            end

            % re-apply red for session bad channels
            if ~isempty(sessionBadIdx)
                for ii = 1:numel(sessionBadIdx)
                    set(h(sessionBadIdx(ii)), 'Color', 'r', 'LineWidth', 1.5);
                end
            end

            % highlight selected (blue) and bring to front
            set(h(selectedIdx), 'Color', 'b', 'LineWidth', 1.8);
            uistack(h(selectedIdx), 'top');
        end
    end

    function toggleBad(~, ~)
        if isempty(selectedIdx)
            return;
        end
        cname = chanLabels{selectedIdx};
        if ismember(cname, sessionBadChans)
            % unmark
            sessionBadChans(strcmp(sessionBadChans, cname)) = [];
        else
            % mark
            sessionBadChans{end+1} = cname;
        end
        % update sessionBadIdx
        sessionBadIdx = find(ismember(chanLabels, sessionBadChans));
        % update listbox
        set(listBox, 'String', sessionBadChans);

        % update visuals: reset all then apply bad and selected styles
        for ii = 1:nChans
            set(h(ii), 'Color', [0.6 0.6 0.6], 'LineWidth', 1);
        end
        if ~isempty(sessionBadIdx)
            for ii = 1:numel(sessionBadIdx)
                set(h(sessionBadIdx(ii)), 'Color', 'r', 'LineWidth', 1.5);
            end
        end
        if ~isempty(selectedIdx)
            set(h(selectedIdx), 'Color', 'b', 'LineWidth', 1.8);
            uistack(h(selectedIdx), 'top');
        end
    end

    function closeOnly(~, ~)
        % Close GUI without removing channels
        userAction = 'close';
        uiresume(f);
        delete(f);
    end

    function closeAndRemove(~, ~)
        % Close GUI and remove bad channels using pop_select
        userAction = 'remove';
        uiresume(f);
        delete(f);
    end

    function close_noaction(~, ~)
        % invoked by window manager close button - treat as cancel
        userAction = '';
        uiresume(f);
        delete(f);
    end

end

