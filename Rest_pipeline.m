% minimal resting state pipeline

% start eeglab and check plug-ins
close all;
clear variables;
rng('default');

% variables
rawdata_path  = '/indirect/staff/cyrilpernet/ds004148';
epoch_length  = 2; % seconds
epoch_overlap = 0.25; % 25% overlap allowed
nsess         = 3; % 3 sessions per subject

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
if ~exist('pop_importbids','file')
    plugin_askinstall('bids-matlab-tools',[],1);
end

if ~exist('pop_loadbv','file')
    plugin_askinstall('bva-io',[],1);
end

if ~exist('pop_zapline_plus','file')
    plugin_askinstall('zapline-plus',[],1);
end

if ~exist('picard','file')
    plugin_askinstall('picard', 'picard', 1);
end

% import data and remove unwanted channels
[STUDY, ALLEEG] = pop_importbids(rawdata_path, ...
    'outputdir',fullfile(rawdata_path,'derivatives'),...
    'bidsevent','off',...
    'bidstask','task-eyesclosed');
STUDY = pop_statparams(STUDY, 'default');
[~,~,AvgChanlocs] = std_prepare_neighbors(STUDY, ALLEEG, 'force', 'on');
channel_info = AvgChanlocs.expected_chanlocs;
save(fullfile(rawdata_path,['derivatives' filesep 'channel_info.mat']),'channel_info')

figure('Name','Channel locations')
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized','outerposition',[0 0 1 1])
topoplot([],AvgChanlocs.expected_chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', ...
    'chaninfo',AvgChanlocs.expected_chanlocs);
title('EEG channels')
drawnow

% for each subject, downsample, clean 50Hz, remove bad channels,
% interpolate, re-reference to the average, run ICA to remove
% eye and muscle artefacts, delete bad segments

for s=1:size(EEG,2)
    try
        % downsample
        if EEG(s).srate ~= 250
            EEG(s) = pop_resample(EEG(s), 250);
        end
        % 50Hz removal
        EEG(s) = pop_zapline_plus(EEG(s),'noisefreqs','line',...
            'coarseFreqDetectPowerDiff',4,'chunkLength',0,...
            'adaptiveNremove',1,'fixedNremove',1,'plotResults',0);
        % remove bad channels
        EEG(s) = pop_clean_rawdata( EEG(s),'FlatlineCriterion',5,'ChannelCriterion',0.87, ...
            'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
            'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian', ...
            'WindowCriterionTolerances',[-Inf 7]);
        % interpolate missing channels and reference
        [~,idx] = setdiff({AvgChanlocs.expected_chanlocs.labels},{EEG(s).chanlocs.labels});
        if ~isempty(idx)
            EEG(s) = pop_interp(EEG(s), AvgChanlocs.expected_chanlocs(idx), 'spherical');
        end
        EEG(s) = pop_reref(EEG(s),[],'interpchan','off');
        
        % ICA cleaning
        EEG(s) = pop_runica(EEG(s), 'icatype','picard','concatcond','on','options',{'pca',EEG(s).nbchan-1});
        EEG(s) = pop_iclabel(EEG(s), 'default');
        EEG(s) = pop_icflag(EEG(s),[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
        EEG(s) = pop_subcomp(EEG(s),[],0);
        
        % clear data using ASR - just the bad segment
        EEG(s) = pop_clean_rawdata(EEG(s),'FlatlineCriterion','off','ChannelCriterion','off',...
            'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
            'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
            'WindowCriterionTolerances',[-Inf 7] );
        
        % epoching -- add random markers to create epochs we can use for power,
        % correlations, etc ...
        EEG(s) = eeg_regepochs(EEG(s),'recurrence',epoch_length * (1-epoch_overlap),...
            'limits',[0 epoch_length * (1-epoch_overlap)],'eventtype','epoch_start','extractepochs','off');
        EEG(s) = pop_epoch(EEG(s),{'epoch_start'},[0 epoch_length],'epochinfo','yes');
    catch pipe_error
        error_report{s} = pipe_error.message; %#ok<SAGROW>
    end
end

% Save study
if exist('error_report','var')
    mask = cellfun(@(x) ~isempty(x), error_report); % which subject/session
    EEG(mask)           = []; % delete from data structure
    % find subject names
    idx2 = 1;
    bad_sub = cell(1,length(find(mask)));
    for idx = find(mask)
        if ~isempty(STUDY.datasetinfo(idx).subject)
            bad_sub{idx2} = STUDY.datasetinfo(idx).subject;
            idx2=idx2+1;
        end
    end
    % delete from STUDY if all 3 sessions are bad
    bad_names = unique(bad_sub);
    for s=1:length(bad_names)
        if sum(strcmp(bad_names{s},bad_sub)) == nsess
            STUDY.subject(strcmp(bad_names{s},STUDY.subject)) = [];
        end
    end
    STUDY.datasetinfo(mask) = [];
end

STUDY = pop_savestudy(STUDY, EEG, ...
    'filename', 'Resting_state', ...
    'filepath', fullfile(rawdata_path,'derivatives'));

% figure; pop_spectopo(ALLEEG(1), 1, [], 'EEG' , 'freq', [2 10 22], 'freqrange',[5 60],'electrodes','off');

% export data
hermes = fullfile(rawdata_path,['derivatives' filesep 'HERMES']);
mkdir(hermes)
for s=1:size(EEG,2)
    tmp = EEG(s).data;
    [~,name] = fileparts(EEG(s).filename);
    save(fullfile(hermes,[name '.mat']),'tmp')
end

