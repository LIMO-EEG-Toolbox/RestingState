% minimal resting state pipeline

% start eeglab and check plug-ins
close all;
clear variables;
rng('default');

% variables
rawdata_path  = '/indirect/staff/cyrilpernet/ds004148';
epoch_length  = 2; % seconds
epoch_overlap = 0.25; % 25% overlap allowed

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
    'subjects', 1:2,...
    'bidsevent','off',...
    'bidstask','task-eyesclosed');
ALLEEG = pop_select( ALLEEG,'nochannel',{'EXG1','EXG2','EXG3','EXG4','EXG5','EXG6','EXG7','EXG8', 'GSR1', 'GSR2', 'Erg1', 'Erg2', 'Resp', 'Plet', 'Temp'});
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

for s=1:size(ALLEEG,2)
    % downsample
    if ALLEEG(s).srate ~= 250
        ALLEEG(s) = pop_resample(ALLEEG(s), 250);
    end
    % 50Hz removal
    ALLEEG(s) = pop_zapline_plus(ALLEEG(s),'noisefreqs','line',...
        'coarseFreqDetectPowerDiff',4,'chunkLength',0,...
        'adaptiveNremove',1,'fixedNremove',1,'plotResults',0);
    % remove bad channels
    ALLEEG(s) = pop_clean_rawdata( ALLEEG(s),'FlatlineCriterion',5,'ChannelCriterion',0.87, ...
        'LineNoiseCriterion',4,'Highpass',[0.25 0.75] ,'BurstCriterion',20, ...
        'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian', ...
        'WindowCriterionTolerances',[-Inf 7]);
    % interpolate missing channels and reference 
    [~,idx] = setdiff({AvgChanlocs.expected_chanlocs.labels},{ALLEEG(s).chanlocs.labels});
    if ~isempty(idx)
        ALLEEG(s) = pop_interp(ALLEEG(s), AvgChanlocs.expected_chanlocs(idx), 'spherical');
    end
    ALLEEG(s) = pop_reref(ALLEEG(s),[],'interpchan','off');
    
    % ICA cleaning
    ALLEEG(s) = pop_runica(ALLEEG(s), 'icatype','picard','concatcond','on','options',{'pca',ALLEEG(s).nbchan-1});
    ALLEEG(s) = pop_iclabel(ALLEEG(s), 'default');
    ALLEEG(s) = pop_icflag(ALLEEG(s),[NaN NaN;0.8 1;0.8 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN]);
    ALLEEG(s) = pop_subcomp(ALLEEG(s),[],0);
    
    % clear data using ASR - just the bad segment
    ALLEEG(s) = pop_clean_rawdata(ALLEEG(s),'FlatlineCriterion','off','ChannelCriterion','off',...
        'LineNoiseCriterion','off','Highpass','off','BurstCriterion',20,...
        'WindowCriterion',0.25,'BurstRejection','on','Distance','Euclidian',...
        'WindowCriterionTolerances',[-Inf 7] );
    
    % epoching -- add random markers to create epochs we can use for power,
    % correlations, etc ...
    ALLEEG(s) = eeg_regepochs(ALLEEG(s),'recurrence',epoch_length * (1-epoch_overlap),...
        'limits',[0 epoch_length * (1-epoch_overlap)],'eventtype','epoch_start','extractepochs','off');
    ALLEEG(s) = pop_epoch(ALLEEG(s),{'epoch_start'},[0 epoch_length],'epochinfo','yes');
end

% Save study
STUDY = pop_savestudy(STUDY, ALLEEG, ...
    'filename', 'Resting_state', ...
    'filepath', fullfile(rawdata_path,'derivatives'));

% figure; pop_spectopo(ALLEEG(1), 1, [], 'EEG' , 'freq', [2 10 22], 'freqrange',[5 60],'electrodes','off');

% export data
hermes = fullfile(rawdata_path,['derivatives' filesep 'HERMES']);
mkdir(hermes)
for s=1:size(ALLEEG,2)
    tmp = ALLEEG(s).data;
    [~,name] = fileparts(ALLEEG(s).filename);
    save(fullfile(hermes,[name '.mat']),'tmp')
end

