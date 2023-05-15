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
[STUDY, EEG] = pop_importbids(rawdata_path, ...
    'outputdir',fullfile(rawdata_path,'derivatives'),...
    'bidsevent','off',...
    'bidstask','task-eyesclosed');
STUDY = pop_statparams(STUDY, 'default');
[~,~,AvgChanlocs] = std_prepare_neighbors(STUDY, EEG, 'force', 'on');
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

parfor s=1:size(EEG,2)
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
        EEG(s) = pop_runica(EEG(s), 'icatype','picard','maxiter',500,'mode','standard','concatcond','on','options',{'pca',EEG(s).nbchan-1});
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
        % EEG(s) = eeg_regepochs(EEG(s),'recurrence',epoch_length * (1-epoch_overlap),...
        %     'limits',[0 epoch_length * (1-epoch_overlap)],'eventtype','epoch_start','extractepochs','off');
        % EEG(s) = pop_epoch(EEG(s),{'epoch_start'},[0 epoch_length],'epochinfo','yes');
        EEG(s) = pop_saveset(EEG(s),'savemode''resave')
    catch pipe_error
        error_report{s} = pipe_error.message; %#ok<SAGROW>
    end
end

% Save study
if exist('error_report','var')
    mask = cellfun(@(x) ~isempty(x), error_report); % which subject/session
    STUDY = std_rmdat(STUDY, EEG, 'datinds', find(mask));
    EEG(find(mask)) = [];
end
[STUDY,EEG] = pop_savestudy(STUDY,EEG,'savemode''resave');

% figure; pop_spectopo(ALLEEG(1), 1, [], 'EEG' , 'freq', [2 10 22], 'freqrange',[5 60],'electrodes','off');

%% connectivity [STUDY, EEG] = pop_loadstudy
% Compute source level FCM on Regions of interest with ROI connect.
vol       = fullfile(fileparts(which("eegplugin_dipfit")),['standard_BEM' filesep 'standard_vol.mat']);
mri       = fullfile(fileparts(which("eegplugin_dipfit")),['standard_BEM' filesep 'standard_mri.mat']);
elec      = fullfile(fileparts(which("eegplugin_dipfit")),['standard_BEM' filesep 'elec' filesep 'standard_1005.elc']);
leadfield = fullfile(fileparts(which("eegplugin_dipfit")),'LORETA-Talairach-BAs.mat');
parfor s=2:size(EEG,2)
    EEG(s)  = pop_dipfit_settings( EEG(s),'hdmfile',vol,'mrifile',mri,'chanfile',elec, ...
        'coordformat','MNI','coord_transform',[4.4114e-06 -1.4064e-05 7.5546e-06 1.2424e-07 -1.614e-07 -1.5708 1 1 1] ,'chansel',1:61);
    % EEG(s)  = pop_multifit( EEG(s),[],'threshold',100,'plotopt',{'normlen','on'});
    EEG(s)  = pop_leadfield( EEG(s),'sourcemodel',leadfield,'sourcemodel2mni',[],'downsample',1);
    EEG(s).roi = [];
    EEG(s) = pop_roi_activity(EEG(s),'model','LCMV','modelparams',{0.05},'atlas','Brain-Regions','nPCA',3);
    EEG(s) = pop_roi_connect(REEG(s),'morder',20,'naccu',[],'methods',{'CS','COH','wPLI','MIM','MIC'});
    EEG(s) = pop_saveset(EEG(s),'savemode''resave');
end

% convert to a csv file 
results = cell(60,3,)
for i = 1:numel(REEG)
    subject: 'sub-01'
    session

    csData = EEG(i).roi.CS;
    cohData = EEG(i).roi.COH;
    wpliData = EEG(i).roi.wPLI;
    mimData = EEG(i).roi.MIM;
    micData = EEG(i).roi.MIC;

    combinedData{i, 1} = csData;
    combinedData{i, 2} = cohData;
    combinedData{i, 3} = wpliData;
    combinedData{i, 4} = mimData;
    combinedData{i, 5} = micData;
end

combinedTable = cell2table(combinedData, 'VariableNames', {'CS', 'COH', 'wPLI', 'MIM', 'MIC'});
writetable(combinedTable, 'ROI.csv');

% export data
hermes = fullfile(rawdata_path,['derivatives' filesep 'HERMES']);
mkdir(hermes)
for s=1:size(EEG,2)
    tmp = EEG(s).data;
    [~,name] = fileparts(EEG(s).filename);
    save(fullfile(hermes,[name '.mat']),'tmp')
end


