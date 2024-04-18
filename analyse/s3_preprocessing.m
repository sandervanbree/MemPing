% For MemPing project; let us take the formatted data and apply
% preprocessing, participant by participant. We'll save every one
% individually, meaning we can for the most part recover what we've done.
%
% Note: from participants 1:14, we use the old caps with old layout.
% from participants 18:33, we use Marios' caps with Marios' layout.
% For participant 15:17, we use Marios' caps with old layout. This is
% accounted for automatically in the scripts:
%
% I have the right layout files for old and new caps. Hence, I just load
% the right layout files for 1:14 and 18:33. For 15:17, I first
% transform the data structure to accommodate Marios' layout.
%
% There are two preprocess modes: res_mode 1 and 2. Both need to be done
% separately.
%
% SvB
clear all; close all; clc;

res_mode = 1; %1 = timelocked ping/no-ping; 2 = timelocked ret cue

%% 0: Before starting
ft_defaults

% Parameters
pp = 1; % run once for every participant (for both res_modes)
sancheck = 0; %sanity check by inserting a sine wave into one electrode (to check layout and stuff)
channame = 'Fpz'; %which channel to sanity check?

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = [eeg_path];

% Set stuff up
if pp < 10
    sind = ['pp0',num2str(pp)];
else
    sind = ['pp',num2str(pp)];
end
rng('default'); rng(pp,'twister');

%% 1: Load data
% EEG data
load([save_path,sind,'_format']);
% Cap files; we have two cap files, one for each set of participants
if pp < 15
    load([dep_path,'cap_old']);
elseif pp > 14
    load([dep_path,'cap_marios']);
end

% We'll re-append for the moment to ensure preprocessing is done
% identically for both conditions
cfg = [];
cfg.keepsampleinfo = 'yes';
data_comb = ft_appenddata(cfg,noping_data,ping_data);
% the first 80 trials are pseudo-pings, the rest are pings. We'll tease
% them back apart later.

if res_mode == 1 % preprocess pinged and non-ping locked data
    data_curr = data_comb;
elseif res_mode == 2 % preprocess ret cue locked data
    data_curr = cuelock_data;
end

% For participant 1, amps seem to have been switched. Repair:
if pp == 1
    data_curr2 = data_curr;
    for tr = 1:numel(data_curr.trial)
        data_curr2.trial{tr} = [data_curr.trial{tr}(33:64,:) ; data_curr.trial{tr}(1:32,:)];
    end
    data_curr = data_curr2;
end

%% 2: Preprocessing: Initial filtering

% 2.1: For participant 15 and 17, fix channel list to accommodate Marios' layout
if pp >= 15 && pp <= 17
    data_curr.label = {'Fp1';'Fp2';'F3';'F4';'C3';'C4';'P3';'P4';'O1';'O2';'F7';...
        'F8';'T7';'T8';'P7';'P8';'Fz';'Cz';'Pz';'Oz';'FC1';'FC2';'CP1';'CP2';...
        'FC5';'FC6';'CP5';'CP6';'TP9';'TP10';'POz';'ECG';'F1';'F2';'C1';...
        'C2';'P1';'P2';'AF3';'AF4';'FC3';'FC4';'CP3';'CP4';'PO3';'PO4';...
        'F5';'F6';'C5';'C6';'P5';'P6';'AF7';'AF8';'FT7';'FT8';'TP7';'TP8'...
        ;'PO7';'PO8';'FT9';'FT10';'Fpz';'CPz'}; 
end

% SANITY CHECK ELECTRODE
if sancheck == 1
    t = linspace(0,1,size(data_curr.trial{1,1},2));
    chan_i = find(strcmp(channame,data_curr.label));
    for i = 1:numel(data_curr.trial)
        data_curr.trial{1,i}(chan_i,:) = sin(2*pi*60*t);
    end
end

% 2.2 Append EOG and ECG to the end
cfg              = [];
cfg.channel      = {'ECG','EOG'};
eyechans         = ft_selectdata(cfg, data_curr);
    
cfg.channel      = {'all', '-ECG', '-EOG'};
data_curr_noeye  = ft_selectdata(cfg, data_curr);

cfg              = [];
data_curr        = ft_appenddata(cfg, data_curr_noeye,eyechans);
 
% 2.2: Filter
cfg = [];
cfg.bsfilter = 'yes';
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.bsfreq = [49 51];
cfg.lpfreq = 80;
cfg.hpfiltord = 4;
cfg.hpfreq = 0.05;
data_filt = ft_preprocessing(cfg,data_curr);

% 2.3: Downsample
cfg = [];
cfg.resamplefs = 250;
cfg.detrend = 'no';
data_filt = ft_resampledata(cfg, data_filt);

% 2.4: Browse through data to prime user for preprocessing
cfg = [];
ft_databrowser(cfg,data_filt);

% 2.5: I like to get an ERP on the raw data, to get a rough clue which
% channels are bad visually
cfg = [];
cfg.channel = {'all', '-ECG', '-EOG'}; 
data_filt_tl_num = ft_timelockanalysis(cfg, data_filt);
data_filt_tl_lab = data_filt_tl_num;
for i = 1:numel(data_filt_tl_num.label)
    data_filt_tl_num.label{i,1} = num2str(i);
end
figure;
cfg = [];
lay2 = lay;
lay2.label(1:numel(data_filt_tl_num.label)) = data_filt_tl_num.label;
cfg.layout = lay2;
cfg.showlabels  = 'yes';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg,data_filt_tl_num);

figure;
cfg.layout = lay;
ft_multiplotER(cfg,data_filt_tl_lab);

%% 3: Preprocessing: Remove bad channels and trials
    cfg             = [];
    cfg.metric      = 'zvalue';  % use by default zvalue method
    data_clean       = ft_rejectvisual(cfg,data_filt);
    %NOTE 1: write down which channels and trials you've removed so you can
    %always re-do it 
    %NOTE 2: Remove trials before channels
    %NOTE 3: Ensure you don't remove EEG/ECG! 
    % For old cap, these are 31 (EOG) and 32 (ECG)
    % For the new cap, this is only 32 (ECG)
    %NOTE 4: Try different metrics
    
%% 4: Preprocessing: use independent component analysis to remove artefacts 
    cfg             = [];
    cfg.latency     = [0 1]; % Keep it to a second
    data_ica        = ft_selectdata(cfg,data_clean);
    
    % demean and detrend data
    cfg             = [];
    cfg.demean      = 'yes';
    cfg.detrend     = 'yes';
    cfg.channel     = {'all','-EOG','-ECG'};
    data_ica        = ft_preprocessing(cfg,data_ica);
    
    % run component analysis
    cfg                 = [];
    cfg.method          = 'runica';
    cfg.runica.maxsteps = 100;
    data_ica            = ft_componentanalysis(cfg,data_ica);
    
    % backproject to original data
    cfg             = [];
    cfg.unmixing    = data_ica.unmixing;
    cfg.topolabel   = data_ica.topolabel;
    data            = ft_componentanalysis(cfg, data_clean);
    
    cfg             = [];
    cfg.layout      = lay;
    cfg.channel     = [1:10]; % components to be plotted
    cfg.viewmode    = 'component';
    cfg.allowoverlap = 'yes';
    ft_databrowser(cfg, data_ica)
    
    % audio heads-up
    load handel.mat
    soundsc(y,50000)
    
    pause
    
    prompt          = {'components (write components to remove with space in between'};
    defaults        = {'xxx'};
    answer          = inputdlg(prompt, 'Components to remove', 1, defaults);
    [components]    = str2num(cell2mat(answer));
    % you don't have to write these down; components is saved at the end
    
    % Component to delete
    cfg                 = [];
    cfg.component       = components; % an array with the bad components
    data_new        = ft_rejectcomponent(cfg, data);
    data_new.reject = cfg.component; %save info of rejected components
    
%% 5: Preprocessing: repair bad channels by means of interpolation
% detect removed channels
missingchannel = [];
ind = 1;
for i = 1:numel(data_curr_noeye.label)
    if sum(strcmp(data_curr_noeye.label(i),data_new.label)) == 0
        missingchannel(ind) = i;
        ind = ind+1;
    end
end

% interpolate removed channels
if ~isempty(missingchannel)
    cfg = [];
    cfg.elec = elec;
    cfg.missingchannel = data_curr_noeye.label(missingchannel);
    cfg.senstype = 'eeg';
    cfg.method = 'spline'; %http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=634
    cfg.neighbours = neighbours;
    data_rep       = ft_channelrepair(cfg,data_new);
else
    data_rep = data_new;
end

mschans = data_curr.label(missingchannel);

% (6: Rereference to common average)
cfg            = [];
cfg.channel    = {'all'};
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.implicitref = 'FCz';
cfg.refchannel = {'EEG'};
cfg.refmethod  = 'avg';
data_rep_comm  = ft_preprocessing(cfg, data_rep);

% 7: Rereference using Laplacian
cfg = [];
cfg.method = 'spline';
cfg.elec = elec;
cfg.neighbours = neighbours;
data_rep_lap = ft_scalpcurrentdensity(cfg,data_rep);

%% 8: Check data visually (ERPs)
cfg = [];
data_rep_tl = ft_timelockanalysis(cfg, data_rep);
data_rep_comm_tl = ft_timelockanalysis(cfg, data_rep_comm);
data_rep_lap_tl = ft_timelockanalysis(cfg, data_rep_lap);

figure;
cfg = [];
cfg.layout = lay;
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg,data_rep_tl,data_rep_comm_tl);

figure;
ft_multiplotER(cfg,data_rep_lap_tl);

%% 9: Tease back apart into pseudo-ping and ping and save results
if res_mode == 1
    cfg = [];
    cfg.trials    = 1:numel(noping_data.trial);
    noping_data_p = ft_selectdata(cfg,data_rep); % _p stands for preprocessed
    noping_data_comm_p = ft_selectdata(cfg,data_rep_comm);
    noping_data_lap_p = ft_selectdata(cfg,data_rep_lap);
    
    cfg.trials  = numel(noping_data.trial)+1:numel(data_new.trial);
    ping_data_p = ft_selectdata(cfg,data_rep); % _p stands for preprocessed
    ping_data_comm_p = ft_selectdata(cfg,data_rep_comm);
    ping_data_lap_p = ft_selectdata(cfg,data_rep_lap);

    save([save_path,sind,'_preproc'],'data_ica','components','mschans','missingchannel','data_rep',...
    'ping_data_p','noping_data_p','ping_data_comm_p','noping_data_comm_p', 'ping_data_lap_p', 'noping_data_lap_p');

elseif res_mode == 2
    cuelock_data_p = data_rep;
    cuelock_data_comm_p = data_rep_comm;
    cuelock_data_lap_p = data_rep_lap;
    save([save_path,sind,'_preproc_cuelock'],'data_ica','components','mschans','data_rep','missingchannel',...
        'cuelock_data_p','cuelock_data_comm_p','cuelock_data_lap_p');
end