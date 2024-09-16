% REVISIONS CODE
%
% Classify on encoding x retrieval but on POWER values '
% Only differences from main script s10 is replacement of value potentials with
% power values.
%
%
% SvB

clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Participant list
pp_list = [1 3:6 8:15 17:22 24:33];

% Parameters
ref_opt = 3;    %Reference set: 1 = default, 2 = common average, 3 = Laplacian
feature_vals = 2; % voltage potentials (1) vs EEG power at frequencies (2)
foi   = [3 7 ; 8 12 ; 13 30 ; 35 80]; % triggered in case of 2 above (theta (3-7 Hz), alpha (8-12 Hz), beta (13-30 Hz), and gamma (35-80 Hz)).
cat_lvl = 1;    % 1 = top; 2 = mid; 3 = bot
soa_opt = 3;    % 0 = 0.5-0.833;
                % 1 = 0.833 1.167;
                % 2 = 1.167 1.5;
                % 3 = ALL SOAs!
                
% Folders
work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
analysis_path = [work_path,'analyse\results\'];
save_path  = [analysis_path];

% Read function that fetches classify parameters
cat_opt = cat_lvl;
cfg_c = classify_params(cat_opt);

% Functions
shuffle    = @(x) x(randperm(length(x))); % create function that shuffles label vector
randpick   = @(x) x(randi(length(x))); % create function that randomly grabs one of the label vector

% Base functions
objanim = 1:4; % OBJ ANIM
objinanim = 5:8; % OBJ INANIM
sceind = 9:12; % SCE IND
sceoutd = 13:16; % SCE OUTD

% Categories
if cat_lvl == 1
    classes = [objanim objinanim ; sceind sceoutd];
elseif cat_lvl == 2
    classes = [objanim ; objinanim ; sceind ; sceoutd];
elseif cat_lvl == 3
    classes = (1:16)';
end

% Bin edges for SOA
bin_edg = linspace(0.5,1.5,4); % create three bins

%% 1: Load data from all participants and generate class labels
ind = 1;
for s = pp_list
    % Set stuff up
    if s < 10
        sind = ['pp0',num2str(s)];
    else
        sind = ['pp',num2str(s)];
    end
    
    disp(['preparing  ',sind]);
    
    % Load data
    if ref_opt == 1
        append_txt = '_data_p';
    elseif ref_opt == 2
        append_txt = '_data_comm_p';
    elseif ref_opt == 3
        append_txt = '_data_lap_p';
    end
    
    load([eeg_path,sind,'_reorder'],['cuelock',append_txt]);
    cue = eval(['cuelock',append_txt]);
    
    load([eeg_path,sind,'_reorder'],['ping',append_txt]);
    ping = eval(['ping',append_txt]);
    
    load([eeg_path,sind,'_reorder'],['noping',append_txt]);
    noping = eval(['noping',append_txt]);
    
    cfg = [];
    cfg.trials = find(cue.trialinfo(:,5) < 99); % Ping
    cue_ping = ft_selectdata(cfg,cue);    
    cfg.trials = find(cue.trialinfo(:,5) > 98); % No ping
    cue_noping = ft_selectdata(cfg,cue);
    
    % remove EOG and ECG
    cfg = [];
    cfg.channel = ft_channelselection({'all' '-ECG', '-EOG'}, noping.label);
    noping = ft_selectdata(cfg, noping);
    ping = ft_selectdata(cfg,ping);
    
    cfg.channel = ft_channelselection({'all' '-ECG', '-EOG'}, cue_ping.label);
    cue_ping = ft_selectdata(cfg, cue_ping);
    cue_noping = ft_selectdata(cfg, cue_noping);
    
    % Fetch soa info only
    cue_ping_soa = cue_ping.trialinfo(:,5);
    cue_noping_soa = cue_noping.trialinfo(:,6);
    ping_soa = ping.trialinfo(:,5);
    noping_soa = noping.trialinfo(:,5);
    
    % Test quickly whether things went right
    errorcheck = unique(cue_noping.trialinfo(:,5));
    if errorcheck ~= 99
        error('error');
    end
    
    % Filter based on SOA
    cfg = [];
    if soa_opt == 0
        min_soa = bin_edg(1);
        max_soa = bin_edg(2);
    elseif soa_opt == 1
        min_soa = bin_edg(2);
        max_soa = bin_edg(3);
    elseif soa_opt == 2
        min_soa = bin_edg(3);
        max_soa = bin_edg(4);
    elseif soa_opt == 3
        min_soa = -999;
        max_soa = 999;
    end

        % Extract cue-locked ping trials
        cfg.trials = cue_ping_soa > min_soa & cue_ping_soa < max_soa;
        cue_ping = ft_selectdata(cfg,cue_ping);
        
        % Extract ping locked trials
        cfg.trials = ping_soa > min_soa & ping_soa < max_soa;
        ping = ft_selectdata(cfg,ping);
    
    %% Optional additional preprocessing
    cfg = [];
    if cfg_c.blcorr_opt == 1
        cfg.demean = 'yes';
        cfg.baselinewindow = cfg_c.blwin;
        cue_ping = ft_preprocessing(cfg,cue_ping);
        cue_noping = ft_preprocessing(cfg,cue_noping);
    end
    if cfg_c.hpf_freq > 0.05
        cfg        = [];
        cfg.hpfiltord = 4;
        cfg.hpfreq = cfg_c.hpf_freq;
        cue_ping = ft_preprocessing(cfg,cue_ping);
        cue_noping = ft_preprocessing(cfg,cue_noping);
        noping = ft_preprocessing(cfg,noping);
        ping = ft_preprocessing(cfg,ping);
    end
    
    %% Create clabels and filter
    temp_cp = cue_ping.trialinfo(:,4);
    temp_cnp = cue_noping.trialinfo(:,4);
    temp_np = noping.trialinfo(:,4);
    temp_p = ping.trialinfo(:,4);
    
    %% Resample and timelock
    cfg = [];
    if ~isempty(cfg_c.new_fs)
        cfg.resamplefs = cfg_c.new_fs;
        ping = ft_resampledata(cfg,ping);
        noping = ft_resampledata(cfg,noping);
        cue_ping = ft_resampledata(cfg,cue_ping);
        cue_noping = ft_resampledata(cfg,cue_noping);
    end
    
    % Timelock them
    cfg.keeptrials = 'yes';
    ping_tl = ft_timelockanalysis(cfg,ping);
    noping_tl = ft_timelockanalysis(cfg,noping);
    cue_ping_tl = ft_timelockanalysis(cfg,cue_ping);
    cue_noping_tl = ft_timelockanalysis(cfg,cue_noping);
    
    % THIS IS NEW:
    % Freq as features if you wish
    if feature_vals == 2
        ping_f = [];
        noping_f = [];
        cue_ping_f = [];
        cue_noping_f = [];
        
        for f_fold = 1:size(foi,1)
            cfg           = [];
            cfg.output    = 'fourier';
            cfg.channel   = 'all';
            cfg.method    = 'hilbert';
            cfg.foi       = foi(f_fold,1):1:foi(f_fold,end);
            cfg.toi       = 'all';
            cfg.filttype  = 'fir';
            cfg.filtorder = (3/8)*512;
            cfg.filtdir   = 'twopass';
            cfg.width     = 2;
            cfg.pad       = 10;
            
            tmp              = ft_freqanalysis(cfg,ping_tl);
            tmp2             = abs(squeeze(mean(tmp.fourierspctrm,3))).^2;
            ping_f           = [ping_f , tmp2];
            
            tmp              = ft_freqanalysis(cfg,noping_tl);
            tmp2             = abs(squeeze(mean(tmp.fourierspctrm,3))).^2;
            noping_f           = [noping_f , tmp2];
            
            tmp              = ft_freqanalysis(cfg,cue_ping_tl);
            tmp2             = abs(squeeze(mean(tmp.fourierspctrm,3))).^2;
            cue_ping_f           = [cue_ping_f , tmp2];
            
            tmp              = ft_freqanalysis(cfg,cue_noping_tl);
            tmp2             = abs(squeeze(mean(tmp.fourierspctrm,3))).^2;
            cue_noping_f           = [cue_noping_f , tmp2];
            
        end
        
        C = cell(size(cue_noping_f, 2),1);
        for i = 1:size(C,1)
            C{i} = num2str(i);
        end
        
        ping_tl.trial = ping_f;
        ping_tl.label = C;
        
        noping_tl.trial = noping_f;
        noping_tl.label = C;

        cue_ping_tl.trial = cue_ping_f;
        cue_ping_tl.label = C;
        
        cue_noping_tl.trial = cue_noping_f;
        cue_noping_tl.label = C;
    end
    
    %% Loop across exemplars
    cp_design_temp = []; % class label design matrix
    cnp_design_temp = [];
    np_design_temp = [];
    p_design_temp = [];
    
    for xmp = 1:size(classes,1)
        
        % CUE NOPING
        cfg =[];
        cfg.trials = find(ismember(temp_cnp,classes(xmp,:)));
         
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        cnp_dat{ind,xmp} = ft_selectdata(cfg,cue_noping_tl);
        cnp_design_temp = [cnp_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
        % CUE PING
        cfg.trials = find(ismember(temp_cp,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        cp_dat{ind,xmp} = ft_selectdata(cfg,cue_ping_tl);
        cp_design_temp = [cp_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
        % NO PING
        cfg = [];
        cfg.trials = find(ismember(temp_np,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        np_dat{ind,xmp} = ft_selectdata(cfg,noping_tl);
        np_design_temp = [np_design_temp ; ones(numel(cfg.trials),1).*xmp];  
        
        % PING
        cfg.trials = find(ismember(temp_p,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        p_dat{ind,xmp} = ft_selectdata(cfg,ping_tl);
        p_design_temp = [p_design_temp ; ones(numel(cfg.trials),1).*xmp];
    end
    
    %% Get the whole shabam
    cfg = [];
    
    cp_design{ind} = cp_design_temp;
    cnp_design{ind} = cnp_design_temp;
    np_design{ind} = np_design_temp;
    p_design{ind} = p_design_temp;
    
    ind = ind+1;
end

%% Classify
time_np = np_dat{1,1}.time;
time_p = p_dat{1,1}.time;
time_cnp = cnp_dat{1,1}.time;
time_cp = cp_dat{1,1}.time;

winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
win_opt = cfg_c.win_opt;
nperm = cfg_c.nperm; % number of times the classifier is run with shuff labels
nrep = cfg_c.nrep;
nrep_p = cfg_c.nrep_p;
nsubj = numel(pp_list);

% Create windows
fs = round(1/(time_np(2) - time_np(1)));
wins = cfg_c.toi(1):winstep_t:cfg_c.toi(end);
wins_cue = cfg_c.toi_cue(1):winstep_t:cfg_c.toi_cue(end);

% Win time in samples
windur = windur_t.*fs;
winstep = winstep_t.*fs;

% times just for plotting
pl_wins = (wins(1:end-1)+winstep_t./2);
pl_wins_cue = (wins_cue(1:end-1)+winstep_t./2);

% setup classifier settings for the whole lot of them
cfg = [] ;
cfg.method          = 'mvpa';
cfg.avgovertime     = 'no';
cfg.features        = 'chan';
cfg.feedback        = 0;
cfg.mvpa            = [];
cfg.mvpa.classifier = 'multiclass_lda';
cfg.mvpa.cv         = 'kfold';
cfg.mvpa.repeat     = cfg_c.nrep;
cfg.mvpa.metric     = 'accuracy'; %this has to be accuracy; auc produces two values?
cfg.mvpa.k          = cfg_c.nk;
cfg.mvpa.feedback   = 'no';
if cfg_c.zscore_opt == 1
    cfg.mvpa.preprocess = 'zscore';
end

% setup classifier settings for the whole lot of them
cfg2 = [];
cfg2.metric      = cfg_c.metric;
cfg2.k           = cfg_c.nk;
cfg2.repeat      = cfg_c.nrep;
cfg2.feedback    = 0;
if cfg_c.zscore_opt == 1
    cfg2.preprocess = 'zscore';
end


%% Now go in sequence for:
% ping > no ping > cue ping > cue no-ping
curr_time = {time_p, time_np, time_cp, time_cnp};
curr_win = {wins wins wins_cue wins_cue};
curr_dat = {p_dat np_dat cp_dat cnp_dat};
curr_design = {p_design np_design cp_design cnp_design};
bin = [];
bin_s = [];
bin_std = [];
bin_std_s = [];
multi = [];
multi_s = [];
multi_std = [];
multi_std_s = [];

for mc = 1:4 % ping > no ping > cue ping > cue no-ping 
    for s = 1:nsubj
        clab_shuff    = [];
        for w = 1:numel(curr_win{mc})-1
            
            fprintf('-- cond %d subject %d @ %d pct --\n', mc,s,round(w./numel(curr_win{mc}).*100))
            fprintf('-- cond %d subject %d @ %d pct --\n', mc,s,round(w./numel(curr_win{mc}).*100))
            fprintf('-- cond %d subject %d @ %d pct --\n', mc,s,round(w./numel(curr_win{mc}).*100))
            
            % cut out the right data
            start_bin = nearest(curr_time{mc},curr_win{mc}(w)); % starting timepoint
            end_bin = start_bin+windur;
            
            % Cut, average, and append class-by-class
            cutdat = [];
            for c = 1:size(curr_dat{mc},2)
                tmp = curr_dat{mc}{s,c}.trial(:,:,start_bin:end_bin); % slice up
                
                if ~isempty(tmp)
                    tmp = apply_win(tmp,win_opt);
                else
                    tmp = [];
                end
                
                cutdat = [cutdat ; tmp];
            end
            
            % put in a structure
            toclass.trial = cutdat;
            toclass.time = 1;
            toclass.label = p_dat{s}.label;
            
            % class labels
            clabel = curr_design{mc}{s}; % Generate clabels
            
            % Ensure the same permutation across windows
            if w == 1
                clab_shuff = zeros(nperm,numel(clabel));
                for p = 1:nperm
                    clab_shuff(p,:) = shuffle(clabel);
                end
                clab_shuff_cur = clab_shuff;
            else
                clab_shuff_cur = clab_shuff;
            end
            
            cfg.design = clabel; % enter in clabel
%             res = ft_timelockstatistics(cfg, toclass);
               [~, res2] = mv_classify(cfg2, toclass.trial, cfg.design);
             
%             multi(mc,s,w) = res.accuracy;
%             multi_std(mc,s,w) = res.accuracy_std;
            
               bin(mc,s,w) = res2.perf;
               bin_std(mc,s,w) = res2.perf_std;
%                bin(mc,s,w) = 1;
%                bin_std(mc,s,w) = 1;
             
            cfg.mvpa.repeat = nrep_p; % set separate nrep for perms
            for p = 1:nperm
                cfg.design = clab_shuff_cur(p,:); % enter in clabel
                
                % SHUFFLE LABELS
%                 res = ft_timelockstatistics(cfg, toclass);
%                 multi_s(mc,s,w,p) = res.accuracy;
%                  multi_std_s(mc,s,w,p) = res.accuracy_std;
                
                   [~, res2] = mv_classify(cfg2, toclass.trial, cfg.design);
                   bin_s(mc,s,w,p) = res2.perf;
                   bin_std_s(mc,s,w,p) = res2.perf_std;
%                   bin_s(mc,s,w,p) = 1;
%                   bin_std_s(mc,s,w,p) = 1;
            end
        end
    end
end

% SAVING
save([save_path,'lda_rev_bl_powval_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(cfg_c.blcorr_opt),'_z',num2str(cfg_c.zscore_opt)],...
    'cfg_c','multi','multi_std','bin','bin_std','multi_s','multi_std_s','bin_s','bin_std_s',...
    'curr_design','curr_win','curr_design','curr_time');