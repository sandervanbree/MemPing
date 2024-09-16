% REVISIONS CODE
%
% Classification; train and test on encoding data
% SvB

clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Participant list
pp_list = [1 3:6 8:15 17:22 24:33];
% 
% pp_list = 1;

% Parameters
ref_opt = 3;    %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 1;    % 1 = top; 2 = mid; 3 = bot
decode_toi = [-0.5 2.5];

% Folders
work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\eeg_data\';
eeg_path2   = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\data\eeg_data\';
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\eeg_data\';

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
        append_txt = '_';
    elseif ref_opt == 2
        append_txt = '_comm';
    elseif ref_opt == 3
        append_txt = '_lap';
    end
    
    load([eeg_path,sind,'_reorder'],['enc_reord',append_txt]);
    enc = eval(['enc_reord',append_txt]);
    
    load([eeg_path2,sind,'_reorder'],['ping_data',append_txt,'_p']);
    ping = eval(['ping_data',append_txt,'_p']);
  

    % remove EOG and ECG
    cfg = [];
    cfg.channel = ft_channelselection({'all' '-ECG', '-EOG'}, enc.label);
    enc = ft_selectdata(cfg, enc);
    
    cfg.channel = ft_channelselection({'all' '-ECG', '-EOG'}, ping.label);
    ping = ft_selectdata(cfg, ping);
    
    %% Optional additional preprocessing
    cfg = [];
    if cfg_c.blcorr_opt == 1
        cfg.demean = 'yes';
        cfg.baselinewindow = cfg_c.blwin;
        enc = ft_preprocessing(cfg,enc);
        ping = ft_preprocessing(cfg,ping);
 
    end
    if cfg_c.hpf_freq > 0.05
        cfg        = [];
        cfg.hpfiltord = 4;
        cfg.hpfreq = cfg_c.hpf_freq;
        enc = ft_preprocessing(cfg,enc);
        ping = ft_preprocessing(cfg,ping);
    end
    
    %% Create clabels and filter
    temp_enc = enc.trialinfo(:,3);
    temp_ping = ping.trialinfo(:,3);
  
    %% Resample and timelock
    cfg = [];
    if ~isempty(cfg_c.new_fs)
        cfg.resamplefs = cfg_c.new_fs;
        enc = ft_resampledata(cfg,enc);
        ping = ft_resampledata(cfg,ret);
    end
    
    % Timelock them
    cfg.keeptrials = 'yes';
    enc_tl{ind} = ft_timelockanalysis(cfg,enc);
    ping_tl{ind} = ft_timelockanalysis(cfg,ping);
    
    %% Loop across exemplars
    enc_design_temp = []; % class label design matrix
    ping_design_temp = [];
    
    for xmp = 1:size(classes,1)
        
        % ENC
        cfg.trials = find(ismember(temp_enc,classes(xmp,:)));
         
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        enc_dat{ind,xmp} = ft_selectdata(cfg,enc_tl{ind});
        enc_design_temp = [enc_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
        % PING
         cfg.trials = find(ismember(temp_ping,classes(xmp,:)));
         
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        ping_dat{ind,xmp} = ft_selectdata(cfg,ping_tl{ind});
        ping_design_temp = [ping_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
    end
    
    %% Get the whole shabam
    cfg = [];
    
    enc_design{ind} = enc_design_temp;
    ping_design{ind} = ping_design_temp;
    
    ind = ind+1;
end

%% Classify
time_enc = enc_dat{1,1}.time;

winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
win_opt = cfg_c.win_opt;
nperm = cfg_c.nperm; % number of times the classifier is run with shuff labels
nrep = cfg_c.nrep;
nrep_p = cfg_c.nrep_p;
nsubj = numel(pp_list);

% Create windows
fs = round(1/(time_enc(2) - time_enc(1)));

wins = decode_toi(1):winstep_t:decode_toi(end);
wins_cue = decode_toi(1):winstep_t:decode_toi(end);

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
curr_time = {time_enc};
curr_win = {wins_cue};
curr_dat = {enc_dat};
curr_design = {enc_design};
bin = [];
bin_s = [];
bin_std = [];
bin_std_s = [];
multi = [];
multi_s = [];
multi_std = [];
multi_std_s = [];

mc = 1;
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
            toclass.label = enc_dat{s}.label;
            
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
            multi(mc,s,w) = 1;
            multi_std(mc,s,w) = 1;
            
               bin(mc,s,w) = res2.perf;
               bin_std(mc,s,w) = res2.perf_std;
%               bin(mc,s,w) = 1;
%               bin_std(mc,s,w) = 1;
             
            cfg.mvpa.repeat = nrep_p; % set separate nrep for perms
            for p = 1:nperm
                cfg.design = clab_shuff_cur(p,:); % enter in clabel
                
                % SHUFFLE LABELS
%                 res = ft_timelockstatistics(cfg, toclass);
%                 multi_s(mc,s,w,p) = res.accuracy;
%                 multi_std_s(mc,s,w,p) = res.accuracy_std;
                multi_s(mc,s,w,p) = 1;
                multi_std_s(mc,s,w,p) = 1;
                
                   [~, res2] = mv_classify(cfg2, toclass.trial, cfg.design);
                   bin_s(mc,s,w,p) = res2.perf;
                   bin_std_s(mc,s,w,p) = res2.perf_std;
%                   bin_s(mc,s,w,p) = 1;
%                   bin_std_s(mc,s,w,p) = 1;
            end
        end
    end

% SAVING
soa_opt = 3;
save([save_path,'lda_rev_enc_enc_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(cfg_c.blcorr_opt),'_z',num2str(cfg_c.zscore_opt)],...
    'cfg_c','multi','multi_std','bin','bin_std','multi_s','multi_std_s','bin_s','bin_std_s',...
    'curr_design','curr_win','curr_design','curr_time');
