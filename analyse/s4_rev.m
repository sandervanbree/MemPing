% REVISIONS CODE
%
% A few possibilities of decoding:
% Train on encoding and test on encoding (1) %% DEPRACATED, USE S4_REV_ENC_ENC
% Train on encoding and test on non-pinged retrieval trials (2)
% Train on encoding and test on ping retrieval trials (3)
% Either cue or ping locked
%
% SvB

clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Participant list
pp_list = [1 3:6 8:15 17:22 24:33];

% Parameters
ref_opt = 3;    %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 1;    % 1 = top; 2 = mid; 3 = bot
decoding_opt = 1; % see header
ping_vs_cuelock = 2; % ping lock (1) vs cue lock (2)
enc_toi  = [0.1 1]; % what time window in encoding phase determines classifier weights
test_toi_cuelock = [-0.5 2.8]; % what time period to evaluate/test performance in retrieval phase when CUE lock
test_toi_pinglock = [-0.5 1.3]; % and when PING lock option is chosen

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
    cfg = [];

    if decoding_opt < 3 % If test on non-pinged retrieval (or train and test on encoding; then we won't actually use "test" but avoid script issues)
        if ping_vs_cuelock == 1
            load([eeg_path2,sind,'_reorder'],['noping_data',append_txt,'_p']);
            test = eval(['noping_data',append_txt,'_p']);

        else
            load([eeg_path2,sind,'_reorder'],['cuelock_data',append_txt,'_p']);
            test = eval(['cuelock_data',append_txt,'_p']);
            
            if decoding_opt == 2
            cfg.trials = find(test.trialinfo(:,5) > 98); % No ping
            test = ft_selectdata(cfg,test);
            end
        end

    elseif decoding_opt == 3 % If test on ping retrieval
        if ping_vs_cuelock == 1
            load([eeg_path2,sind,'_reorder'],['ping_data',append_txt,'_p']);
            test = eval(['ping_data',append_txt,'_p']);
        else
            load([eeg_path2,sind,'_reorder'],['cuelock_data',append_txt,'_p']);
            test = eval(['cuelock_data',append_txt,'_p']);
            cfg.trials = find(test.trialinfo(:,5) < 99); % Ping
            test = ft_selectdata(cfg,test);
        end
    end

    % remove EOG and ECG
    cfg = [];
    cfg.channel = ft_channelselection({'all' '-ECG', '-EOG'}, enc.label);
    enc = ft_selectdata(cfg, enc);

    cfg.channel = ft_channelselection({'all' '-ECG', '-EOG'}, test.label);
    test = ft_selectdata(cfg, test);

    %% Optional additional preprocessing
    cfg = [];
    if cfg_c.blcorr_opt == 1
        cfg.demean = 'yes';
        cfg.baselinewindow = cfg_c.blwin;
        enc = ft_preprocessing(cfg,enc);
        test = ft_preprocessing(cfg,test);

    end
    if cfg_c.hpf_freq > 0.05
        cfg        = [];
        cfg.hpfiltord = 4;
        cfg.hpfreq = cfg_c.hpf_freq;
        enc = ft_preprocessing(cfg,enc);
        test = ft_preprocessing(cfg,test);
    end

    %% Create clabels and filter
    temp_enc = enc.trialinfo(:,3);
    temp_test = test.trialinfo(:,4);

    %% Resample and timelock
    cfg = [];
    if ~isempty(cfg_c.new_fs)
        cfg.resamplefs = cfg_c.new_fs;
        enc = ft_resampledata(cfg,enc);
        test = ft_resampledata(cfg,test);
    end

    % Timelock them
    cfg.keeptrials = 'yes';
    enc_tl{ind} = ft_timelockanalysis(cfg,enc);
    test_tl{ind} = ft_timelockanalysis(cfg,test);

    %% Loop across exemplars
    enc_design_temp = []; % class label design matrix
    test_design_temp = [];

    for xmp = 1:size(classes,1)

        % ENC
        cfg.trials = find(ismember(temp_enc,classes(xmp,:)));

        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        enc_dat{ind,xmp} = ft_selectdata(cfg,enc_tl{ind});
        enc_design_temp = [enc_design_temp ; ones(numel(cfg.trials),1).*xmp];

        % TEST (ping or no ping)
        cfg.trials = find(ismember(temp_test,classes(xmp,:)));

        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        test_dat{ind,xmp} = ft_selectdata(cfg,test_tl{ind});
        test_design_temp = [test_design_temp ; ones(numel(cfg.trials),1).*xmp];

    end

    %% Get the whole shabam
    cfg = [];

    enc_design{ind} = enc_design_temp;
    test_design{ind} = test_design_temp;

    ind = ind+1;
end

%% Classify
time_enc = enc_dat{1,1}.time;
time_test = test_dat{1,1}.time;

winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
win_opt = cfg_c.win_opt;
nperm = cfg_c.nperm; % number of times the classifier is run with shuff labels
nrep = cfg_c.nrep;
nrep_p = cfg_c.nrep_p;
nsubj = numel(pp_list);

% Create windows
fs = round(1/(time_test(2) - time_test(1)));

if ping_vs_cuelock == 1 % ping lock (1) vs cue lock (2)
    test_toi = test_toi_pinglock;
else
    test_toi = test_toi_cuelock;
end

wins_rev = test_toi(1):winstep_t:test_toi(end);

% Win time in samples
windur = windur_t.*fs;
winstep = winstep_t.*fs;


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
curr_time = {time_test};
curr_win = {wins_rev};
curr_dat = {test_dat};

bin = [];
bin_s = [];
bin_std = [];
bin_std_s = [];
multi = [];
multi_s = [];
multi_std = [];
multi_std_s = [];

mc = 1; % In REVISIONS, only use one mc condition; loops are handled manually
for s = 1:nsubj


    % For decoding opt 2 and 3, there is only one train window (use the full encoding window)
    clab_shuff    = [];

    cut_encdat = [];
    for c = 1:size(test_dat,2)
        start_bin_enc = nearest(time_enc,enc_toi(1));
        end_bin_enc = nearest(time_enc,enc_toi(end));

        temp = enc_dat{s,c}.trial(:,:,start_bin_enc:end_bin_enc); % slice up
        
        if ~isempty(tmp)
            tmp = apply_win(temp,win_opt);
        else
            tmp = [];
        end
        cut_encdat = [cut_encdat ; tmp];
    end


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
        totrain.trial = cut_encdat;
        totrain.time = 1;
        totrain.label = enc_dat{s}.label;

        totest.trial = cutdat;
        totest.time  = 1;
        totest.label = test_dat{s}.label;

        % class labels
        clabel_train = enc_design{s}; % Generate clabels
        clabel_test  = test_design{s};

        % Ensure the same permutation across windows
        if w == 1
            clab_shuff = zeros(nperm,numel(clabel_train));
            for p = 1:nperm
                clab_shuff(p,:) = shuffle(clabel_train);
            end
            clabel_train_shuff_curr = clab_shuff;
        else
            clabel_train_shuff_curr = clab_shuff;
        end
 
     
        % res = ft_timelockstatistics(cfg, toclass);
           if decoding_opt == 1
               [~, res2] = mv_classify(cfg2, totrain.trial, clabel_train);
           else
               [~, res2] = mv_classify(cfg2, totrain.trial, clabel_train, totest.trial, clabel_test);
           end

        % multi(mc,s,w) = res.accuracy;
        % multi_std(mc,s,w) = res.accuracy_std;
        multi(mc,s,w) = 1; % don't need for binary classifier
        multi_std(mc,s,w) = 1;

        bin(mc,s,w) = res2.perf;
        bin_std(mc,s,w) = res2.perf_std;
        %               bin(mc,s,w) = 1;
        %               bin_std(mc,s,w) = 1;

        for p = 1:nperm
            clabel_train_sh = clabel_train_shuff_curr(p,:); % enter in clabel

            % SHUFFLE LABELS

            % res = ft_timelockstatistics(cfg, toclass);
           if decoding_opt == 1
               [~, res2] = mv_classify(cfg2, totrain.trial, clabel_train);
           else
               [~, res2] = mv_classify(cfg2, totrain.trial, clabel_train, totest.trial, clabel_test);
           end

            % multi_s(mc,s,w,p) = res.accuracy;
            % multi_std_s(mc,s,w,p) = res.accuracy_std;
            multi_s(mc,s,w,p) = 1; % don't need for binary classifier
            multi_std_s(mc,s,w,p) = 1;

            bin_s(mc,s,w,p) = res2.perf;
            bin_std_s(mc,s,w,p) = res2.perf_std;
            %               bin_s(mc,s,w,p) = 1;
            %               bin_std_s(mc,s,w,p) = 1;
        end
    end
end

% SAVING
soa_opt = 3;

curr_design = enc_design;
curr_win = wins_rev;
curr_time = time_enc;

% Naming regime
if ping_vs_cuelock == 1
    savename2 = '_pinglock';
elseif ping_vs_cuelock == 2
    savename2 = '_cuelock';
end

if decoding_opt == 1
    savename = 'enc_enc';
    savename2 = [];
elseif decoding_opt == 2
    savename = 'enc_ret_noping';
elseif decoding_opt == 3
    savename = 'enc_ret_ping';
end


save([save_path,'lda_rev_',savename,savename2,'_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(cfg_c.blcorr_opt),'_z',num2str(cfg_c.zscore_opt)],...
    'cfg_c','multi','multi_std','bin','bin_std','multi_s','multi_std_s','bin_s','bin_std_s',...
    'curr_design','curr_win','curr_design','curr_time');