function [cfg] = classify_params(cat_opt)
% These are parameters that are read by the classify_top, classify_mid, and
% classify_bot scripts. This means you don't have to copy over settings
% across scripts.

% Meta-parameters
cfg.new_fs = []; % what sampling rate to implement?
cfg.zscore_opt = 1; %1 should the window of data to-be-classified be z-scored?
cfg.blcorr_opt = 0; %0 Currently z-scoring only makes sense for cue locking
cfg.blwin = [-0.2 0]; % start to end time for baseline correction
cfg.win_opt = 1; % 0 = mean across window; 1 = gaussian-weighted mean across window
cfg.hpf_freq = 0.05; % The data was filtered at 0.05Hz to begin with, do you want additional hpf? Only activates if >0.05.

% Classification for LDA
cfg.metric = {'auc'};
cfg.windur_t = 0.14; %0.14; % window duration in s
cfg.winstep_t = 0.02; %0.02; % by how much time should it move every step
cfg.toi = [-0.5 1]; % time window to analyze for ping/no-ping
cfg.toi_cue = [-0.5 2.8]; % time window to analyze for cue lock

if cat_opt == 1 % TOP (obj sce)
    cfg.nperm = 10; % 50 number of times the classifier is run with shuff labels
    cfg.nk = 5; % 5 number of folds for k-fold crossvalidation
    cfg.nrep = 5; %20 number of repetitions within a classifier run
    cfg.nrep_p = 2; %2 number of repetitions within a shuffled classifier run
    
elseif cat_opt == 2 % MID (anim vs. inanim & ind vs outd & anim vs ind & inanim outd)
    cfg.nperm = 10; %100
    cfg.nk = 5; %5
    cfg.nrep = 5; %10
    cfg.nrep_p = 2; %2
    
elseif cat_opt == 3 % BOT (all examplars)
    cfg.nperm = 200; %50
    cfg.nk = 5; %2
    cfg.nrep = 50; %10
    cfg.nrep_p = 4; %2
end