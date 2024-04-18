% Let us take the classification results and compare
% classification between ping and no-ping. If classification after the ping
% is better than classification after a pseudo-ping, this constitutes
% support for our main hypothesis.
%
% Select the category level which gave the best classification in the
% previous step (in no ping cond).
%
% Note: we continue with the binary classifier
%
% SvB

close all; clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 3; % 1 = top; 2 = mid; 3 = bot
zscore_opt = 1;
blcorr_opt = 0;
cue_or_ping_lock = 2; %lock analysis to cue(1) or ping(2)?
soa_opt = 2; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs!
             
% Stats params
nperm2 = 10^6; % 10^5 2nd level permutations (for FDR)
cluster_perms = 10^6; % 10^5 number of permutations for cluster stats 
stats_toi = [0 0.5]; % time window to perform stats over for pseudo-ping locked
cuestats_toi = [0.5 2]; % time window to perform stats over for cue-locked

% Folders
work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
analysis_path = [work_path,'analyse\results\'];
save_path  = [analysis_path];

% Functions
randpick   = @(x) x(randi(length(x))); % create function that randomly grabs one of the label vector

%% Load relevant class of results
load([save_path,'lda_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(blcorr_opt),'_z',num2str(zscore_opt)]);

% cue or ping lock?
if cue_or_ping_lock == 1
    dim1 = 4;
    dim2 = 3;
    toi = cuestats_toi;
elseif cue_or_ping_lock == 2
    dim1 = 2;
    dim2 = 1;
    toi = stats_toi;
end

tvec = curr_win{dim2}; 
tofilt = 1:numel(tvec)-1;
np_dat        = squeeze(bin(dim1,:,tofilt));
np_dat_s      = squeeze(bin_s(dim1,:,tofilt,:));
p_dat        = squeeze(bin(dim2,:,tofilt));
p_dat_s      = squeeze(bin_s(dim2,:,tofilt,:));

%% Go stats
% Get some parameters
winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
nsubj = size(np_dat,1);
nshuff = size(np_dat_s,3);
wins = cfg_c.toi_cue(1):winstep_t:cfg_c.toi_cue(end);

% Win time in samples
cfg_c.new_fs = 250; %hard-code
windur = windur_t.*cfg_c.new_fs;
winstep = winstep_t.*cfg_c.new_fs;

% % Filter participants with failed classifier attempts (cause of too few
% % trials)
% filtany   = @(x) x(any(x,2),:);
% np_dat = filtany(np_dat);
% np_dat_s = filtany(np_dat_s);
% p_dat = filtany(p_dat);
% p_dat_s = filtany(p_dat_s);

% Cut data
np_avg = squeeze(nanmean((np_dat),1)); % Average across participants
p_avg = squeeze(nanmean((p_dat),1));

startbin = nearest(wins,toi(1));
endbin = nearest(wins,toi(end))-1; % minus one to avoid errors
p_cut = p_dat(:,startbin:endbin,:);
np_cut = np_dat(:,startbin:endbin,:);

p_cut = p_cut'; % flip dimensions
np_cut = np_cut';

% Calculate SEMs
for w = 1:size(p_dat,2)  % For each window
ping_SEM(w,:) = std(p_dat(:,w))/sqrt(length(p_dat(:,w)));
noping_SEM(w,:) = std(np_dat(:,w))/sqrt(length(np_dat(:,w)));
end


%% Wilcoxon version
for w = 1:size(p_cut,1)  % For each window
    [wilc_pval(w),~] = ranksum(p_cut(w,:),np_cut(w,:));
end

% Prepare vector in which we slice out what we need
fdr2 = nan(size(np_avg,2),1);

% Franklin D. Roosevelt it
[~,~,~,wilc_pval] = fdr_bh(wilc_pval,0.05,'dep');

% Put it at the right moment for plotting (i.e. "2" means it's plottable)
fdr2(startbin:endbin) = wilc_pval;

% Mark significant bins
dat_sigind = find(fdr2<=0.05);

figure; hold on; % ADD MARKS IN OTHER FIGURE

% No ping
plot(tvec(1:end-1),np_avg,'LineWidth',4,'Color',[0.25 0.25 0.25]);
shadedErrorBar(tvec(1:end-1),np_avg,noping_SEM,'lineProps',{'LineWidth',0.01,'Color',[0.55 0.55 0.55]},'patchSaturation',0.15); hold on;

% Ping
plot(tvec(1:end-1),p_avg,'LineWidth',4,'Color',[0.89 0.75 0.31]); hold on;
shadedErrorBar(tvec(1:end-1),p_avg,ping_SEM,'lineProps',{'LineWidth',0.01,'Color',[0.89 0.75 0.31]},'patchSaturation',0.45); hold on;
xline(0);
title(['Ping vs no-ping']); hold on;
ylim([0.44 0.62]);
xlim([-0.5 2]);

p3 = plot(tvec(dat_sigind),p_avg(dat_sigind),'o','MarkerSize',4,'LineWidth',2,'color','b'); hold on;
try
    p3.MarkerFaceColor = [0.5 1 0];
    p3.MarkerSize = 4;
end
    legend({'Ping','No ping', 'Event onset'});
    
    
    %% Cluster version
% 
% [clus, pval, tval, pdis] = permutest(p_cut, np_cut, 1, 0.05, cluster_perms, 0);
% 
% % fix results
% for i = 1:numel(clus)
%     clus{i} = clus{i}(:) + startbin;
% end
% 
% % Find clusters
% pos_clust = find(pval < 0.05);
% clust_vec = [];
% 
% for i = 1:numel(pos_clust)
%     clust_vec = [clust_vec ; clus{i}];
% end
% 
% % PLOTTING NO PING

% figure; hold on; % ADD MARKS IN OTHER FIGURE
% 
% % No ping
% plot(tvec(1:end-1),np_avg,'LineWidth',4,'Color',[0.25 0.25 0.25]);
% shadedErrorBar(tvec(1:end-1),np_avg,noping_SEM,'lineProps',{'LineWidth',0.01,'Color',[0.55 0.55 0.55]},'patchSaturation',0.15); hold on;
% 
% % Ping
% plot(tvec(1:end-1),p_avg,'LineWidth',4,'Color',[0.89 0.75 0.31]); hold on;
% shadedErrorBar(tvec(1:end-1),p_avg,ping_SEM,'lineProps',{'LineWidth',0.01,'Color',[0.89 0.75 0.31]},'patchSaturation',0.45); hold on;
% xline(0);
% title(['Ping vs no-ping']); hold on;
% ylim([0.44 0.62]);
% xlim([-0.5 2]);
% % xlim([-0.5 1]);
% 
% p3 = plot(tvec(clust_vec),p_avg(clust_vec),'o','MarkerSize',4,'LineWidth',2,'color','b'); hold on;
% try
%     p3.MarkerFaceColor = [0.5 1 0];
%     p3.MarkerSize = 4;
% end
%     legend({'Ping','No ping', 'Event onset'});

