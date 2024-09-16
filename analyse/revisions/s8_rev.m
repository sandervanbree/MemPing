% REVISIONS
%
% Plot stats of specifically encoding to encoding classifier results
%
% SvB

clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 1; % 1 = top; 2 = mid; 3 = bot
zscore_opt = 1;
blcorr_opt = 1; % 0 originally
soa_opt = 3; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs!
bin_or_multi = 1; % binary (1) or multi (2) classifier?

% Stats params
nperm2 = 10^6; % 10^5 2nd level permutations (for FDR)
cluster_perms = 10^6; % 10^5 number of permutations for cluster stats 
stats_toi = [0 2]; 
wins_range = [-0.5 2.5];

% Folders
work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'rev\data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
analysis_path = [work_path,'analyse\results\'];
save_path  = [analysis_path];

% Functions
randpick   = @(x) x(randi(length(x))); % create function that randomly grabs one of the label vector

%% Load relevant class of results
% load([save_path,'lda_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(blcorr_opt),'_z',num2str(zscore_opt)]);
load([eeg_path,'lda_rev_enc_c1_r3_soa3_b1_z1']);

mc = 1;
tvec = curr_win{mc}; 
tofilt = 1:numel(tvec)-1;
    dat        = squeeze(bin(mc,:,tofilt));
    dat_s      = squeeze(bin_s(mc,:,tofilt,:));


%% Go analysis
% Get some parameters
winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
nsubj = size(dat,1);
nshuff = size(dat_s,4);


    wins = wins_range(1):winstep_t:wins_range(end);
    
cfg_c.new_fs = 250; %hard-code

% Win time in samples
windur = windur_t.*cfg_c.new_fs;
winstep = winstep_t.*cfg_c.new_fs;

dat_avg = squeeze(nanmean(dat,1)); % Average across participants 
dat_avg_p = squeeze(nanmean(squeeze(nanmean(dat_s,1)),2)); % [PERM] Average across participants 


% mean for wilcoxon
emp = dat;
shuf = squeeze(nanmean(dat_s,3));

% Plot
for w = 1:numel(wins)-1  % For each window
    fprintf('-- Working on window %d --\n', w)
   
    [dat_pval(w),~] = ranksum(emp(:,w),shuf(:,w));
    
emp_SEM(w,:) = std(emp(:,w))/sqrt(length(emp(:,w)));
shuf_SEM(w,:) = std(shuf(:,w))/sqrt(length(shuf(:,w)));
end

%% P-value stuff
% Prepare vector in which we slice out what we need
fdr2 = nan(numel(dat_pval),1);

startbin = nearest(wins,stats_toi(1));
endbin = nearest(wins,stats_toi(end))-1; % minus one to avoid errors
pval_cut = dat_pval(startbin:endbin);

% Franklin D. Roosevelt it
[~,~,~,dat_fdr] = fdr_bh(pval_cut,0.05,'dep');

% Put it at the right moment for plotting (i.e. "2" means it's plottable)
fdr2(startbin:endbin) = dat_fdr;

% Mark significant bins
dat_sigind = find(fdr2<=0.05);

%% PLOTTING
res = dat_avg;
res_p = dat_avg_p;

sigind = dat_sigind;
condlist = {'ping','no ping','cue ping','cue no ping'};
tvec_plot = tvec(1:numel(res_p));

norm_perf = res-mean(res_p);

f1 = figure;hold on;
plot(tvec_plot,norm_perf,'LineWidth',4,'Color',[0.50 0.25 0.50]);
errb = [emp_SEM];
shadedErrorBar(tvec_plot,norm_perf,errb,'lineProps',{'LineWidth',0.01,'Color',[0.50 0.25 0.50]}); hold on;
title(condlist{mc}); hold on;
xline(0);
yline(0,'linewidth',2);
xlim([-0.5 2]);
% ylim([0.44 0.6]);
 ylim([-0.05 0.1]);

% xlabel('time [s]'); 
% p1 = plot(tvec_plot(sigind),norm_perf(sigind),'o','MarkerSize',6,'LineWidth',6,'color','r'); hold on;
% try
%     p1.MarkerFaceColor = [1 0.5 0];
%     p1.MarkerSize = 3;
% end

for i = 1:length(sigind)
    x = tvec_plot(sigind(i));
    y = norm_perf(sigind(i));
    line([x-0.01 x+0.01], [0.00 0.00], 'Color', 'r', 'LineWidth', 1.5); % Adjust the length of the line as needed
end

ylabel(cfg_c.metric{1});
l2 = legend({'Average'});
set(l2,'Location','NorthEast');
set(gca,'FontName','Arial');
set(gca,'FontSize',13); hold off;


%% Cluster version
% % NO PING
% [clus_np, p_np, t_np, pdist_np] = permutest(np_avg, np_avg_p', 1, 0.05, cluster_perms, 0);
% 
% % PLOTTING NO PING
% tvec = wins(1:end-1);   
% res = np_avg;
% 
% % Find clusters
% pos_clust = find(p_np < 0.05);
% clust_vec = [];
% for i = 1:numel(pos_clust)
%     clust_vec = [clust_vec clus_np{i}];
% end
% figure(f1); hold on; % ADD MARKS IN OTHER FIGURE
% p3 = plot(tvec(clust_vec),res(clust_vec),'o','MarkerSize',4,'LineWidth',2,'color','b'); hold on;
% try
%     p3.MarkerFaceColor = [0.5 1 0];
%     p3.MarkerSize = 4;
% end
% 
% % CUE
% [clus_cue, p_cue, t_cue, pdist_cue] = permutest(cue_avg, cue_avg_p', 1, 0.05, cluster_perms, 0);
% 
% % PLOTTING CUE
% tvec = wins_cue(1:end-1);
% res = cue_avg;
% 
% % Find clusters
% pos_clust = find(p_cue < 0.05);
% clust_vec = [];
% for i = 1:numel(pos_clust)
%     clust_vec = [clust_vec clus_cue{i}];
% end
% 
% figure(f2); hold on; % ADD MARKS IN OTHER FIGURE
% p4 = plot(tvec(clust_vec),res(clust_vec),'o','MarkerSize',4,'LineWidth',2,'color','b'); hold on;
% try
%     p4.MarkerFaceColor = [0.5 1 0];
%     p4.MarkerSize = 4;
% end
