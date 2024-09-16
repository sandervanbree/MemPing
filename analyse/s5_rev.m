% REVISIONS
%
% Plot classifier results
% 
% A few possibilities of decoding:
% Train on encoding and test on encoding (1)
% Train on encoding and test on non-pinged retrieval trials (2)
% Train on encoding and test on ping retrieval trials (3)
% Either cue or ping locked
%
% SvB

clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 1; % 1 = top; 2 = mid; 3 = bot
zscore_opt = 1;
blcorr_opt = 1;
decoding_opt    = 3; % see header
ping_vs_cuelock = 1; % ping lock (1) vs cue lock (2)
soa_opt = 3; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs!
bin_or_multi = 1; % binary (1) or multi (2) classifier?

% Stats params
nperm2 = 10^5; % 10^5 2nd level permutations (for FDR)
cluster_perms = 10^5; % 10^5 number of permutations for cluster stats 
stats_toi = [0 1]; % time window to perform stats over for pseudo-ping locked
cuestats_toi = [0.5 2]; % time window to perform stats over for cue-locked


% Naming regime
if ping_vs_cuelock == 1
    loadname2 = '_pinglock';
elseif ping_vs_cuelock == 2
    loadname2 = '_cuelock';
end

if decoding_opt == 1
    loadname = 'enc_enc';
    loadname2 = [];
elseif decoding_opt == 2
    loadname = 'enc_ret_noping';
elseif decoding_opt == 3
    loadname = 'enc_ret_ping';
end
% Folders
work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\eeg_data\';
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = eeg_path;

% Functions
randpick   = @(x) x(randi(length(x))); % create function that randomly grabs one of the label vector

%% Load relevant class of results
load([save_path,'lda_rev_',loadname,loadname2,'_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(blcorr_opt),'_z',num2str(zscore_opt)])

tvec = curr_win; 
toi = [-0.5 1];
toi_start = nearest(tvec,toi(1));
toi_end   = nearest(tvec,toi(end));

if bin_or_multi == 1
    dat        = squeeze(bin(1,:,toi_start:toi_end));
    dat_s      = squeeze(bin_s(1,:,toi_start:toi_end,:));
elseif bin_or_multi == 2
    dat        = squeeze(multi(1,:,toi_start:toi_end));
    dat_s      = squeeze(multi_s(1,:,toi_start:toi_end,:));
end

%% Go analysis
% Get some parameters

winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = toi(end)-toi(1);
nsubj = size(dat,1);
nshuff = size(dat_s,4);

wins = toi(1):winstep_t:toi(end);
cfg_c.new_fs = 250; %hard-code

% Win time in samples
windur = windur_t.*cfg_c.new_fs;
winstep = winstep_t.*cfg_c.new_fs;

dat_avg = squeeze(nanmean(dat,1)); % Average across participants 
dat_avg_p = squeeze(nanmean(squeeze(nanmean(dat_s,1)),2)); % [PERM] Average across participants 

% mean for wilcoxon
emp = dat;
emp_mn = mean(dat);
shuf = squeeze(nanmean(dat_s,3));

% Plot
for w = 1:numel(wins)  % For each window
    fprintf('-- Working on window %d --\n', w)
    
    if w == 1 % REUSE SHUFFLED CLASSIFIER LABELS
        for s = 1:nsubj
            for p2 = 1:nperm2
                shuff_i(s,p2) = randi(size(dat_s,3));
            end
        end
    end
    
    for p2 = 1:nperm2  % For each 2nd level perm
        for s = 1:nsubj % Grab a random shuff of each subject
            dat_temp(s) = dat_s(s,w,shuff_i(s,p2));
        end
        dat_sh(w,p2) = nanmean(dat_temp); % average across participants
    end
    dat_pval(w) = numel(find(dat_sh(w,:)>=emp_mn(w)))/nperm2; % p-value
    dat_errb_2nd(w,:) = ([prctile(dat_sh(w,:),5) prctile(dat_sh(w,:),95)]); % error bar
end


%% P-value stuff
% Prepare vector in which we slice out what we need
fdr2 = nan(numel(dat_pval),1);

% for ping/no-ping
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
res_p = dat_avg_p';

sigind = dat_sigind;
condlist = {'cue'};
tvec_plot = wins;

norm_perf = res-mean(res_p);
norm_perf_s = res_p-mean(res_p);

f1 = figure;hold on;
plot(tvec_plot,norm_perf,'LineWidth',4,'Color',[0.50 0.25 0.50]);
errb = (dat_errb_2nd(:,2)-dat_errb_2nd(:,1))./2;
shadedErrorBar(tvec_plot,norm_perf_s,errb,'lineProps',{'LineWidth',4,'Color',[0.7 0.7 0.7 0.3]}); hold on;
hold on;
xline(0);
yline(0,'linewidth',2);
 ylim([-0.05 0.1]);
 xlim([-0.5 2]);
xlabel('time [s]');
p1 = plot(tvec(sigind),norm_perf(sigind),'o','MarkerSize',6,'LineWidth',2,'color','r'); hold on;
try
    p1.MarkerFaceColor = [1 0.5 0];
    p1.MarkerSize = 3;
end
ylabel(cfg_c.metric{1});
l2 = legend({'empirical (HA)','shuffle (H0)','onset'});
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