% REVISIONS
%
% Compute variance in the data across conditions
% 
% SvB
clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pplist = [1 3:6 8:15 17:22 24:33];

min_trials = 30; % Only proceed if at least n trials have a "forgotten" response
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
soa_opt = 3; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs!

smoothbins  = 30; %How many samples to apply a Gaussian smoothing kernel over?
analysis_opt = 1; % variance across channels (1) or trials (2)
toi_plot    = [-0.3 2];
toi_stats    = [0.5 2];

% Channels for ERP that are common to both caps
   common_chans = {'EEG','-EOG','-ECG','-AFZ','-FT9','-FT10'};
%  common_chans = {'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'Cz', 'Pz', 'Oz', 'CP1', 'CP2', 'C1',...
%      'C2',' P1', 'P2', 'CP3', 'CP4', 'PO3', 'PO4', 'PO7', 'PO8', 'CPz', 'POz'};

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\data\eeg_data\';
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = eeg_path;

% Start looping
pp_ind = 1;
for pp = pplist
    disp(['Working on participant ',num2str(pp)]);
    
    if pp < 10
        sind = ['pp0',num2str(pp)];
    else
        sind = ['pp',num2str(pp)];
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
    
    % No-ping & ping-locked
        load([eeg_path,sind,'_reorder'],['cuelock',append_txt]);
        cue = eval(['cuelock',append_txt]);
        
        load([eeg_path,sind,'_reorder'],['ping',append_txt]);
        ping = eval(['ping',append_txt]);
        
        load([eeg_path,sind,'_reorder'],['noping',append_txt]);
        noping = eval(['noping',append_txt]);
        
        % Cap files; we have two cap files, one for each set of participants
        if pp < 15
            load cap_old
        elseif pp > 14
            load cap_marios
        end
        
    cfg = [];
    cfg.trials = find(cue.trialinfo(:,5) < 99); % Ping
    cue_ping = ft_selectdata(cfg,cue);    
    cfg.trials = find(cue.trialinfo(:,5) > 98); % No ping
    cue_noping = ft_selectdata(cfg,cue);
    
    % Fetch soa info only
    cue_ping_soa = cue_ping.trialinfo(:,5);
    cue_noping_soa = cue_noping.trialinfo(:,5);
    ping_soa = ping.trialinfo(:,5);
    noping_soa = noping.trialinfo(:,5);
    
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
    
    % Cap files; we have two cap files, one for each set of participants
    if pp < 15
        load cap_old
    elseif pp > 14
        load cap_marios
    end
    
%     % Clean up a bit
%     cfg                 = [];
%     cfg.channel         = common_chans;            % read all EEG channels except eye elecs
%     cfg.detrend         = 'yes';
%     cfg.hpfilter        = 'yes';                              % apply high pass filter for clean ERPs
%     cfg.hpfreq          = 0.2;
%     cfg.lpfilter        = 'yes';
%     cfg.lpfilter        = 40;                                 % apply high pass filter for clean ERPs
%     cfg.demean          = 'yes';
%     cfg.baselinewindow  = [-0.2 0];
%   
   cfg = [];
   cfg.keeptrials = 'yes';
  cfg.latency = toi_plot;
   cue_ping = ft_timelockanalysis(cfg,cue_ping);
   cue_noping = ft_timelockanalysis(cfg,cue_noping);
   
   for t = 1:size(cue_ping.trial,3) % Time is all the same so we can loop
       for tr = 1:size(cue_ping.trial,1) % across channels, no ping
           ping_fetch = squeeze(cue_ping.trial(tr,:,t));
           intervar_ping_ch(pp_ind,t) = var(ping_fetch);
       end
          
       for tr = 1:size(cue_noping.trial,1) % across channels, ping
           noping_fetch = squeeze(cue_noping.trial(tr,:,t));
           intervar_noping_ch(pp_ind,t) = var(noping_fetch);
       end
       
       for ch = 1:size(cue_ping.trial,2)
           ping_fetch = squeeze(cue_ping.trial(:,ch,t));
           intervar_ping_tr(pp_ind,t) = var(ping_fetch);   
       end
       
        for ch = 1:size(cue_noping.trial,2)
           noping_fetch = squeeze(cue_noping.trial(:,ch,t));
           intervar_noping_tr(pp_ind,t) = var(noping_fetch);   
       end
   end
   pp_ind = pp_ind+1;
end
   
tvec = cue_noping.time;

% Cut data
if analysis_opt == 1
np_dat = intervar_noping_ch;
p_dat = intervar_ping_ch;
elseif analysis_opt == 2
np_dat = intervar_noping_tr;
p_dat = intervar_ping_tr;
end

np_avg = squeeze(nanmean(np_dat,1));
p_avg = squeeze(nanmean(p_dat,1));

startbin = nearest(tvec,toi_stats(1));
endbin = nearest(tvec,toi_stats(end))-1; % minus one to avoid errors
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
[~,~,~,wilc_pval_corr] = fdr_bh(wilc_pval,0.05,'dep');

% Put it at the right moment for plotting (i.e. "2" means it's plottable)
fdr2(startbin:endbin) = wilc_pval;

% Mark significant bins
dat_sigind = find(fdr2<=0.05);

figure; hold on; % ADD MARKS IN OTHER FIGURE

% No ping
plot(tvec,np_avg,'LineWidth',4,'Color',[0.25 0.25 0.25]);
shadedErrorBar(tvec,np_avg,noping_SEM,'lineProps',{'LineWidth',0.01,'Color',[0.55 0.55 0.55]},'patchSaturation',0.15); hold on;

% Ping
plot(tvec,p_avg,'LineWidth',4,'Color',[0.89 0.75 0.31]); hold on;
shadedErrorBar(tvec,p_avg,ping_SEM,'lineProps',{'LineWidth',0.01,'Color',[0.89 0.75 0.31]},'patchSaturation',0.45); hold on;
xline(0);
title(['Ping vs no-ping']); hold on;

xlim([-0.5 2]);

p3 = plot(tvec(dat_sigind),p_avg(dat_sigind),'o','MarkerSize',4,'LineWidth',2,'color','b'); hold on;
try
    p3.MarkerFaceColor = [0.5 1 0];
    p3.MarkerSize = 4;
end
    legend({'Ping','No ping', 'Event onset'});

set(gca,'FontName','Arial');
set(gca,'FontSize',15);