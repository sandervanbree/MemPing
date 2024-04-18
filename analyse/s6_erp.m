% Let us take the preprocessed data and perform
% basic event-related potential style analyses and visualizations.
% This plots the ERP for cue-lock, ping, and no-ping locked data,
% collapsing all ping timings where relevant.
%
% Note: if you are interested in whole-brain topographies, run with all
% channels. If you want to plot visual ERPs, run the relevant channels.
%
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pplist = [1 3:6 8:15 17:22 24:33];
toi_cue = [-0.5 2.5]; % time range for plotting (cue condition)
toi_pnp = [-0.5 1]; % time range for plotting (ping no ping cond)
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
smoothbins = 20; %How many samples to apply a Gaussian smoothing kernel over?
soa_opt = 3; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs! 
             % For this script, this filters only ping-locked data
             
% Channels for ERP (make sure they are common to both caps!)
  common_chans = {'EEG','-EOG','-ECG','-AFZ','-FT9','-FT10'}; % ALL
% common_chans = {'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'Cz', 'Pz', 'Oz', 'CP1', 'CP2', 'C1',...
%      'C2' ,' P1', 'P2', 'CP3', 'CP4', 'PO3', 'PO4', 'PO7', 'PO8', 'CPz','POz'}; % VISUAL ONLY
 
work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
analysis_path = [work_path,'analyse\results\'];
save_path  = [analysis_path];

% Generate class-specific trigger sets
cat_trigs{1} = [0 1 2 3]; % everything
cat_trigs{2} = [0 1]; % object
cat_trigs{3} = [2 3]; % scene
cat_trigs{4} = [0];   % anim
cat_trigs{5} = [1];   % inanim
cat_trigs{6} = [2];   % ind
cat_trigs{7} = [3];   % outd
% 0 = object animate
% 1 = object inanimate
% 2 = scene indoor
% 3 = scene outdoor

% Bin edges for SOA
bin_edg = linspace(0.5,1.5,4); % create three bins

    % Start looping
    pp_ind = 1;
    for pp = pplist
        disp(['Working on participant ',num2str(pp)]);
        
        % Set stuff up
        if pp < 10
            sind = ['pp0',num2str(pp)];
        else
            sind = ['pp',num2str(pp)];
        end
        rng('default'); rng(pp,'twister');
        
        % Load data
        if ref_opt == 1
            append_txt = '_data_p';
        elseif ref_opt == 2
            append_txt = '_data_comm_p';
        elseif ref_opt == 3
            append_txt = '_data_lap_p';
        end
        
        
        load([eeg_path,sind,'_reorder'],['cuelock',append_txt]);
        dat_cue = eval(['cuelock',append_txt]);
        
        load([eeg_path,sind,'_reorder'],['ping',append_txt]);
        dat_ping = eval(['ping',append_txt]);
        
        load([eeg_path,sind,'_reorder'],['noping',append_txt]);
        dat_noping = eval(['noping',append_txt]);
        
        
        % Cap files; we have two cap files, one for each set of participants
        if pp < 15
            load cap_old
        elseif pp > 14
            load cap_marios
        end
        
        % Clean up a bit
        cfg                 = [];
        cfg.channel         = common_chans;                       % read all EEG channels except eye elecs
        cfg.detrend         = 'yes';
        cfg.hpfilter        = 'yes';                              % apply high pass filter for clean ERPs
        cfg.hpfreq          = 0.2;
        cfg.lpfilter        = 'yes';
        cfg.lpfilter        = 40;                                 % apply high pass filter for clean ERPs
        
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.2 0];
        dat_noping          = ft_preprocessing(cfg,dat_noping);
        dat_ping            = ft_preprocessing(cfg,dat_ping);
        dat_cue             = ft_preprocessing(cfg,dat_cue);

        % Filter to one ping SOA condition only, if needed
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

        %Fetch ping SOA
        ping_soa = dat_ping.trialinfo(:,5);
    
        % Extract ping locked trials
        cfg.trials = ping_soa > min_soa & ping_soa < max_soa;
        dat_ping = ft_selectdata(cfg,dat_ping);
                
        % sort out all conditions of interest
        for i = 1:numel(cat_trigs)
            trls_cue{i} = ismember(dat_cue.trialinfo(:,3),cat_trigs{i});
            trls_ping{i} = ismember(dat_ping.trialinfo(:,3),cat_trigs{i});
            trls_noping{i} = ismember(dat_noping.trialinfo(:,3),cat_trigs{i});
        end
        
        % extract data
        cfg = [];
        for i = 1:numel(cat_trigs)
            cfg.trials        = trls_cue{i};
            dat_sel_cue{i}    = ft_selectdata(cfg,dat_cue);
            
            cfg.trials        = trls_ping{i};
            dat_sel_ping{i}   = ft_selectdata(cfg,dat_ping);
            
            cfg.trials        = trls_noping{i};
            dat_sel_noping{i} = ft_selectdata(cfg,dat_noping);
        end
        
        % Timelock
        cfg               = [];
        for i = 1:numel(cat_trigs)
            cfg.latency                  = toi_cue;
            tl_cue{pp_ind,i}             = ft_timelockanalysis(cfg,dat_sel_cue{i});
            
            cfg.latency                   = toi_pnp;
            tl_ping{pp_ind,i}             = ft_timelockanalysis(cfg,dat_sel_ping{i});
            tl_noping{pp_ind,i}           = ft_timelockanalysis(cfg,dat_sel_noping{i});
        end
        
        %Average per individual
        all_indiv_cue(pp_ind,:)    = mean(tl_cue{pp_ind,1}.avg,1); % Fetch average across chans for this participant (all cats)
        all_indiv_ping(pp_ind,:)    = mean(tl_ping{pp_ind,1}.avg,1); % Fetch average across chans for this participant (all cats)
        all_indiv_noping(pp_ind,:)    = mean(tl_noping{pp_ind,1}.avg,1); % Fetch average across chans for this participant (all cats)
        
        pp_ind = pp_ind+1;
        clearvars dat_sel trls
    end
    
    % Grand average
    
    % CUE
    cfg = [];
    for i = 1:numel(cat_trigs)
        temp = tl_cue(:,i);
        ga_cue{i} = ft_timelockgrandaverage(cfg,temp{:}); % Get one grand average (GA) structure per category
        ga_3d_cue(i,:,:) = ga_cue{i}.avg; % Move out of FieldTrip structure (stimcat x chans x time)
        ga_2d_cue(i,:) = mean(ga_3d_cue(i,:,:),2); % Further collapse across channels
    end
    
    % PING
    for i = 1:numel(cat_trigs)
        temp = tl_ping(:,i);
        ga_ping{i} = ft_timelockgrandaverage(cfg,temp{:}); % Get one grand average (GA) structure per category
        ga_3d_ping(i,:,:) = ga_ping{i}.avg; % Move out of FieldTrip structure (stimcat x chans x time)
        ga_2d_ping(i,:) = mean(ga_3d_ping(i,:,:),2); % Further collapse across channels
    end
    
    % NO PING
    for i = 1:numel(cat_trigs)
        temp = tl_noping(:,i);
        ga_noping{i} = ft_timelockgrandaverage(cfg,temp{:}); % Get one grand average (GA) structure per category
        ga_3d_noping(i,:,:) = ga_noping{i}.avg; % Move out of FieldTrip structure (stimcat x chans x time)
        ga_2d_noping(i,:) = mean(ga_3d_noping(i,:,:),2); % Further collapse across channels
    end
    
    time_cue                   = tl_cue{1,1}.time;
    time_pnp                   = tl_ping{1,1}.time;
    

%% Plotting: Individual plots (only all categories)
% Remember the ordering across columns:
% 1 = everything
% 2 = object
% 3 = scene
% 4 = anim
% 5 = inanim
% 6 = ind
% 7 = outd
% The following analysis plots all trials (1)

conds = {'_cue','_ping','_noping'}; 
cond_opt = 2; % which cond to plot?

figure; hold on; % Just plot individual ERPs across all cats
if cond_opt == 1
    tvec = time_cue;
else
    tvec = time_pnp;
end

for pp = 1:numel(pplist)
    subplot(5,6,pp);
    toplot = strcat('all_indiv', conds(cond_opt), '(pp,:)');
    plot(tvec,smoothdata(eval(toplot{1}),'gaussian',smoothbins),'LineStyle','-','LineWidth',1.5,'Color',[0.33 0.03 0.33]);
    title(['pp' num2str(pplist(pp))]);
    xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
    yline(0,'LineWidth',0.5,'Color',[0.9 0.9 0.9 0.5],'LineStyle','-');
    xlim([tvec(1) tvec(end)]);
    xlabel('time [s]')
    ylabel('amplitude')
    
    if i == 1
        lgd = legend({'ERP','t=0','Zero'});
        lgd.FontSize = 10;
    end
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',10);
end

%% Plotting: Grand Average
% cat = 1; % Which dimension do you want to plot?
% % 1 = everything
% % 2 = object
% % 3 = scene
% % 4 = anim
% % 5 = inanim
% % 6 = ind
% % 7 = outd
%
% figure;hold on;
% title(['Cat: ',num2str(cat),' | Cond: ',num2str(cond), ' | Ref: ',num2str(ref_opt),' | Chans: ', num2str(chan_set)]);
% plot(time,smoothdata(ga_2d(cat,:),'gaussian',smoothbins),'LineStyle','-','LineWidth',4,'Color',[0.33 0.03 0.33]);
% xline(0,'LineWidth',3,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
% yline(0,'LineWidth',2.5,'Color',[0.9 0.9 0.9 0.5],'LineStyle','-');
% YL = get(gca, 'YLim');
% maxlim = max(abs(YL));
% set(gca, 'YLim', [-maxlim maxlim]);
%
% legend({'ERP','t = 0','amp = 0'});
% xlim([toi(1) toi(end)]);
% xlabel('time [s]')
% ylabel('amplitude')
% set(gca,'FontName','Arial');
% set(gca,'FontSize',16);

%% Plot time-resolved topography of ERP (for all)
% FIRST: CUE-LOCKED
load cap_marios

% define parameters for plotting
figure; hold on;
tvec = time_cue;
sample_count = length(time_cue);
sampling_rate = round(1/(time_cue(2)-time_cue(1)));

timestep      = 0.15; %(in seconds)
j = [tvec(1):timestep:tvec(end)];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

figure;
suptitle(['Cue-locked | Ref: ',num2str(ref_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(4,5,k);
    cfg.xlim       = [j(k) j(k+1)];
    cfg.zlim       = [min(ga_cue{1}.avg(:)).*0.6 max(ga_cue{1}.avg(:)).*0.6];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, ga_cue{1});
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end

% NOW PING vs NO PING
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
ga_pingvsnoping = ft_math(cfg,ga_ping{1},ga_noping{1});

tvec = time_pnp;
sample_count = length(time_pnp);
sampling_rate = round(1/(time_pnp(2)-time_pnp(1)));

timestep      = 0.15; %(in seconds)
j = [0:timestep:2];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

figure;
suptitle(['Red: Ping > No-ping | Ref: ',num2str(ref_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(5,5,k);
    cfg.xlim       = [j(k) j(k+1)];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.zlim       = [min(ga_pingvsnoping.avg(:)).*0.6 max(ga_pingvsnoping.avg(:)).*0.6];
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, ga_pingvsnoping);
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end

%% Stats (https://www.fieldtriptoolbox.org/tutorial/eventrelatedstatistics/)
% Perform topo ERPs with stats for all trials
cfg = [];
cfg.channel     = common_chans;
cfg.latency     = [0.2 0.4];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'bonferroni';
cfg.correcttail = 'prob';
cfg.numrandomization = 10^5;  

Nsub = numel(pplist);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

dim1 = tl_ping(:,1);
dim2 = tl_noping(:,1);

stat = ft_timelockstatistics(cfg, dim1{:}, dim2{:});

% make the plot
figure;
suptitle(['Red: Ping > No-ping | Ref: ',num2str(ref_opt)]);
cfg = [];
% plot
cfg.xlim       = [0.2 0.4];
cfg.comment    = 'xlim';
cfg.commentpos = 'title';
cfg.layout     = lay;
cfg.highlight = 'on';
cfg.highlightchannel = find(stat.mask);
cfg.highlightsymbol    = '*';
cfg.highlightcolor     = [0.1 0.8 0.1];
cfg.highlightsize      = 8;
cfg.colormap   = flipud(brewermap([],'RdBu'));
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, ga_pingvsnoping)
colorbar;
set(gca,'FontName','Arial');
set(gca,'FontSize',13);