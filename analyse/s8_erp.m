% ERP check; here you can check the ERP for a specific stimulus level
% category, for a specific ping time (soa) option
%
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pp_list = [1 3:6 8:15 17:22 24:33];

min_trials = 30; % Only proceed if at least n trials have a "forgotten" response
toi_cue = [-0.5 2.5]; % time range for plotting (cue condition)
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 1; % 1 = top; 2 = mid; 3 = bot
soa_opt = 3; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs!
zscore_opt = 0; 

smoothbins = 30; %How many samples to apply a Gaussian smoothing kernel over?
saveorload = 1; % 1 = saving mode; 2 = loading + plotting mode

% Channels for ERP that are common to both caps
   common_chans = {'EEG','-EOG','-ECG','-AFZ','-FT9','-FT10'};
%  common_chans = {'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'Cz', 'Pz', 'Oz', 'CP1', 'CP2', 'C1',...
%      'C2',' P1', 'P2', 'CP3', 'CP4', 'PO3', 'PO4', 'PO7', 'PO8', 'CPz', 'POz'};

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
analysis_path = [work_path,'analyse\results\'];
save_path  = [analysis_path];

% Bin edges for SOA
bin_edg = linspace(0.5,1.5,4); % create three bins

% Base functions
objanim = 1:4; % OBJ ANIM
objinanim = 5:8; % OBJ INANIM
sceind = 9:12; % SCE IND
sceoutd = 13:16; % SCE OUTD

% Categories
if cat_lvl == 1
    classes = [objanim objinanim ; sceind sceoutd];
elseif cat_lvl ~= 1
    error('only obj/sce is implemented right now');
end

% Start looping
pp_ind = 1;
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
    
    % Read function that fetches classify parameters
    cat_opt = cat_lvl;
    cfg_c = classify_params(cat_opt);
    
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
    cfg.channel = common_chans;
    noping = ft_selectdata(cfg, noping);
    ping = ft_selectdata(cfg,ping);
    
    cfg.channel = common_chans;
    cue_ping = ft_selectdata(cfg, cue_ping);
    cue_noping = ft_selectdata(cfg, cue_noping);
    
    % Fetch soa info only
    cue_ping_soa = cue_ping.trialinfo(:,5);
    ping_soa = ping.trialinfo(:,5);
    
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
    
    cfg_c.blcorr_opt = 1;
    cfg_c.blwin = [0.5 0];
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
    ping_tl{ind} = ft_timelockanalysis(cfg,ping);
    noping_tl{ind} = ft_timelockanalysis(cfg,noping);
    cue_ping_tl{ind} = ft_timelockanalysis(cfg,cue_ping);
    cue_noping_tl{ind} = ft_timelockanalysis(cfg,cue_noping);
    

    %% Loop across exemplars
    cp_design_temp = []; % class label design matrix
    cnp_design_temp = [];
    np_design_temp = [];
    p_design_temp = [];
    
    for xmp = 1:size(classes,1)
        
        % CUE NOPING
        cfg.trials = find(ismember(temp_cnp,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        cnp_dat{ind,xmp} = ft_selectdata(cfg,cue_noping_tl{ind});
        cnp_design_temp = [cnp_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
        cfg = [];
        cnp_dat{ind,xmp} = ft_timelockanalysis(cfg,cnp_dat{ind,xmp});
        
        % CUE PING
        cfg.trials = find(ismember(temp_cp,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        cp_dat{ind,xmp} = ft_selectdata(cfg,cue_ping_tl{ind});
        cp_design_temp = [cp_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
        cfg = [];
        cp_dat{ind,xmp} = ft_timelockanalysis(cfg,cp_dat{ind,xmp});
        
        % NO PING
        cfg = [];
        cfg.trials = find(ismember(temp_np,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        np_dat{ind,xmp} = ft_selectdata(cfg,noping_tl{ind});
        np_design_temp = [np_design_temp ; ones(numel(cfg.trials),1).*xmp];  
        
        cfg = [];
        np_dat{ind,xmp} = ft_timelockanalysis(cfg,np_dat{ind,xmp});
        
        % PING
        cfg.trials = find(ismember(temp_p,classes(xmp,:)));
        
        if numel(cfg.trials) < 2 % If less than two classes, enter blank
            cfg.trials = [];
        end
        p_dat{ind,xmp} = ft_selectdata(cfg,ping_tl{ind});
        p_design_temp = [p_design_temp ; ones(numel(cfg.trials),1).*xmp];
        
        cfg = [];
        p_dat{ind,xmp} = ft_timelockanalysis(cfg,p_dat{ind,xmp});
    end
    
    %% Get the whole shabam
    cfg = [];
    
    cp_design{ind} = cp_design_temp;
    cnp_design{ind} = cnp_design_temp;
    np_design{ind} = np_design_temp;
    p_design{ind} = p_design_temp;
    
    ind = ind+1;
end

% Grand average
cfg = [];
p_ga_c1 = ft_timelockgrandaverage(cfg,p_dat{:,1});
np_ga_c1 = ft_timelockgrandaverage(cfg,np_dat{:,1});
cp_ga_c1 = ft_timelockgrandaverage(cfg,cp_dat{:,1});
cnp_ga_c1 = ft_timelockgrandaverage(cfg,cnp_dat{:,1});

p_ga_c2 = ft_timelockgrandaverage(cfg,p_dat{:,2});
np_ga_c2 = ft_timelockgrandaverage(cfg,np_dat{:,2});
cp_ga_c2 = ft_timelockgrandaverage(cfg,cp_dat{:,2});
cnp_ga_c2 = ft_timelockgrandaverage(cfg,cnp_dat{:,2});

% Difference ERP
% Cue-lock
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
p_diff = ft_math(cfg,p_ga_c1,p_ga_c2);
np_diff = ft_math(cfg,np_ga_c1,np_ga_c2);
cp_diff = ft_math(cfg,cp_ga_c1,cp_ga_c2);
cnp_diff = ft_math(cfg,cnp_ga_c1,cnp_ga_c2);

% All the same timevec so fetch one
tvec = p_dat{1,1}.time;

% Put into 1D
p_c1_1d = mean(p_ga_c1.avg,1);
np_c1_1d = mean(np_ga_c1.avg,1);
cp_c1_1d = mean(cp_ga_c1.avg,1);
cnp_c1_1d = mean(cnp_ga_c1.avg,1);

p_c2_1d = mean(p_ga_c2.avg,1);
np_c2_1d = mean(np_ga_c2.avg,1);
cp_c2_1d = mean(cp_ga_c2.avg,1);
cnp_c2_1d = mean(cnp_ga_c2.avg,1);

% Plot 1D results
% C1 & C2 Ping
figure; title(['Ping - obj & sce - ping-lock | SOA: ',num2str(soa_opt)]); hold on;
plot(tvec,smoothdata(p_c1_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.75 0.15 0.15]);
hold on;
plot(tvec,smoothdata(p_c2_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.15 0.15 0.75]);

xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
xlim([tvec(1) tvec(end)]);
% ylim([-3.0000e-04 4.0000e-04]);
% ylim([-5*(10^-3) 6*(10^-3)]);
% ylim([-0.02 0.035]);
   
xlabel('time [s]')
ylabel('amplitude')
lgd = legend({'obj','sce'});
lgd.FontSize = 10;

set(gca,'FontName','Arial');
set(gca,'FontSize',10);

% CUE LOCK
figure; title(['Ping - obj & sce - cue-lock | SOA: ',num2str(soa_opt)]); hold on;
plot(tvec,smoothdata(cp_c1_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.75 0.15 0.15]);
hold on;
plot(tvec,smoothdata(cp_c2_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.15 0.15 0.75]);

xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
xlim([tvec(1) tvec(end)]);
% ylim([-5*(10^-3) 6*(10^-3)]);
% ylim([-0.02 0.035]);

xlabel('time [s]')
ylabel('amplitude')
lgd = legend({'obj','sce'});
lgd.FontSize = 10;

set(gca,'FontName','Arial');
set(gca,'FontSize',10);

% C1 & C2 Ping
figure; title(['No Ping - obj & sce - pseudo-lock | SOA: ',num2str(soa_opt)]); hold on;
plot(tvec,smoothdata(np_c1_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.75 0.15 0.15]);
hold on;
plot(tvec,smoothdata(np_c2_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.15 0.15 0.75]);

xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
xlim([tvec(1) tvec(end)]);
% ylim([-5*(10^-3) 6*(10^-3)]);
% ylim([-0.02 0.035]);

   
xlabel('time [s]')
ylabel('amplitude')
lgd = legend({'obj','sce'});
lgd.FontSize = 10;

set(gca,'FontName','Arial');
set(gca,'FontSize',10);

% CUE LOCK
figure; title(['No Ping - obj & sce - cue-lock | SOA: ',num2str(soa_opt)]); hold on;
plot(tvec,smoothdata(cnp_c1_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.75 0.15 0.15]);
hold on;
plot(tvec,smoothdata(cnp_c2_1d),'LineStyle','-','LineWidth',1.5,'Color',[0.15 0.15 0.75]);

xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
xlim([tvec(1) tvec(end)]);
% ylim([-5*(10^-3) 6*(10^-3)]);
% ylim([-0.02 0.035]);

xlabel('time [s]')
ylabel('amplitude')
lgd = legend({'obj','sce'});
lgd.FontSize = 10;

set(gca,'FontName','Arial');
set(gca,'FontSize',10);

%% Topo stuff
% % Set up timelock stats
load cap_marios

% define parameters for plotting
time_cue = cue_ping.time{1};

figure; hold on;
tvec = time_cue;
sample_count = length(time_cue);
sampling_rate = 250;

timestep      = 0.15; %(in seconds)
j = [-0.5:timestep:2];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

figure;
suptitle(['Ping - obj & sce - ping-lock | SOA: ',num2str(soa_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(4,5,k);
    cfg.xlim       = [j(k) j(k+1)];
     cfg.zlim       = [-max(abs(p_diff.avg(:))).*0.6 max(abs(p_diff.avg(:))).*0.6];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, p_diff);
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end

figure;
suptitle(['Ping - obj & sce - cue-lock | SOA: ',num2str(soa_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(4,5,k);
    cfg.xlim       = [j(k) j(k+1)];
     cfg.zlim       = [-max(abs(cp_diff.avg(:))).*0.3 max(abs(cp_diff.avg(:))).*0.3];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, cp_diff);
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end

figure;
suptitle(['No Ping - obj & sce - pseudo-lock | SOA: ',num2str(soa_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(4,5,k);
    cfg.xlim       = [j(k) j(k+1)];
     cfg.zlim       = [-max(abs(np_diff.avg(:))).*0.3 max(abs(np_diff.avg(:))).*0.3];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, np_diff);
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end

figure;
suptitle(['No Ping - obj & sce - cue-lock | SOA: ',num2str(soa_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(4,5,k);
    cfg.xlim       = [j(k) j(k+1)];
     cfg.zlim       = [-max(abs(cnp_diff.avg(:))).*0.3 max(abs(cnp_diff.avg(:))).*0.3];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, cnp_diff);
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end

figure;
subplot(2,1,1);
plot(time_cue,cp_diff.avg); title('CUE-LOCK ping O-S');hold on;
subplot(2,1,2);
plot(time_cue,cnp_diff.avg); title('CUE-LOCK no ping O-S');

figure;
subplot(2,1,1);
plot(time_cue,p_diff.avg); title('PING-LOCK ping O-S');hold on;
subplot(2,1,2);
plot(time_cue,np_diff.avg); title('PING-LOCK no ping O-S');