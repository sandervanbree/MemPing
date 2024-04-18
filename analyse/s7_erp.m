% ERP check; similar as the previous script but here you can choose one
% specific ping timing condition (soa_opt).
% 
% SvB
clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pplist = [1 3:6 8:15 17:22 24:33];

min_trials = 30; % Only proceed if at least n trials have a "forgotten" response
toi_cue = [-0.5 2.5]; % time range for plotting (cue condition)
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
soa_opt = 0; % 0 = 0.5-0.833;
             % 1 = 0.833 1.167;
             % 2 = 1.167 1.5;
             % 3 = ALL SOAs!

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
    
    % Clean up a bit
    cfg                 = [];
    cfg.channel         = common_chans;            % read all EEG channels except eye elecs
    cfg.detrend         = 'yes';
    cfg.hpfilter        = 'yes';                              % apply high pass filter for clean ERPs
    cfg.hpfreq          = 0.2;
    cfg.lpfilter        = 'yes';
    cfg.lpfilter        = 40;                                 % apply high pass filter for clean ERPs
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.2 0];
        
   ping = ft_preprocessing(cfg,ping);
   noping = ft_preprocessing(cfg,noping);
   cue_ping = ft_preprocessing(cfg,cue_ping);
   cue_noping = ft_preprocessing(cfg,cue_noping);

   cfg = [];
   p{pp_ind} = ft_timelockanalysis(cfg,ping);
   np{pp_ind} = ft_timelockanalysis(cfg,noping);
   cp{pp_ind} = ft_timelockanalysis(cfg,cue_ping);
   cnp{pp_ind} = ft_timelockanalysis(cfg,cue_noping);
   
   pp_ind = pp_ind+1;
end

% Grand average
cfg = [];
p_ga = ft_timelockgrandaverage(cfg,p{:});
np_ga = ft_timelockgrandaverage(cfg,np{:});
cp_ga = ft_timelockgrandaverage(cfg,cp{:});
cnp_ga = ft_timelockgrandaverage(cfg,cnp{:});

% Difference ERP
% Cue-lock
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
p_np_ga = ft_math(cfg,p_ga,np_ga);

% Ping/pseudo lock
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
cp_cnp_ga = ft_math(cfg,cp_ga,cnp_ga);

% All the same timevec so fetch one
tvec = p{1}.time;

% Put into 1D
p_1d = mean(p_ga.avg,1);
np_1d = mean(np_ga.avg,1);
cp_1d = mean(cp_ga.avg,1);
cnp_1d = mean(cnp_ga.avg,1);

% Plot 1D results

% PING/NO-PING LOCK
figure; title(['Ping/no-ping locked ERP | SOA: ',num2str(soa_opt)]); hold on;
plot(tvec,smoothdata(p_1d,'gaussian',smoothbins),'LineStyle','-','LineWidth',1.5,'Color',[0.50 0.25 0.50]);
hold on;
plot(tvec,smoothdata(np_1d,'gaussian',smoothbins),'LineStyle','-','LineWidth',1.5,'Color',[0.25 0.50 0.25]);

xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
xlim([tvec(1) tvec(end)]);
% ylim([-3.0000e-04 4.0000e-04]);
% ylim([-1*(10^-4) 3.5*(10^-4)]);
   
xlabel('time [s]')
ylabel('amplitude')
lgd = legend({'ping','no-ping','ping/no-ping'});
lgd.FontSize = 10;

set(gca,'FontName','Arial');
set(gca,'FontSize',10);

% CUE LOCK
figure; title(['Cue-locked ERP | SOA: ',num2str(soa_opt)]); hold on;
plot(tvec,smoothdata(cp_1d,'gaussian',smoothbins),'LineStyle','-','LineWidth',1.5,'Color',[0.50 0.25 0.50]);
hold on;
plot(tvec,smoothdata(cnp_1d,'gaussian',smoothbins),'LineStyle','-','LineWidth',1.5,'Color',[0.25 0.50 0.25]);

xline(0,'LineWidth',1.5,'Color',[0.4 0.4 0.4 0.5],'LineStyle','-');
xlim([tvec(1) tvec(end)]);
xlabel('time [s]')
ylabel('amplitude')
lgd = legend({'ping','no-ping','cue'});
lgd.FontSize = 10;

set(gca,'FontName','Arial');
set(gca,'FontSize',10);

% % Set up timelock stats
load cap_marios

% define parameters for plotting
time_cue = cp{1}.time;

figure; hold on;
tvec = time_cue;
sample_count = length(time_cue);
sampling_rate = 250;

timestep      = 0.15; %(in seconds)
j = [-0.5:timestep:2];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples
% 
% figure;
% suptitle(['Cue-locked ping | SOA: ',num2str(soa_opt)]);
% cfg = [];
% % plot
% for k = 1:numel(j)-1
%     cfg.figure     = subplot(4,5,k);
%     cfg.xlim       = [j(k) j(k+1)];
%     cfg.zlim       = [-max(abs(cp_ga.avg(:))).*0.6 max(abs(cp_ga.avg(:))).*0.6];
%     cfg.comment    = 'xlim';
%     cfg.commentpos = 'title';
%     cfg.layout     = lay;
%     cfg.figure     = 'gca';
%     cfg.colormap   = flipud(brewermap([],'RdBu'));
%     ft_topoplotER(cfg, cp_ga);
%     
%     set(gca,'FontName','Arial');
%     set(gca,'FontSize',13);
% end
% 
% figure;
% suptitle(['Cue-locked pseudo-ping | SOA: ',num2str(soa_opt)]);
% cfg = [];
% % plot
% for k = 1:numel(j)-1
%     cfg.figure     = subplot(4,5,k);
%     cfg.xlim       = [j(k) j(k+1)];
%     cfg.zlim       = [-max(abs(cnp_ga.avg(:))).*0.3 max(abs(cnp_ga.avg(:))).*0.3];
%     cfg.comment    = 'xlim';
%     cfg.commentpos = 'title';
%     cfg.layout     = lay;
%     cfg.figure     = 'gca';
%     cfg.colormap   = flipud(brewermap([],'RdBu'));
%     ft_topoplotER(cfg, cnp_ga);
%     
%     set(gca,'FontName','Arial');
%     set(gca,'FontSize',13);
% end

figure;
suptitle(['Ping-locked ping | SOA: ',num2str(soa_opt)]);
cfg = [];
% plot
for k = 1:numel(j)-1
    cfg.figure     = subplot(4,5,k);
    cfg.xlim       = [j(k) j(k+1)];
    cfg.xlim       = [0.2 0.4];
%     cfg.zlim       = [-max(abs(p_ga.avg(:))).*0.3 max(abs(p_ga.avg(:))).*0.3];
    cfg.comment    = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout     = lay;
    cfg.figure     = 'gca';
    cfg.colormap   = flipud(brewermap([],'RdBu'));
    ft_topoplotER(cfg, p_ga);
    
    set(gca,'FontName','Arial');
    set(gca,'FontSize',13);
end
% 
% figure;
% suptitle(['Pseudo ping-locked no ping | SOA: ',num2str(soa_opt)]);
% cfg = [];
% % plot
% for k = 1:numel(j)-1
%     cfg.figure     = subplot(4,5,k);
%     cfg.xlim       = [j(k) j(k+1)];
%     cfg.zlim       = [-max(abs(np_ga.avg(:))).*0.3 max(abs(np_ga.avg(:))).*0.3];
%     cfg.comment    = 'xlim';
%     cfg.commentpos = 'title';
%     cfg.layout     = lay;
%     cfg.figure     = 'gca';
%     cfg.colormap   = flipud(brewermap([],'RdBu'));
%     ft_topoplotER(cfg, np_ga);
%     
%     set(gca,'FontName','Arial');
%     set(gca,'FontSize',13);
% end


