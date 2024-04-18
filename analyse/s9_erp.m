% ERP check; here we visualize ERPs across all ping timings (soa) for
% cue-locked data.
% SvB

clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pp_list = [1 3:6 8:15 17:22 24:33];

ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
zscore_opt = 1;
blcorr_opt = 0;
append_txt = '_data_lap_p';
common_chans = {'C3', 'C4', 'P3', 'P4', 'O1', 'O2', 'Cz', 'Pz', 'Oz', 'CP1', 'CP2', 'C1',...
    'C2',' P1', 'P2', 'CP3', 'CP4', 'PO3', 'PO4', 'PO7', 'PO8', 'CPz', 'POz'};
% common_chans = 'all';

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

% 2) ERP results
ind = 1;
for s = pp_list
    % Set stuff up
    if s < 10
        sind = ['pp0',num2str(s)];
    else
        sind = ['pp',num2str(s)];
    end
    
    % Load data
        if ref_opt == 1
            append_txt = '_data_p';
        elseif ref_opt == 2
            append_txt = '_data_comm_p';
        elseif ref_opt == 3
            append_txt = '_data_lap_p';
        end
        
       % Cue-locked
    load([eeg_path,sind,'_reorder'],['cuelock',append_txt]);
    cue = eval(['cuelock',append_txt]);
    
    % Clean up a bit
    cfg                 = [];
    cfg.channel         = common_chans;            % read all EEG channels except eye elecs
    cfg.detrend         = 'yes';
    cfg.hpfilter        = 'no';                              % apply high pass filter for clean ERPs
    cfg.hpfreq          = 0.2;
    cfg.lpfilter        = 'yes';
    cfg.lpfilter        = 40;                                 % apply high pass filter for clean ERPs
    cfg.demean          = 'yes';
    cfg.baselinewindow  = [-0.2 0];
    cue                 = ft_preprocessing(cfg,cue);
    
    % Loop over SOAs
    for soa_opt = 0:2
        
        cfg.trials = find(cue.trialinfo(:,1) == soa_opt); % Ping
        temp_cp   = ft_selectdata(cfg,cue);
        
        cfg.trials = find(cue.trialinfo(:,1) == 3); % no Ping
        temp_cnp   = ft_selectdata(cfg,cue);
        
        cfg = [];
        cfg.keeptrials = 'yes';
        temp_cp   = ft_timelockanalysis(cfg,temp_cp);
        temp_cnp  = ft_timelockanalysis(cfg,temp_cnp);
        
        p_erp(soa_opt+1,ind,:) = squeeze(mean(squeeze(mean(temp_cp.trial,1)),1));
        np_erp(soa_opt+1,ind,:) = squeeze(mean(squeeze(mean(temp_cnp.trial,1)),1));
        
    end
    
    ind = ind+1;
    tvec_erp = cue.time{1};
end

%% Set up parameters
cat_opt = 1;
cfg_c = classify_params(cat_opt);

% Get some parameters
winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
cfg_c.new_fs = 250; %hard-code

% Get some more info
wins_cue   = cfg_c.toi_cue(1):winstep_t:cfg_c.toi_cue(end);
wins       = wins_cue;

toi_new    = [0 2];

% Cut Classifier
startbin_cl  = nearest(wins,toi_new(1));
endbin_cl    = nearest(wins,toi_new(end))-1; % minus one to avoid errors
tvec_cl   = wins(startbin_cl:endbin_cl);

% Cut ERP
startbin_erp  = nearest(tvec_erp,toi_new(1));
endbin_erp    = nearest(tvec_erp,toi_new(end))-1; % minus one to avoid errors
tvec_erp      = tvec_erp(startbin_erp:endbin_erp);

p_erp         = p_erp(:,:,startbin_erp:endbin_erp);
np_erp        = np_erp(:,:,startbin_erp:endbin_erp);

% Plot the ERP
np_erp_mn = squeeze(mean(np_erp,2));
p_erp_mn = squeeze(mean(p_erp,2));

ef1 = figure;hold on;plot(tvec_erp,mean(np_erp_mn,1),'linewidth',3,'color',[0.2 0.2 0.2]);
title('no ping ERP');legend;
xlabel('time [s]');
ylabel('amp');
set(gca,'FontName','Arial');
set(gca,'FontSize',13); hold off;

ef2 = figure;hold on;plot(tvec_erp,p_erp_mn,'linewidth',3);
title('ping ERP');legend;
xlabel('time [s]');
ylabel('amp');
legend({'soa1','soa2','soa3'});
set(gca,'FontName','Arial');
set(gca,'FontSize',13); hold off;
