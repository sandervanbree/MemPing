% For MemPing project; let us take the classification results and compare
% classification peaks across ping conditions
%
% This analysis is cue-locked by default (because the question concerns
% whether classifier peaks arise according to their condition)
%
% SvB

clear all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pp_list = [1 3:6 8:15 17:22 24:33];
ref_opt = 3; %Reference set: 1 = default, 2 = common average, 3 = Laplacian
cat_lvl = 1; % 1 = top; 2 = mid; 3 = bot
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

% 1) classifier results
np_all  = [];
p_all   = [];
nps_all = [];
ps_all  = [];

for soa_opt = 0:2
    load([save_path,'lda_c',num2str(cat_lvl),'_r',num2str(ref_opt),'_soa',num2str(soa_opt),'_b',num2str(blcorr_opt),'_z',num2str(zscore_opt)]);

    np_all(soa_opt+1,:,:)        = bin(4,:,:);
    p_all(soa_opt+1,:,:)         = bin(3,:,:);
    
    nps_all(soa_opt+1,:,:,:)     = bin_s(4,:,:,:);
    ps_all(soa_opt+1,:,:,:)      = bin_s(3,:,:,:);
end

% 2) ERP results
ind = 1;
for s = pp_list
    % Set stuff up
    if s < 10
        sind = ['pp0',num2str(s)];
    else
        sind = ['pp',num2str(s)];
    end
    
    % Cue-locked
    load([eeg_path,sind,'_reorder'],['cuelock',append_txt]);
    cue = eval(['cuelock',append_txt]);
    
    
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
% Get some parameters
winstep_t = cfg_c.winstep_t; % fetch from config file
windur_t = cfg_c.windur_t;
cfg_c.new_fs = 250; %hard-code

% Get some more info
wins_cue   = cfg_c.toi_cue(1):winstep_t:cfg_c.toi_cue(end);
wins       = wins_cue;

toi_new    = [0.5 2];

% Cut Classifier
startbin_cl  = nearest(wins,toi_new(1));
endbin_cl    = nearest(wins,toi_new(end))-1; % minus one to avoid errors
tvec_cl   = wins(startbin_cl:endbin_cl);

np_all  = np_all(:,:,startbin_cl:endbin_cl);
p_all   = p_all(:,:,startbin_cl:endbin_cl);
nps_all = nps_all(:,:,startbin_cl:endbin_cl,:);
ps_all  = ps_all(:,:,startbin_cl:endbin_cl,:);

% Cut ERP
startbin_erp  = nearest(tvec_erp,toi_new(1));
endbin_erp    = nearest(tvec_erp,toi_new(end))-1; % minus one to avoid errors
tvec_erp      = tvec_erp(startbin_erp:endbin_erp);

p_erp         = p_erp(:,:,startbin_erp:endbin_erp);
np_erp        = np_erp(:,:,startbin_erp:endbin_erp);

% Plot the ERP
np_erp_mn = squeeze(mean(np_erp,2));
p_erp_mn = squeeze(mean(p_erp,2));

ef1 = figure;hold on;plot(tvec_erp,mean(np_erp_mn,1),'linewidth',3,'color',[0.25 0.25 0.25]);
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

% Plot the classifier
np_class_mn = squeeze(mean(np_all,2));
p_class_mn = squeeze(mean(p_all,2));
np_class_mn = squeeze(median(np_all,2,'omitnan'));
p_class_mn = squeeze(median(p_all,2,'omitnan'));

 cf1 = figure;hold on;plot(tvec_cl,mean(np_class_mn,1),'linewidth',3,'color',[0.25 0.25 0.25]);
% cf1 = figure;hold on;plot(tvec_cl,median(np_class_mn,1,'omitnan'),'linewidth',3,'color',[0.2 0.2 0.2]);
title('no ping classifier');legend;
xlabel('time [s]');
ylabel('perf');
set(gca,'FontName','Arial');
set(gca,'FontSize',13); hold off;

 cf2 = figure;hold on;plot(tvec_cl,p_class_mn,'linewidth',3);
% cf2 = figure;hold on;plot(tvec_cl,p_class_mn,'linewidth',3);
title('ping classifier');legend;
xlabel('time [s]');
ylabel('perf');
legend({'soa1','soa2','soa3'});
set(gca,'FontName','Arial');
set(gca,'FontSize',13); hold off;

%% Start analysis 1: PING PEAK DETECTION
% Params
smooth_opt = 1;
smoothbins = 4;
% Find maximum peak on a participant-by-participant and shuffle-by-shuffle

locs_emp = nan(3,size(p_all,2));
locs_shuff = nan(3,size(p_all,2),size(ps_all,4));

% basis
% PING
for soa = 1:3
    disp(['currently working on SOA ',num2str(soa)]);
    
    for pp = 1:size(p_all,2)
        curr_emp = squeeze(p_all(soa,pp,:));
        
        if smooth_opt == 1
            curr_emp = imgaussfilt(curr_emp,smoothbins);
        end
        
        csum = cumsum(curr_emp);
        csum_diff = diff(csum);
        ctrapz = cumtrapz(csum_diff);
        [a,b] = findpeaks(diff(ctrapz));
        [~,c] = max(a);
        
        if ~isempty(c)
            locs_emp(soa,pp) = b(c);
        end
                
        for s = 1:size(ps_all,4)
            
            curr_shuff = squeeze(ps_all(soa,pp,:,s));
            
            if smooth_opt == 1
                curr_shuff = imgaussfilt(curr_shuff,smoothbins);
            end
            
            csum = cumsum(curr_shuff);
            csum_diff = diff(csum);
            ctrapz = cumtrapz(csum_diff);
            [a,b] = findpeaks(diff(ctrapz));
            [~,c] = max(a);
            
            if ~isempty(c)
            locs_shuff(soa,pp,s) = b(c);
            end
        end
    end
end

emp_pk = nanmean(locs_emp,2);
shuff_pk = squeeze(nanmean(squeeze(nanmean(locs_shuff,3)),2));

% emp_pk = median(locs_emp,2,'omitnan');
% shuff_pk = squeeze(median(squeeze(median(locs_shuff,3,'omitnan')),2,'omitnan'));

emp_tv = tvec_cl(round(emp_pk));
shuff_tv = tvec_cl(round(shuff_pk));

cols = lines(3);

figure(cf2)
for soa = 1:3
    xline(emp_tv(soa),'color',cols(soa,:),'LineWidth',3);
end

% figure(cf2)
% for soa = 1:3
%     xline(shuff_tv(soa),'color',[0.2 0.2 0.2],'LineWidth',3);
% end

currsoa = 1;

cf3 = figure;hold on;plot(tvec_cl,squeeze(mean(p_all(currsoa,:,:),1)),'linewidth',1); hold on;
plot(tvec_cl,p_class_mn(currsoa,:),'color',[0.2 0.2 0.2],'linewidth',5);
title('ping classifier');legend;
ylim([0.4 0.7]);
xlabel('time [s]');
ylabel('perf');
set(gca,'FontName','Arial');
set(gca,'FontSize',13); hold off;
legend off


%% Start analysis 2: PING ORDER DETECTION
for pp = 1:numel(pp_list)
    for soa = 1:3
        temp1(soa)= locs_emp(soa,pp);
        for s = 1:size(ps_all,4)
            temp2(s,soa)= locs_shuff(soa,pp,s);
        end
    end
    
    [~, empord] = sort(temp1);
    
    for s = 1:size(ps_all,4)
        [~, shufford(s,:)] = sort(temp2(s,:));
    end
    
    for soa = 1:3
        empscore(soa) = abs(soa-empord(soa));
        
        for s = 1:size(ps_all,4)
            shuffscore(s,soa) = abs(soa-shufford(s,soa));
        end
    end
    
    dist_empscore(pp,:) = sum(empscore)./4; %scale to between 1 (perfect) and 0 (furthest away)
    
    for s = 1:size(ps_all,4)
    dist_shuffscore(pp,s,:) = sum(shuffscore(s,:))./4;
    end
end


% STATS time

nperm2 = 10^6;

for p2 = 1:nperm2
    
    if mod(p2,10000) == 0
        fprintf('-- Working on perm2 %d --\n', p2)
    end
    
    for pp = 1:numel(pp_list)
        sh(pp) = randpick(dist_shuffscore(pp,:));
    end
    ord_shuff(p2) = nanmean(sh);
end

pval = numel(find(nanmean(dist_empscore) >= ord_shuff(:)))/nperm2;
mn_emp = nanmean(dist_empscore);

% Plot results
% Create figure
figure;

% Create KDE using ksdensity
[f,xi] = ksdensity(ord_shuff);

% Plot KDE  
plot(xi,f,'LineWidth',2,'Color',[0.5 0.5 0.5]);

% Add vertical line at H1
hold on;
xline(mn_emp,'r','LineWidth',2); 

% Add p-value text box
text(0.35, 8, ['p = ' num2str(round(pval,3))], 'FontSize',23); 

% Set font to LaTeX 
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 14)


% Shaded area 
h = area([mn_emp, xi(xi>=mn_emp)], [0, f(xi>=mn_emp)]); 
set(h, 'FaceColor', [0.8 0.8 0.8]);

h = area([mn_emp, xi(xi<=mn_emp)], [0, f(xi<=mn_emp)]); 
set(h, 'FaceColor', [0.8 0.5 0.5]);


% Label axes  
xlabel('Error distance');
ylabel('Density');

% Add title
title('Distance error scores across participants');

% Tighten plot padding
set(gca, 'LooseInset', max(get(gca,'TightInset'), 0.02));

% Set legend
legend('Shuffle', 'Empirical', 'Location', 'NorthWest')