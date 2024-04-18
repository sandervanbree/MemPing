% For MemPing project;
% Let us take the preprocessed, trialinfo-updated data and correct
% a few mistakes.
% 
% First correction: The trialinfo of all preprocessed data has an error; the 4th column
% is wrong (it assumes 12 types but there are 16). Correct the 4th column of trialinfo based
% on information derived in the 2nd column.
%
% Second correction: the randomized no-ping SOAs currently are not matched to stimulus index;
% unlike in the real pinging condition, here each stimulus has a range of
% SOAs associated with it. Fix this: one stimulus index <-> one random SOA.
%
% Third correction: Some of the trials in the ping and no-ping locked
% condition are erroneously from the other condition. Fix this by taking
% the cue-locked trialinfo as a ground truth and removing any ping trial
% from non-pinged dataset and any non-pinged trial from the pinged dataset.
%
% Fourth correction: reorder the channel to comport to the same order
% across participants, this improves ERP plotting
% 
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pplist = [1 3:6 8:15 17:22 24:33];

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = [eeg_path];

shuffle    = @(x) x(randperm(length(x))); % create function that shuffles label vector
randpick   = @(x) x(randi(length(x))); % create function that randomly grabs one of the label vector

% Generate base set of stimulus identities
    base_set = ones(1,12); % TWELVE TOKENS
    base_examp_sort = [];
    for i = 1:16 % 16 TYPES
        curr_set = base_set.*i;
        base_examp_sort = [base_examp_sort curr_set];
    end

% Start looping
for pp = pplist
    disp(['Working on participant ',num2str(pp)]);
    
    % Set stuff up
    if pp < 10
        sind = ['pp0',num2str(pp)];
    else
        sind = ['pp',num2str(pp)];
    end
    rng('default'); rng(pp,'twister');
    
    % Load EEG data
    load([eeg_path,sind,'_preproc_cuelock']);
    load([eeg_path,sind,'_preproc']);
    
    cue_trlinfo = cuelock_data_p.trialinfo;
    
    % First we check whether all the ping/no-ping locked trials really are
    % in the right dataset by crosschecking with cue-locked trialinfo
    noping_stims_ind = find(cue_trlinfo(:,1) == 3);
    temp_noping = cue_trlinfo(:,2);
    noping_stims = unique(temp_noping(noping_stims_ind)); % These are pinged stims
    
    ping_stims_ind = find(cue_trlinfo(:,1) < 3);
    temp_ping = cue_trlinfo(:,2);
    ping_stims = unique(temp_ping(ping_stims_ind)); % These are pinged stims
    
    % Check if non-overlapping
    errorcheck = sum(intersect(noping_stims,ping_stims)) == 0;
    if errorcheck ~= 1
        error('this should be non-overlapping');
    end
        
    % Filter no-ping locked trials to ONLY the right ones
    cfg = [];
    cfg.trials = ismember(noping_data_p.trialinfo(:,2),noping_stims);
    noping_data_p = ft_selectdata(cfg,noping_data_p);
    
    % Filter ping locked trials to ONLY the right ones
    cfg = [];
    cfg.trials = ismember(ping_data_p.trialinfo(:,2),ping_stims);
    ping_data_p = ft_selectdata(cfg,ping_data_p);
    
    % Correct cue-locked data's sixth column
    for i = 1:size(cuelock_data_p.trialinfo,1)
       if cuelock_data_p.trialinfo(i,1) < 3
           cuelock_data_p.trialinfo(i,6) = 0;
       end
    end
        
    % Repair fourth column of all data
    noping_data_p.trialinfo(:,4) = base_examp_sort(noping_data_p.trialinfo(:,2));

    cuelock_data_p.trialinfo(:,4) = base_examp_sort(cuelock_data_p.trialinfo(:,2));

    ping_data_p.trialinfo(:,4) = base_examp_sort(ping_data_p.trialinfo(:,2));
    
    % Fix noping
    %% Generate stimulus-to-pseudoSOA distribution using pinged trials
    soa_dist = unique(ping_data_p.trialinfo(:,end)); % what are the ping SOAs?
    noping_stim = unique(noping_data_p.trialinfo(:,2)); % what are the unpinged stimuli?
    noping_stim_shuff = shuffle(noping_stim); % shuffle just to allocate randomly next
    bin_edg = linspace(0.5,1.5,4); % create three bins

    % Step 1: Fetch the distribution of real pings
    for i = 1:numel(bin_edg)-1 % move through bins
        soa_bins{i} = soa_dist(find(soa_dist >= bin_edg(i) & soa_dist <= bin_edg(i+1))); % find within range
    end 
    
    % Step 2: Generate new list with no-pinged stimuli in column 1, and a
    % pseudo-SOA grabbed from the list of real ping SOas in column 2.
    
    triplecycle = repmat([1 2 3],[1 99]); % we need to equally grab from all three distribution
    pseudo_dist = zeros(numel(noping_stim_shuff),2);
    ind = 1;
    
    for i = 1:numel(noping_stim_shuff)
        pseudo_dist(ind,1) = noping_stim_shuff(ind); % stimulus of no-ping trial
        rand_soa = randpick(soa_bins{triplecycle(i)}); % random SOA
        
        while ismember(rand_soa,pseudo_dist(:,2)) % regrab if already grabbed
            rand_soa = randpick(soa_bins{triplecycle(i)});
        end
        
        pseudo_dist(ind,2) = rand_soa; % insert that SOA
        ind = ind+1;
    end
    
    % Step 3: Adapt the trial info of cue-locked data to incorporate these
    % pseudo-pings
    cue_trlinfo = cuelock_data_p.trialinfo;
    
    for i = 1:size(pseudo_dist,1)
        loc = find(pseudo_dist(i,1) == cue_trlinfo(:,2));
        cue_trlinfo(loc,6) = pseudo_dist(i,2); % Let's tack on another dimension so we don't override the ping or no-ping coding
    end
        
    % Step 4: Adapt the trial info of no ping-locked data
    np_trlinfo = noping_data_p.trialinfo;
    for i = 1:size(pseudo_dist,1)
        loc = find(pseudo_dist(i,1) == np_trlinfo(:,2));
        np_trlinfo(loc,5) = pseudo_dist(i,2);
    end   
    
    % Step 5: Fix dichotomy code
    for i = 1:size(np_trlinfo,1) % 
        if np_trlinfo(i,end) < bin_edg(2)
            np_trlinfo(i,1) = 0;
        elseif np_trlinfo(i,end) > bin_edg(2) && np_trlinfo(i,end) < bin_edg(3)
            np_trlinfo(i,1) = 1;
        elseif np_trlinfo(i,end) > bin_edg(3) && np_trlinfo(i,end) < bin_edg(4)
            np_trlinfo(i,1) = 2;
        elseif np_trlinfo(i,end) > bin_edg(4)
            np_trlinfo(i,1) = 3;
        end
    end
    
    % Step 6: Put the new trialinfos back in their box
    noping_data_p.trialinfo = np_trlinfo;
    cuelock_data_p.trialinfo = cue_trlinfo;
    
    % Step 7: Double check to reinterpolate missing channels
    % Load layout
    if pp < 15
        load cap_old
    elseif pp > 14
        load cap_marios
    end

    ping_data_p.trialinfo(:,end+1) = 1;
    noping_data_p.trialinfo(:,end+1) = 99;
    
    cfg = [];
    comb = ft_appenddata(cfg,ping_data_p,noping_data_p);
    cue = cuelock_data_p;
        
    % detect and fix removed channels
    ms_comb = [];
    ms_cue = [];
        
    % COMB
    for i = 1:numel(lay.label)
        if sum(strcmp(lay.label(i),comb.label)) == 0
            ms_comb = [ms_comb ; i];
        end
    end
    
    % CUE
    for i = 1:numel(lay.label)
        if sum(strcmp(lay.label(i),cue.label)) == 0
            ms_cue = [ms_cue ; i];
        end
    end
    
    % interpolate removed channels
    cfg = [];
    cfg.elec = elec;
    cfg.senstype = 'eeg';
    cfg.method = 'spline'; %http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=634
    cfg.neighbours = neighbours;
    
    % PING
    if ~isempty(ms_comb)
        cfg.missingchannel = lay.label(ms_comb);
        comb       = ft_channelrepair(cfg,comb);
    end
        
    % CUE
    if ~isempty(ms_cue)
        cfg.missingchannel = lay.label(ms_cue);
        cue       = ft_channelrepair(cfg,cue);
    end
     
% OK now reord
 comb_reord = comb;
 comb_reord = reordchan(comb_reord,lay);  
 cue_reord = reordchan(cue,lay);   
      
% Rereference to common average
cfg            = [];
cfg.channel    = {'all'};
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.implicitref = 'FCz';
cfg.refchannel = {'EEG'};
cfg.refmethod  = 'avg';
cue_reord_comm  = ft_preprocessing(cfg, cue_reord);
comb_reord_comm  = ft_preprocessing(cfg, comb_reord);

% Rereference using Laplacian
cfg = [];
cfg.method = 'spline';
cfg.elec = elec;
cfg.neighbours = neighbours;
cue_reord_lap  = ft_preprocessing(cfg, cue_reord);
comb_reord_lap  = ft_preprocessing(cfg, comb_reord);

% Tease apart
  cfg = [];
    cfg.trials    = find(comb_reord.trialinfo(:,end)==99);
    noping_data_p = ft_selectdata(cfg,comb_reord); % _p stands for preprocessed
    noping_data_comm_p = ft_selectdata(cfg,comb_reord_comm);
    noping_data_lap_p = ft_selectdata(cfg,comb_reord_lap);
    
    cfg.trials  = find(comb_reord.trialinfo(:,end)==1);
    ping_data_p = ft_selectdata(cfg,comb_reord); % _p stands for preprocessed
    ping_data_comm_p = ft_selectdata(cfg,comb_reord_comm);
    ping_data_lap_p = ft_selectdata(cfg,comb_reord_lap);

    cuelock_data_p = cue_reord;
    cuelock_data_comm_p = cue_reord_comm;
    cuelock_data_lap_p = cue_reord_lap;
    
    save([save_path,sind,'_reorder'],'noping_data_p','noping_data_comm_p','noping_data_lap_p'...
        ,'ping_data_p','ping_data_comm_p','ping_data_lap_p'...
        ,'cuelock_data_p','cuelock_data_comm_p','cuelock_data_lap_p',...
        'ms_comb','ms_cue');
    

end
