% REVISIONS CODE
%
% Update trial_info to aid decoding
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pplist = [1 3:6 8:15 17:22 24:33];

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\eeg_data\';
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\eeg_data\';

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
    load([eeg_path,sind,'_preproc_encoding']);
    
    enc_trlinfo = encoding_data_p.trialinfo;
    
    % Repair fourth column of all data
    enc_trlinfo(:,3) = base_examp_sort(encoding_data_p.trialinfo(:,1));
    
    % Step 6: Put the new trialinfos back in their box
    encoding_data_p.trialinfo = enc_trlinfo;
    
    nchans_start = size(encoding_data_p.trial{1},1);
    
    % Step 7: Double check to reinterpolate missing channels
    % Load layout
    if pp < 15
        load cap_old
    elseif pp > 14
        load cap_marios
    end

    enc = encoding_data_p;
        
    % detect and fix removed channels

    ms_enc = [];

    % CUE
    for i = 1:numel(lay.label)
        if sum(strcmp(lay.label(i),enc.label)) == 0
            ms_enc = [ms_enc ; i];
        end
    end
    
    % interpolate removed channels
    cfg = [];
    cfg.elec = elec;
    cfg.senstype = 'eeg';
    cfg.method = 'spline'; %http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=634
    cfg.neighbours = neighbours;
    

    % CUE
    if ~isempty(ms_enc)
        cfg.missingchannel = lay.label(ms_enc);
        encoding       = ft_channelrepair(cfg,enc);
    end
     
% OK now reord
enc_reord = encoding;
enc_reord = reordchan(enc_reord,lay);  
      
% Rereference to common average
cfg            = [];
cfg.channel    = {'all'};
cfg.demean     = 'yes';
cfg.reref      = 'yes';
cfg.implicitref = 'FCz';
cfg.refchannel = {'EEG'};
cfg.refmethod  = 'avg';
enc_reord_comm  = ft_preprocessing(cfg, enc_reord);


% Rereference using Laplacian
cfg = [];
cfg.method = 'spline';
cfg.elec = elec;
cfg.neighbours = neighbours;
enc_reord_lap  = ft_preprocessing(cfg, enc_reord);

nchans_end = size(enc_reord_lap.trial{1},1);

N_insert = nchans_end - nchans_start;
fprintf(['inserted: ',num2str(N_insert)]);

save([save_path,sind,'_reorder'],'enc_reord','enc_reord_comm','enc_reord_lap','ms_enc');
  
end
