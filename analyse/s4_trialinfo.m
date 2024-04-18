% For MemPing project; let us take the preprocessed data and retroactively
% add generate a trial info structure with trial number (i.e., the original
% trial index before trial removal) as well as subjective memory performance.
% This is only relevant for cue-locked data.
%
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pplist = [1 3:15 17:22 24:33];

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = [eeg_path];

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
    
    % Get preprocessed data
    load([eeg_path,sind,'_preproc_cuelock'],'cuelock_data_p');
    dat = eval('cuelock_data_p');

    % Get behavior
    behav_data = [behav_path,'behav_res_',sind,'.mat'];
    load(behav_data);
    disp(phs3_ans); % 1 = forgotten; 2 = remembered; 3 = no response
    
    % Find out which trials where included/removed
    % Thanks FieldTrip for your cfg.previous structures
    try
    triallist = dat.cfg.previous.previous{1}.previous.trials;
    catch error
    triallist = dat.cfg.previous{1,1}.previous.trials;
    end    
    
    % Get the behavior
    if pp > 5
        mem_resp = phs3_ans(triallist);
    else
        mem_resp = nan;
    end
    
    % Append these to trial structure
    dat.trialinfo(:,end+1) = mem_resp;
    dat.trialinfo(:,end+1) = triallist; 
    
    % OK now here's what's in trialinfo
    % For ping AND cue-lock:
    % 1st column = pinging SOA condition dichotomized (3 = no ping)
    % 2nd column = stimulus ID
    % 3rd column = stim category mid/top
    % 4th column = stim category bot
    % 5th column = pinging timing (99 = no-ping)
    % 6th column = memory response (1 = forgot; 2 = remember; 3 = no ans)
    % 7th column = original trial number before removing trials
    
    % For no_ping:
    % idem except:
    % 1st column = no-ping (SOA conditon) dichotomized
    % 5th column = no-ping SOA timing
   
   % Save updated trialinfo back into the preprocessed data files
        cuelock_trialinfo = dat.trialinfo;
        save([eeg_path,sind,'_preproc_cuelock'],'cuelock_trialinfo','-append');
end
