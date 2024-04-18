% Let us turn to the EEG data analysis. First,
% let us structure the data into a digestible format
% before we turn to preprocessing and analysis.
% In contrast to the behaviour, here we go participant by participant,
% one run at a time.
%
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pp = 1; % run once for every participant from 1 until the last (except 7 & 16).

% NOTE: As of 18, we swap caps; taken care of by scripts

checkdata_opt = 1; % [1] ft_databrowser over each output condition to see if trigs are
                   % where they ought to be

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = [eeg_path];

% Set stuff up
if pp < 10
    sind = ['pp0',num2str(pp)];
else
    sind = ['pp',num2str(pp)];
end
rng('default'); rng(pp,'twister'); % Set seed based on current participant

% Load data
dataset = [eeg_path,sind,'.vhdr'];
cfg = [];                               
cfg.dataset = dataset;
data_eeg    = ft_preprocessing(cfg);
fsample     = data_eeg.fsample;

%% 1: Retrieve behavioural results and triggers
load([behav_path,'behav_res_',sind,'.mat']);
disp(set_ret)

% Collapse stimulus info
clps = @(x) reshape(x.',1,[]);
soa_ret_dich_c = clps(soa_ret_dich);
soa_ret_c = clps(soa_ret);
set_ret_c = clps(set_ret);

% We also need a shuffled version of soa_ret_c to draw from for our
% non-ping condition
soa_ret_shuff_c = soa_ret_c(soa_ret_c < 98); % only pinged stims
soa_ret_shuff_c = soa_ret_shuff_c(randperm(length(soa_ret_shuff_c))); % shuffle it

% Make a dichotomized option just for reference sake
bin_edg = linspace(0.5,1.5,4);
for i = 1:numel(soa_ret_shuff_c)
    soa_temp = soa_ret_shuff_c(i);
    if soa_temp < bin_edg(2)
        soa_ret_shuff_dich_c(i) = 0;
    elseif soa_temp > bin_edg(2) && soa_temp < bin_edg(3)
        soa_ret_shuff_dich_c(i) = 1;
    elseif soa_temp > bin_edg(3) && soa_temp < bin_edg(4)
        soa_ret_shuff_dich_c(i) = 2;
    elseif soa_temp > bin_edg(4)
        soa_ret_shuff_dich_c(i) = 3;
    end
end

%% 2: Check "triggers" (event values) from data
cfg = [];
cfg.dataset             = dataset;
cfg.trialdef.eventtype  = '?';
dummy                   = ft_definetrial(cfg); % this shows all the triggers
 
% Here is the meaning of each trigger:
% 10 = start of encoding block (https://i.imgur.com/W5Wjqni.png)
% 11 = start of working memory test (https://i.imgur.com/IJsfY0N.png)
% 12 MISSING?? supposed to be start retrieval
% 13 = start of retrieval test (https://i.imgur.com/zvJKkZa.png)
% 14 = start of new block
% 20 = a stimulus (word of a word-image pair) during retrieval that will soon be pinged
% 21 = a stimulus (word of a word-image pair) during retrieval that will not be pinged
% 22 = the actual ping itself
% >22 = these are stimulus indices (e.g. 44-40 = stimulus 4) but they are not reliable
% so please ignore them.

%% 3: Get trigger info and basic checks
ind = 1;
for i = 2:size(dummy.event,2)
    if ~isempty(dummy.event(i).value)
        trig(ind) = sscanf(dummy.event(i).value,'S%d');
    else
        trig(ind) = 999;
    end
    samp(ind) = dummy.event(i).sample;
    
    ind = ind+1;
end

% check ping SOAs
ind = 1;
for i = 1:size(trig,2)-1
    if trig(i) == 20 && trig(i+1) == 22
        soa_check(ind) = (samp(i+1) - samp(i))./fsample;
        ind = ind+1;
    end
end

% plot results
bin_edg = linspace(0.5,1.5,4);
figure;histogram(soa_check,bin_edg); % should be roughly the same for each condition
title('Ping SOA distribution check');
xlabel('SOA [s]');
ylabel('frequency [n trials]');

% soa_check will simultaneously serve as a pool to draw random SOAs from
% for unpinged conditions (yielding pseudo-ping SOAs).

%% 4: Create trialinfo vector
% start with top and mid level;
base_set = ones(1,48);
base_set_sort = [base_set*0 base_set*1 base_set*2 base_set*3];
set_ret_dich_c = base_set_sort(set_ret_c); % 0 = object animate
                                          % 1 = object inanimate
                                          % 2 = scene indoor
                                          % 3 = scene outdoor

% move on to exemplars
base_set = ones(1,16);
base_examp_sort = [];
for i = 1:16 % create base set to grab from
    curr_set = base_set.*i;
base_examp_sort = [base_examp_sort ; curr_set];                                         
end
base_examp_sort = base_examp_sort';
base_examp_sort = base_examp_sort(:);
set_examp_c = base_examp_sort(set_ret_c');

% Append to data
soa_ret_dich_c; % 0 = 0.5 to 0.83
                % 1 = 0.83 to 1.16
                % 2 = 1.16 to 1.5
                % 3 = no ping
                
% One participant has a little oopsie in the data; fix automatically:
if pp == 13
    soa_ret_dich_c(164) = [];
    set_ret_dich_c(164) = [];
    set_ret_c(164) = [];
    set_examp_c(164) = [];
    soa_ret_c(164) = [];
end

% Generate trial infos that are unique to condition
trlinfo_ping = [soa_ret_dich_c; set_ret_c; set_ret_dich_c; set_examp_c'; soa_ret_c];
% in order: pinging condition dichotomized (3 = no-ping); stim ID; stim
% category mid/top; stim category bot; pinging condition (99 = no-ping); 
trlinfo_noping = [soa_ret_shuff_c; soa_ret_shuff_dich_c]; 
% no-ping condition; no-ping condition dichotomized

%% 5: Extract data structure of interest based on triggers
% Note: if you get an error for any of these, go into the respective
% trialfun and uncomment the "if ind == 319" code.

cfg = [];  
cfg.dataset             = dataset;
cfg.trig                = trig;
cfg.samp                = samp;
cfg.trlinfo_ping        = trlinfo_ping; %this...
cfg.trlinfo_noping      = trlinfo_noping; %...and this will be used for both trial info and trial definition
cfg.trialdef.prestim    = 1;     % Cut out some time before [s]
cfg.trialdef.poststim   = 3;     % and some time after [s]
cfg.trialfun            = 'trialdef_ping';
ping_cfg                = ft_definetrial(cfg);
ping_data               = ft_preprocessing(ping_cfg);

cfg.trialfun            = 'trialdef_noping';
noping_cfg              = ft_definetrial(cfg);
noping_data             = ft_preprocessing(noping_cfg);

cfg.trialfun            = 'trialdef_cuelock';
cuelock_cfg             = ft_definetrial(cfg);
cuelock_data            = ft_preprocessing(cuelock_cfg);

%% 6: Check data!
pause(0.01);
if checkdata_opt == 1
cfg = [];
ft_databrowser(cfg, ping_data) % at 0 should be ping, preceded by a cue S20 (which may not always be visible depending on prestim)
ft_databrowser(cfg, noping_data) % at 0 should be pseudo-ping [not labelled], preceded by a cue S21 ("")
ft_databrowser(cfg, cuelock_data); % at 0 should be a cue S20 or S21
end

%% 7: Save ping and no-ping data
save([save_path,sind,'_format'],'ping_data','noping_data','cuelock_data');
