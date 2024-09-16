% REVISIONS CODE
%
% Structure data into digestible format; this time for ENCODING data
% because we will run a decoder using the encoding
%
% SvB
clear all; close all; clc;

%% 0: Before starting
ft_defaults

pp_list = [1 3:6 8:15 17:22 24:33];

% Parameters
for pp = pp_list(1:end)

    % NOTE: As of 18, we swap caps; taken care of by scripts

checkdata_opt = 0; % [1] ft_databrowser over each output condition to see if trigs are
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
load(['\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\behav_data\behav_res',sind]);
disp(set_enc)

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

%% 4: Create trialinfo vector
% start with top and mid level;
% Collapse stimulus info
clps = @(x) reshape(x.',1,[]);
set_enc_c = clps(set_enc);

base_set = ones(1,48);
base_set_sort = [base_set*0 base_set*1 base_set*2 base_set*3];
set_enc_dich_c = base_set_sort(set_enc_c); % 0 = object animate
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
set_examp_c = base_examp_sort(set_enc_c');


%% 5: Extract data structure of interest based on triggers
% Note: if you get an error for any of these, go into the respective
% trialfun and uncomment the "if ind == 319" code.

cfg = [];  
cfg.dataset             = dataset;
cfg.trig                = trig;
cfg.samp                = samp;
cfg.encinfo             = [set_enc_c ;set_enc_dich_c];
cfg.trialdef.prestim    = 1;     % Cut out some time before [s]
cfg.trialdef.poststim   = 3;     % and some time after [s]
cfg.trialfun            = 'trialdef_encoding';
enc_cfg                = ft_definetrial(cfg);
enc_data               = ft_preprocessing(enc_cfg);


%% 6: Check data!
pause(0.01);
if checkdata_opt == 1
cfg = [];
ft_databrowser(cfg, enc_data)
end

%% 7: Save ping and no-ping data
save_path = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\eeg_data\';
save([save_path,sind,'_format'],'enc_data');

end
