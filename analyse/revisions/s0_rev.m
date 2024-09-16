% REVISIONS CODE
%
% Revisions analogue of s0_extractdata; i.e., load and format data to
% enable preprocessing.
%
% Some of the comments and variable names might not completely map onto
% what is being done, this is because the code was copy-pasted from the
% original scripts.
% SvB

clear all; close all; clc;

% Params
slist = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33]; % participants
nblocks = 8; % number of blocks

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\';
addpath(genpath(work_path));
eeg_path   = ['\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\data\eeg_data\'];
behav_path = ['\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\data\behav_data\'];
dep_path   = ['\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\dependencies\'];
save_path  =['\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\behav_data\'];

for pp = slist % loop across participants
    
    % Load params
    if pp <10
        pp_string = ['pp0',num2str(pp)];
    else
        pp_string = ['pp',num2str(pp)];
    end
    
    behav_file = ['behav_',pp_string,'.csv'];
    encstim_file = ['encstim_',pp_string,'.csv'];
    
    % Read in data
    tabl = readtable([behav_path,behav_file]);
    tabl_enc = readtable([behav_path,encstim_file]);
    
    % Find right columns
    findCol = @(name) find(string(tabl.Properties.VariableNames) == name);
    i_enc_rt = findCol("key1_rt"); % reaction time during encoding
    i_enc_avail = findCol("key1_keys"); % just to index encoding trials
    i_enc_stim = find(string(tabl_enc.Properties.VariableNames) == "curr_set_enc"); % stimulus used during encoding phase (separate csv)
   
    % Transform into workable structure
    getCol = @(i) table2cell(tabl(:,i));
    i_enc_rt_temp = getCol(i_enc_rt);
    i_enc_avail_temp = getCol(i_enc_avail);
    i_enc_stim_temp = table2cell(tabl_enc(:,i_enc_stim));
    
    % SET_RET and SOA_RET   
    set_enc_temp2 = i_enc_stim_temp(~cellfun('isempty',i_enc_stim_temp));
    enc_rt_temp2 = i_enc_rt_temp(~cellfun('isempty',i_enc_avail_temp)); % odd one out note: filtering enc_rt using rt_avail

    clear set_enc
    for i = 1:nblocks % Patch things up
        t3 = set_enc_temp2{i}(2:end-1);
        curr_set_enc = strsplit(t3);

        curr_set_enc(strcmp('',curr_set_enc)) = [];
 
        curr_set_enc = str2double(curr_set_enc);
        set_enc(i,:) = curr_set_enc;
    end
     

save(['\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\rev\data\behav_data\behav_res',pp_string],'set_enc');
end
