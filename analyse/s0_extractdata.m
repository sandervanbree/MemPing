% Simple behavioural data extraction, calculating recognition
% performance per block, and reaction time in the encoding phase for
% correct and incorrectly recognized trials.
% Note: Behavioral data are unavailable for participant 1 to 4
% SvB

clear all; close all; clc;

% Params
slist = [5:15 17:33]; % participants
nblocks = 8; % number of blocks

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = behav_path;

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
    i_set_ret = findCol("curr_set_ret"); % retrieved stim indices
    i_soa_ret = findCol("curr_soa_ret"); % retrieved ping soas
    i_phs3_ans = findCol("key_resp_2_keys"); % retrieval: forgotten or remembered?
    i_phs4_ans = findCol("key_resp_5_keys"); % recognition: answer left or right?
    i_phs4_corr = findCol("corr_ans"); % recognition: truth left or right?
    i_enc_rt = findCol("key1_rt"); % reaction time during encoding
    i_enc_avail = findCol("key1_keys"); % just to index encoding trials
    i_enc_stim = find(string(tabl_enc.Properties.VariableNames) == "curr_set_enc"); % stimulus used during encoding phase (separate csv)
    i_rettest_stim = find(string(tabl_enc.Properties.VariableNames) == "rettest_trls"); % stimulus used during encoding phase (separate csv)
     
    % Transform into workable structure
    getCol = @(i) table2cell(tabl(:,i));
    set_ret_temp = getCol(i_set_ret);
    soa_ret_temp = getCol(i_soa_ret);
    phs3_ans_temp = getCol(i_phs3_ans);
    phs4_ans_temp = getCol(i_phs4_ans);
    phs4_corr_temp = getCol(i_phs4_corr);
    i_enc_rt_temp = getCol(i_enc_rt);
    i_enc_avail_temp = getCol(i_enc_avail);
    i_enc_stim_temp = table2cell(tabl_enc(:,i_enc_stim));
    i_rettest_stim_temp = table2cell(tabl_enc(:,i_rettest_stim));
    
    % SET_RET and SOA_RET
    set_ret_temp2 = set_ret_temp(~cellfun('isempty',set_ret_temp));
    soa_ret_temp2 = soa_ret_temp(~cellfun('isempty',soa_ret_temp));    
    set_enc_temp2 = i_enc_stim_temp(~cellfun('isempty',i_enc_stim_temp));
    enc_rt_temp2 = i_enc_rt_temp(~cellfun('isempty',i_enc_avail_temp)); % odd one out note: filtering enc_rt using rt_avail
    set_rettest_temp2 = unique(i_rettest_stim_temp(~cellfun('isempty',i_rettest_stim_temp)),'stable');

    for i = 1:nblocks % Patch things up
        t1 = soa_ret_temp2{i}(2:end-1);
        t2 = set_ret_temp2{i}(2:end-1);
        t3 = set_enc_temp2{i}(2:end-1);
        t4 = set_rettest_temp2{i}(2:end-1);
            
        curr_soa_ret = strsplit(t1);
        curr_set_ret = strsplit(t2);
        curr_set_enc = strsplit(t3);
        curr_set_rettest = strsplit(t4);
        
        curr_soa_ret(strcmp('',curr_soa_ret)) = [];
        curr_set_ret(strcmp('',curr_set_ret)) = [];
        curr_set_enc(strcmp('',curr_set_enc)) = [];
        curr_set_rettest(strcmp('',curr_set_rettest)) = [];
        
        curr_soa_ret = str2double(curr_soa_ret);
        curr_set_ret = str2double(curr_set_ret);
        curr_set_enc = str2double(curr_set_enc);
        curr_set_rettest = str2double(curr_set_rettest);
        
        soa_ret(i,:) = curr_soa_ret;
        set_ret(i,:) = curr_set_ret;
        set_enc(i,:) = curr_set_enc;
        set_rettest(i,:) = curr_set_rettest;
    end
    
    if pp > 4
        % PHS3_ANS
        phs3_ans = filt_ans(phs3_ans_temp);
        
        % PHS4_ANS
        phs4_ans = filt_ans(phs4_ans_temp);
        phs4_ans = (phs4_ans*-1)+3; % flip coding scheme
        
        % PHS4_CORR
        phs4_corr = filt_ans(phs4_corr_temp);

    elseif pp < 5
        phs3_ans = phs3_ans_temp;
        phs4_ans = phs4_ans_temp;
        phs4_corr = [];
    end
    
    % DICHOTOMIZE SOA_RET
    bin_edg = linspace(0.5,1.5,4);
    for i = 1:size(set_ret,1)
        for j = 1:size(set_ret,2)
            soa_temp = soa_ret(i,j);
            if soa_temp < bin_edg(2)
                soa_ret_dich(i,j) = 0;
            elseif soa_temp > bin_edg(2) && soa_temp < bin_edg(3)
                soa_ret_dich(i,j) = 1;
            elseif soa_temp > bin_edg(3) && soa_temp < bin_edg(4)
                soa_ret_dich(i,j) = 2;
            elseif soa_temp > bin_edg(4)
                soa_ret_dich(i,j) = 3;
            end
        end
    end
    
    if pp > 4
        % Score it
        for i = 1:nblocks %score across blocks
            ind = (i*10)-9;
            score_xblocks(i) = sum(phs4_ans(ind:ind+9) == phs4_corr(ind:ind+9))/10;
            score_xblocks(i) = sum(phs4_ans(ind:ind+9) == phs4_corr(ind:ind+9))/10;           
        end
        
        score = mean(score_xblocks); %average score
        
        % Display
        disp(score_xblocks); % every element is one block
        disp(score); % average across blocks
        
        
        enc_rt = cell2mat(enc_rt_temp2);
        
        corr_rt = [];
        false_rt = [];
        for i = 1:nblocks %reaction time as a function of correctness
            ind = (i*10)-9;
            vec = (phs4_ans(ind:ind+9) == phs4_corr(ind:ind+9));
            
            rettest_corr = set_rettest(i,vec);
            rettest_false = set_rettest(i,~vec);
                        
            [~,idx_corr] = intersect(rettest_corr,set_enc(i,:));
            [~,idx_false] = intersect(rettest_false,set_enc(i,:));
            
            curr_rt = enc_rt(ind:ind+9);
            
            corr_rt = [corr_rt ; curr_rt(idx_corr)];
            false_rt = [false_rt ; curr_rt(idx_false)];
        end

    else
        score_xblocks = [];
        score = [];
        corr_rt = [];
        false_rt = [];
    end
    
        % Display
        disp(corr_rt); 
        disp(false_rt); 
    
    % Save
    save([behav_path,'behav_res_',pp_string],'score_xblocks','score','set_ret','soa_ret',...
        'soa_ret_dich','set_rettest','phs3_ans','phs4_ans','phs4_corr','corr_rt','false_rt');
end
