% Let's start by plotting behavoural information --
% specifically with regard to within and cross-block performance, and
% reaction time for encoding as a function of whether trials were correctly
% recognized in phase 4 or not.
% SvB

clear all; close all; clc;

%% 0: Before starting
ft_defaults

% Parameters
pp_range = [5:15 17:33]; % loop over every participant from 5 until the last, except 16.

work_path  = '\\analyse4.psy.gla.ac.uk\project0318\Sander\memflash\memping\';
addpath(genpath(work_path));
eeg_path   = [work_path,'data\eeg_data\'];
behav_path = [work_path,'data\behav_data\'];
dep_path   = [work_path,'dependencies\'];
save_path  = [work_path,'analyse\results\'];

% Init
corr_rt_dist = [];
false_rt_dist = [];

ind = 1;
for pp = pp_range
    
    % Set stuff up
    if pp < 10
        sind = ['pp0',num2str(pp)];
    else
        sind = ['pp',num2str(pp)];
    end
    
    % 1: Load data
    behav_data = [behav_path,'behav_res_',sind,'.mat'];
    load(behav_data);
    
    % 2: Calculate current participant's score per block
    for block = 1:8
        curr_score = phs4_ans((block-1)*10+1:block*10) == phs4_corr((block-1)*10+1:block*10);
        score_overall(ind,block) = mean(curr_score);
    end
    
    % 3: Reaction time
    corr_rt_dist = [corr_rt_dist ; corr_rt];
    false_rt_dist = [false_rt_dist ; false_rt];
    
    % 4: Subjective score
     subj_rem = sum(find(phs3_ans==2));
     subj_forg = sum(find(phs3_ans==1));
     subj_rat(ind) = subj_rem./(subj_forg+subj_rem);
    
    % 5: Ping vs No-ping score

     ind_ping = find(soa_ret(:) < 99);
     ind_noping = find(soa_ret(:) == 99);
     
     set_ret_flat = set_ret(:);
     stim_id_ping = set_ret_flat(ind_ping);
     stim_id_noping = set_ret_flat(ind_noping);
     
     [~,loc_ping] = intersect(set_rettest,stim_id_ping);
     [~,loc_noping] = intersect(set_rettest,stim_id_noping);
     
     score_ping(ind) = mean(phs4_ans(loc_ping) == phs4_corr(loc_ping));
     score_noping(ind) = mean(phs4_ans(loc_noping) == phs4_corr(loc_noping));
     
    ind = ind+1;
end

%% Visualize 1: Average across participants for each block
avg_per_block = nanmean(score_overall,1);
sd_blocks = std(avg_per_block);

% Set figure properties
f1 = figure('Position',[100 100 600 400]);
subplot(211);hold on;
 bar(1:length(avg_per_block),avg_per_block,'BarWidth', 0.8,'FaceColor',[0.3 0.5 0.5]);

% Tighter y-axis 
ylim([0.75 1])
yticks([0.75 0.8 0.85 0.9 0.95 1]) 

% Add x-axis labels
xlabel('Block number','FontSize',14) 

% Add axis labels  
ylabel('Recognition accuracy [%]','FontSize',14)
title('Average per block','FontSize',14)

% Font
set(gca,'FontSize',12) 
set(gca,'FontName','Arial'); hold on;

%% Visualize 2: Average across blocks for each participant
avg_per_participant = mean(score_overall,2);
sd_pp = std(avg_per_participant);

% Set figure properties
subplot(212)
bar(1:numel(avg_per_participant),avg_per_participant,'BarWidth', 0.8,'FaceColor',[0.4 0.4 0.4]);

% Tighter y-axis 
ylim([0.75 1])
yticks([0.75 0.8 0.85 0.9 0.95 1]) 

% Add x-axis labels
xlabel('Participant','FontSize',14) 

% Add axis labels  
ylabel('Recognition accuracy [%]','FontSize',14)
title('Average per participant','FontSize',14)

% Font
set(gca,'FontSize',12) 
set(gca,'FontName','Arial'); hold off;

%% Visualize 3: RT distribution
avg_rt_corr = nanmean(corr_rt_dist); % mean
std_rt_corr = nanstd(corr_rt_dist); % st deviation
sem_rt_corr = std_rt_corr./sqrt(numel(corr_rt_dist)); %standard error

avg_rt_false = nanmean(false_rt_dist);
std_rt_false = nanstd(false_rt_dist);
sem_rt_false = std_rt_false./sqrt(numel(false_rt_dist));

% Transform to millisecs
avg_rt_corr = avg_rt_corr.*1000; % mean
std_rt_corr = std_rt_corr.*1000; % st deviation
sem_rt_corr = sem_rt_corr.*1000; %standard error

avg_rt_false = avg_rt_false.*1000;
std_rt_false = std_rt_false.*1000;
sem_rt_false = sem_rt_false.*1000;

% Plot
xpos1 = 1.1;
xpos2 = 1.9;

f3 = figure('Position',[100 100 600 400]);
b1 = bar(xpos1,avg_rt_corr,'FaceColor',[0.4 0.2 0.4],'BarWidth',0.5);
hold on
b2 = bar(xpos2,avg_rt_false,'FaceColor',[0.4 0.4 0.4],'BarWidth',0.5);
ylim([0 3400]);
xlim([0.5 2.5]);
xticks([xpos1 xpos2])
xticklabels({'Recognized','Forgotten'})
ylabel('Reaction time [ms]')

errorbar(xpos1,avg_rt_corr,sem_rt_corr,'k','linestyle','none');  
errorbar(xpos2,avg_rt_false,sem_rt_false,'k','linestyle','none');

% Add x-axis labels  
xticks([xpos1 xpos2])
xticklabels({'Recognized','Incorrect'})

% Add n values
text(xpos1+0.13,avg_rt_corr+100,sprintf('n=%d',length(corr_rt_dist)), 'Horiz','center')
text(xpos2+0.13,avg_rt_false+100,sprintf('n=%d',length(false_rt_dist)), 'Horiz','center')

% Font
set(gca,'FontSize',12) 
set(legend,'FontSize',12)
set(gca,'FontName','Arial') 
set(legend,'FontName','Arial') 

% Add axis labels
ylabel('Reaction time [ms]')
title('Reaction time during encoding')

% Error bars
errorbar(xpos1,avg_rt_corr,sem_rt_corr,'k','linestyle','none','linewidth',2);
errorbar(xpos2,avg_rt_false,sem_rt_false,'k','linestyle','none','linewidth',2);

% Legend 
legend([b1 b2],{'Recognized','Incorrect'})

% Font
set(gca,'FontSize',12) 
set(legend,'FontSize',12)
set(gca,'FontName','Arial') 
set(legend,'FontName','Arial') 

%% Visualize 3: Raincloud
f4 = figure('Position',[100 100 600 400]);
subplot(211);
b1 = raincloud_plot(corr_rt_dist.*1000,'box_on',1,'color',[0.5 0.3 0.5],'line_width',1.2,'band_width',300);
xlim([-1000 6000]);
legend([b1{1}], {'Recognized'});
set(gca,'FontSize',12) 
set(gca,'FontName','Arial') 

subplot(212);
b2 = raincloud_plot(false_rt_dist.*1000,'box_on',1,'color',[0.3 0.3 0.5],'line_width',1.2,'band_width',400);
xlim([-1000 6000]);
legend([b2{1}], {'Forgotten'});
set(gca,'FontSize',12) 
set(gca,'FontName','Arial') 

%% Descriptive stats subjective
disp(mean(subj_rat));
disp(var(subj_rat));

%% Ping vs No ping
% t-test
f5 = figure('Position', [100, 100, 600, 400]);

% calculating means and standard errors
mean_ping = mean(score_ping);
mean_noping = mean(score_noping);
std_err_ping = std(score_ping) / sqrt(length(score_ping));
std_err_noping = std(score_noping) / sqrt(length(score_noping));

% bar plot with error bars
bar(1, mean_ping, 'FaceColor', [0.2 0.2 0.2]);
hold on;
errorbar(1, mean_ping, std_err_ping, 'k');
bar(2, mean_noping, 'FaceColor', [0.8 0.8 0.8]);
errorbar(2, mean_noping, std_err_noping, 'k');
% Adding individual points
x_jitter = 0.1; % jitter to avoid overlap
scatter(ones(size(score_ping)) * 1 - x_jitter + rand(size(score_ping)) * 2 * x_jitter, score_ping, 'ko');
scatter(ones(size(score_noping)) * 2 - x_jitter + rand(size(score_noping)) * 2 * x_jitter, score_noping, 'ko');

xlim([0.5 2.5]);
ylim([0.8 1]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'Ping', 'No Ping'}, 'FontSize', 12, 'FontName', 'Arial');
ylabel('recognition %');

% t-test
[h, p, ci, stats] = ttest(score_ping, score_noping);

% display results
fprintf('T-test p-value: %f\n', p);
fprintf('Ping Mean ± SE: %f ± %f\n', mean_ping, std_err_ping);
fprintf('No-Ping Mean ± SE: %f ± %f\n', mean_noping, std_err_noping);
