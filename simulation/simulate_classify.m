% Stimulation demonstrating how Ntrials and Nclasses influence the variance
% of shuffled and empirical results, and how this is only consequential 
% toward p-values if there is a true effect.
%
% The "shuffle" classifier result is just a random permutation of the ground
% truth vector.
% The "empirical" classifier result is also a random permutation of ground
% truth vector, but with a portion of the ground truth manually inserted.

% Clean
clear all; clc

% Params
ntrials = 1000; % number of classified trials
nrep = 10^3; % number of classifier repetitions (Nshuffles & Nempirical)
nclass_top = 2; % Ntop level categories
nclass_bot = 16; % Nbot level categories
emp_power  = 0.6; % Fraction of trials accurately detected by empirical classifier for both levels
                  % 0 = no effect (empirical = shuffle)
                  % 0.95 = near-full effect (empirical = ground truth)

% Init variables
true_top   = randi(nclass_top,ntrials,1); % generate random ground truth seq (top)
true_bot   = randi(nclass_bot,ntrials,1); % same for bot (for simplicity, the same
                                          % structure holds across trials)
emp_insert = floor(ntrials*emp_power);    % how many trials of ground truth to insert into emp

for r = 1:nrep % loop across classifier repetitions
    
    % GENERATE NEW SHUFFLE
    shuff_top  = randi(nclass_top,ntrials,1); % generate random shuffle (top)
    shuff_bot  = randi(nclass_bot,ntrials,1); % same for bot
    
    % GENERATE NEW EMPIRICAL
    emp_top                 = randi(nclass_top,ntrials,1); % random pattern as with shuffle
    emp_bot                 = randi(nclass_bot,ntrials,1);
    emp_top(1:emp_insert)   = true_top(1:emp_insert); % But insert some of the true structure
    emp_bot(1:emp_insert)   = true_bot(1:emp_insert);
        
    % EVALUATE
    for k = 1:ntrials          % Derive accuracy by comparing for each trial whether
                               % estimated class matches true class
        top_shuff_tmp(k) = shuff_top(k)==true_top(k);
        bot_shuff_tmp(k) = shuff_bot(k)==true_bot(k);
        
        top_emp_tmp(k) = emp_top(k)==true_top(k);
        bot_emp_tmp(k) = emp_bot(k)==true_bot(k);
    end
    
    % To fetch accuracy, calculate proportion of accurately guessed trials
    top_shuff_score(r) = mean(top_shuff_tmp);
    bot_shuff_score(r) = mean(bot_shuff_tmp);
    
    top_emp_score(r) = mean(top_emp_tmp);
    bot_emp_score(r) = mean(bot_emp_tmp);
end

%% Plot results
% TOP EMPIRICAL
figure; hold on; subplot(221); 
raincloud_plot(top_emp_score,'color',[0.8 0.2 0.2]); title('top emp');
xlim([0 1]);
xlabel('accuracy'); set(gca, 'yTickLabel', [])
xline(1/nclass_top,'LineWidth',4); % chance level line

% BOT EMPIRICAL
subplot(222)
raincloud_plot(bot_emp_score,'color',[0.2 0.2 0.8]); title('bot emp');
xlim([0 1]);
xlabel('accuracy'); set(gca, 'yTickLabel', [])
xline(1/nclass_bot,'LineWidth',4); % chance level line

% TOP SHUFFLE
subplot(223)
raincloud_plot(top_shuff_score,[0.4 0.2 0.2]); title('top shuff');
xlim([0 1]);
xlabel('accuracy'); set(gca, 'yTickLabel', [])
xline(1/nclass_top,'LineWidth',4); % chance level line

% BOT SHUFFLE
subplot(224)
raincloud_plot(bot_shuff_score,[0.2 0.2 0.4]); title('bot shuff');
xlim([0 1]);
xlabel('accuracy'); set(gca, 'yTickLabel', [])
xline(1/nclass_bot,'LineWidth',4); % chance level line

% Note: the vertical line is chance level (= shuffled mean)

%% Findings
% Note how with increasing trials AND classes, the spread of shuffle and
% empirical reduces.
%
% Now, if there is a true effect (emp_power > 0), Ntrial & Nclass will increase
% the distance between empirical and shuffle. 
% --> thus also reducing p-vals

% If there is NO true effect (emp_power = 0), Ntrial & Nclass do not
% influence the distance between empirical and shuffle because they are
% both at chance, so narrowing of spread helps both equally
% --> thus no reduction of p-vals with Ntrial & Nclass

% Only if there is a true effect will higher Ntrials and Nclasses reduce
% p-values.

