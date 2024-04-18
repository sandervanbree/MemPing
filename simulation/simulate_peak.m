% This script serves to simulate the best way to detect a peak in time
% series with a single peak and noise. The best way to derive ground truth
% peaks will be used to detect empirical peaks

clear all; clc

% Load data
load classdata

% Simulation parameters
t = tvec_cl;  % Time values for generated signal
npp = size(ps_all,2); % Number of participants
snr = 1.15; % Signal-to-noise ratio for peak
autocorr_lvl = 0.4; % Autocorrelation strength parameter
ntrials = 1000; % Number of trials to simulate
pk_loc = 1; % Location of peak in seconds
peakwidth = 0.04; %width of peak

smoothbins = 10; % Smoothing window size if smoothing enabled
smoothmeth = 1; % 0 = no filt ; 1 = gaussfilt, 2 = sgolayfilt, 3 = medfilt1

% Preallocate
peakLocs = nan(npp,ntrials,8);
errdist = nan(npp,ntrials,8);

% Loop through participants
for pp = 1:npp
    
    disp(['preparing participant ', num2str(pp)]);
    
    % Get signal information for current participant
    dat = squeeze(p_all(3,pp,:)); % do for all SOAs together

    % Generate signal
    for tr = 1:ntrials
        % autocorrelation
        rx_rand = zeros(numel(dat),1);
        for i = 2:numel(dat)
            step = randi([0 1],1);
            if step == 0
                rx_rand(i) = (rx_rand(i-1) - randn);
            else
                rx_rand(i) = (rx_rand(i-1) + randn);
            end
        end
        
        noise = randn(size(t));
        peak = snr*exp(-((t-pk_loc).^2)/peakwidth);
        sig = noise + peak + (rx_rand.*autocorr_lvl)';
        
        % smooth
        if smoothmeth == 0
            sig = sig;
        elseif smoothmeth == 1
            sig = imgaussfilt(sig,smoothbins);
        elseif smoothmeth == 2
            sig = sgolayfilt(sig,3,11);
        elseif smoothmeth == 3
            sig = medfilt1(sig,smoothbins);
        end
        
        
        %% Peak detection approaches
        % Approach 1 - find peak
        filtsig = lowpass(sig,10,numel(sig));
        [a,b] = findpeaks(filtsig);
        [~,c] = max(a);
        locs1 = b(c);
        
        % Approach 2 - find Max
        [~, locs2] = max(sig);
        locs2 = locs2(1); % Convert to scalar
        
        % Approach 3 - Finding peaks in the cumulative sum
        csum = cumsum(sig);
        csum_diff = diff(csum);
        [a,b] = findpeaks(csum_diff);
        [~,c] = max(a);
        locs3 = b(c);
        
        % Approach 4 - Area under curve
        ctrapz = cumtrapz(sig);
        [a,b] = findpeaks(diff(ctrapz));
        [~,c] = max(a);
        locs4 = b(c);
        
        % Approach 5 - Area under curve (cumulative sum)
        ctrapz = cumtrapz(csum_diff);
        [a,b] = findpeaks(diff(ctrapz));
        [~,c] = max(a);
        locs5 = b(c);
        
        % Approach 6 - Wavelet transform
        [cfs,frs] = wavedec(sig,3,'db4');
        wv = waverec(cfs,frs,'db4'); % reconstruct detail coef
        [a,b] = findpeaks(wv);
        [~,c] = max(a);
        locs6 = b(c);
        
        % Approach 7 - Hilbert transform
        hil_sig = hilbert(sig);
        amp_env = abs(hil_sig);
        [a,b] = findpeaks(wv);
        [~,c] = max(a);
        locs7 = b(c);
        
        % Approach 8 - xCorr template matching
        % Create Gaussian peak
        t_temp = linspace(-0.75,0.75,75);
        template = exp(-(t_temp-0).^2 /peakwidth);
        [corr_sig,corr_bin] = xcorr(sig,template);
        [a,b] = findpeaks(abs(corr_sig));
        [~,c] = max(a);
        toindex = b(c);
        shiftbin = corr_bin(toindex);
        zeroind = find(t_temp==0);
        locs8 = zeroind+shiftbin;
        
        % Store results
        locsvec = nan(1,8);
        locstim = nan(1,8);
        
        for i = 1:8
            if ~isempty(eval(['locs',num2str(i)]))
                locsvec(i) = eval(['locs',num2str(i)]); % the answers we got
                
                locstim(i) = t(locsvec(i)); % what is each method's guess?
                errdist(pp,tr,:) = abs(pk_loc-locstim);
            end
        end
        
        peakLocs(pp,tr,:) = locsvec;                     % store them
        
    end
    
    % Save in dist
    ns_dist(pp,:) = noise;
    pk_dist(pp,:) = peak;
    rx_dist(pp,:) = rx_rand;
    sig_dist(pp,:) = sig;
end

%% Plot example data
% last participant
figure; subplot(121); hold on; title('example participant');
plot(t,noise,'LineWidth',1,'Color',[0.8 0.2 0.8]);
plot(t,rx_rand,'LineWidth',1,'Color',[0.2 0.2 0.2]);
plot(t,peak,'LineWidth',1,'Color',[0.8 0.8 0.2]);
plot(t,sig,'LineWidth',3,'Color',[0.2 0.8 0.8]);
legend({'noise','autocorr','peak','signal'});
xlabel('time');
ylabel('amp');

% Plot results
% average participants
subplot(122); hold on; title('average participant');
plot(t,nanmean(ns_dist),'LineWidth',1,'Color',[0.8 0.2 0.8]);
plot(t,nanmean(rx_dist),'LineWidth',1,'Color',[0.2 0.2 0.2]);
plot(t,nanmean(pk_dist),'LineWidth',1,'Color',[0.8 0.8 0.2]);
plot(t,nanmean(sig_dist),'LineWidth',3,'Color',[0.2 0.8 0.8]);
legend({'noise','autocorr','peak','signal'});
xlabel('time');
ylabel('amp');

%% Score methods
errscore = squeeze(nanmean(squeeze(nanmean(errdist,1)),1));

methods = {'Peakfind','Max','Csum','Deriv', 'Csum Deriv','Wavelet','Hilbert','Xcorr'};

[~, minIdx] = sort(errscore);

% BAR PLOT

figure; bar(errscore,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1);

hold on
plot(minIdx(1), errscore(minIdx(1)), 'r.','MarkerSize',60)
plot(minIdx(2), errscore(minIdx(2)), 'y.','MarkerSize',60)

text(minIdx(1),errscore(minIdx(1)) + 0.05,num2str(errscore(minIdx(1)),'%0.3f'),'FontSize',10,'HorizontalAlignment','center')
text(minIdx(2),errscore(minIdx(2)) + 0.05,num2str(errscore(minIdx(2)),'%0.3f'),'FontSize',10,'HorizontalAlignment','center')

hold off
title(['method results: snr',num2str(snr),' | autocorr', num2str(autocorr_lvl),' | smoothmeth',num2str(smoothmeth)]);
set(gca,'xticklabel',methods)
xlabel('Method','FontSize',12,'Interpreter','latex')
ylabel('Error (ms)','FontSize',12,'Interpreter','latex')


% VIOLIN PLOT
err2d =  reshape(errdist, [], size(errdist,3));

kd = ksdensity(err2d(:,5));

figure;plot(kd,'linewidth',4,'color',[0.2 0.2 0.2]);
xlabel('Error (ms)','FontSize',12,'Interpreter','latex')
ylabel('Occurence density','FontSize',12,'Interpreter','latex')


figure;h= histogram(err2d(:,5),'Normalization','probability'); morebins(h); morebins(h); morebins(h);
% Label axes
xlabel('Error [s]'); 
ylabel('Percentage');

% Make plot look nicer
set(gca,'linewidth',1.5,'fontsize',14,'fontname','times');
box off;

% 
% 
% 
% 
% figure;plot(scatter(1,err2d(:,1)));
% 
% figure; violinplot(err2d,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1);
% 
% 
% violin(Y,'xlabel',{'a','b','c','d'},'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','b',...