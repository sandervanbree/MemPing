function [trl, event] = trialdef_noping(cfg)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% fetch info
samp = cfg.samp;
trig = cfg.trig;
soa_ret_dich_c  = cfg.trlinfo_ping(1,:); %ping SOA dichotomized
set_ret_c       = cfg.trlinfo_ping(2,:); %stim ID
set_ret_dich_c  = cfg.trlinfo_ping(3,:); %stim category mid/top
set_examp_c     = cfg.trlinfo_ping(4,:); %stim category bot
soa_ret_c       = cfg.trlinfo_ping(5,:); %pinging condition
soa_ret_shuff_c = cfg.trlinfo_noping(1,:); %no ping SOA
soa_ret_shuff_dich_c = cfg.trlinfo_noping(2,:); %no ping SOA dichotomized

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.prestim  * hdr.Fs);
posttrig =  round(cfg.trialdef.poststim * hdr.Fs);

% Loop
trl = [];
ind = 1; %Trigger ticker
ind2 = 1; %Shuffled SOA ticker
for i = 1:size(trig,2)-1
    
    if ind == 319
        break
    end
    
    if trig(i) == 21 % to not-be-pinged trial
        trlbegin = samp(i)+pretrig;
        trlbegin = trlbegin+round(soa_ret_shuff_c(ind2).*hdr.Fs);  % add random SOA grabbed from real SOAs
        trlend   = samp(i)+posttrig;
        trlend   = trlend+round(soa_ret_shuff_c(ind2).*hdr.Fs);  % add random SOA grabbed from real SOAs
        offset   = pretrig;
        newtrl   = [trlbegin trlend offset soa_ret_shuff_dich_c(ind2) set_ret_c(ind)...
            set_ret_dich_c(ind) set_examp_c(ind) round(soa_ret_shuff_c(ind2),2)];
        % NOTE: here the fourth and ninth column are NON-PING rather than PING soa
        
        if soa_ret_dich_c(ind) ~= 3
            error('there should be nothing except 3 here');
            % further context: this is a very important error check that
            % tracks the ping and no-ping conditions in the triggers with
            % those from the original Psychopy script. I.e., it should be exactly
            % those post-recording trials that do not have a ping as the
            % pre-recording script intended it (3 = no-ping).
        end
        
        % begin trial sample; end sample; offset to left; retrieved stim; soa of ping
        trl = [trl;newtrl];
        ind      = ind+1;
        ind2     = ind2+1;
    elseif trig(i) == 20 && trig(i+1) == 22 % If to-be-pinged trial followed by ping
        ind      = ind+1;
    elseif trig(i) == 20 && trig(i+1) == 20 % this happens a couple of times
        ind      = ind+1;
    elseif trig(i) == 20 && trig(i+1) == 21 % this happens a couple of times
        ind      = ind+1;
    end
end

figure;histogram(trl(:,5),200);xlabel('stim_idx');ylabel('frequency');
title('not pinged stims');
