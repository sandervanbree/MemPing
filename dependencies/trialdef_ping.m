function [trl, event] = trialdef_ping(cfg)

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
ind = 1;
for i = 1:size(trig,2)-1
    if trig(i) == 20 && trig(i+1) == 22 % If to-be-pinged trial followed by ping

        trlbegin = samp(i+1)+pretrig;
        trlend   = samp(i+1)+posttrig;
        offset   = pretrig;
        newtrl   = [trlbegin trlend offset soa_ret_dich_c(ind) set_ret_c(ind)...
            set_ret_dich_c(ind) set_examp_c(ind) round(soa_ret_c(ind),2)];
        
        if newtrl(4) == 3
            error('there should be no 3 here');
        end

        % begin trial sample; end sample; offset to left; retrieved stim; soa of ping
        trl = [trl;newtrl];
        ind      = ind+1;        
    elseif trig(i) == 21 % to be not pinged
        ind      = ind+1;
    elseif trig(i) == 20 && trig(i+1) == 20 % this happens a couple of times
        ind      = ind+1;
    elseif trig(i) == 20 && trig(i+1) == 21 % this happens a couple of times
        ind      = ind+1;
    end
    
    if ind == 319
        break
    end
end

figure;histogram(trl(:,5),200);xlabel('stim_idx');ylabel('frequency');
title('pinged stims');

if any(trl(:,4)==3)
    error('oops; there should be 3s here');
end