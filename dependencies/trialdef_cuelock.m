function [trl, event] = trialdef_cuelock(cfg)

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
    if trig(i) == 20 || trig(i) == 21 % to-be-pinged or not-to-be-pinged trial
        trlbegin = samp(i)+pretrig;
        trlend   = samp(i)+posttrig;
        offset   = pretrig;
        newtrl   = [trlbegin trlend offset soa_ret_dich_c(ind) set_ret_c(ind)...
            set_ret_dich_c(ind) set_examp_c(ind) round(soa_ret_c(ind),2)];
        
        % begin trial sample; end sample; offset to left; dichotomized SOA
        % of ping; stimulus ID; top/mid level cat; bot level cat; soa
        % precise; stim;
        % soa of ping; 
        trl = [trl;newtrl];
        ind      = ind+1;
        
        if ind == 319
            break
        end
     end
end