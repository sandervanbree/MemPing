% Sander van Bree -- Inspiration taken from Juan Linde-Domingo's function
% Input = 2D data (trial x chan) & win_opt = 0 (mean) or 1 (gauss-weighted)

function [output] = apply_win(data,win_opt)
n_trial = size(data,1);
n_chans = size(data,2);

if win_opt == 1 % Create a gaussian to perform weighted mean with
    x = -2:.01:2;
    norm = (normpdf(x,0,1))*2.5;
    norm = imresize(norm,[1 size(data,3)])';
end

for tr=1:n_trial
    for ch=1:n_chans
        if win_opt == 0
            output(tr,ch,:) = mean(squeeze(data(tr,ch,:)));
        elseif win_opt  == 1
            output(tr,ch,:) = wmean(squeeze(data(tr,ch,:)),norm,1);
        end
    end
end

return