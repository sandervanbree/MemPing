function y = filt_ans(x)
    for i = 1:numel(x)
        if contains(x(i),'left')
            tmp(i) = 1;
        elseif contains(x(i),'right')
            tmp(i) = 2;
        elseif contains(x(i),'None')
            tmp(i) = 3;
        else tmp(i) = 0;
        end
    end
    y = tmp(tmp>0);   
end