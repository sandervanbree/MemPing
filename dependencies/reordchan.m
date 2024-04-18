% This function reorders FieldTrip channels according to a template

function [y] = reordchan(x,lay)

y = x;

l = lay.label;

[a,b] = setdiff(l,y.label);
l(b) = [];

for i = 1:numel(x.label)
    ord = find(strcmp(x.label{i},l));
    y.label{ord} = x.label{i};
    
    for i2 = 1:numel(x.trial)
        y.trial{i2}(ord,:) = x.trial{i2}(i,:);
    end
end
isq = isequal(l{numel(x.label)},y.label{numel(x.label)});
if isq == false
    error('reordering went wrong');
end

isq = isequal(l{numel(x.label)-5},y.label{numel(x.label)-5});
if isq == false
    error('reordering went wrong');
end
