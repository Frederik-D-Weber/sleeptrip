function v = roundsensible(v)

w = abs(v);
if numel(v)> 1
    w = w - min(w);
end
w = w(w~=0);
if ~isempty(w)
om = floor(log10(max(w)));
nowm0 = 2;
for iV = 1:numel(v)
    n = max(0,nowm0-om);
    v(iV) = round(v(iV), n);
end
end
end