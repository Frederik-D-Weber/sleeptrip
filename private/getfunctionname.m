function functionname = getfunctionname()
st = dbstack;
names = {st.name};
if numel(names)>1
    functionname = names{2};
else
    functionname = 'unknown';
end
end