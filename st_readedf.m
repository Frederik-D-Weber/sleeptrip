function data_edf=st_readedf(edf_file,varargin)

%construct FT data struct
data_edf=struct;

%----get header
fprintf('reading header...\n')
edf_info=edfinfo(edf_file);
edf_chanLabels=convertStringsToChars(edf_info.SignalLabels);

if nargin>1
    useChans=varargin{1};
    [~,chanInds]=ismember(useChans,edf_chanLabels);

else
    useChans=edf_chanLabels;
    chanInds=1:edf_info.NumSignals;
end

%srates of requested channels
srates=edf_info.NumSamples(chanInds)./seconds(edf_info.DataRecordDuration);

if ~all(srates==srates(1))
    fprintf('unequal sample rates, aborting\n')
    return
end
srate=srates(1);
numSamples=edf_info.NumDataRecords*edf_info.NumSamples(1);

time_vect = linspace(0, numSamples-1, numSamples).'/srate;

%---get data
fprintf('reading data...\n')
edf_data=edfread(edf_file);
edf_data=cell2mat(edf_data{:,chanInds})';

if ~iscell(useChans)
    useChans={useChans};
end

%populate data struct
data_edf.label=useChans(:);
data_edf.time={time_vect'};
data_edf.trial={edf_data};
data_edf.fsample=srate;
data_edf.sampleinfo=[1 numSamples];

%construct header from available data
data_edf.hdr=ft_fetch_header(data_edf);
data_edf.hdr.orig=edf_info;

