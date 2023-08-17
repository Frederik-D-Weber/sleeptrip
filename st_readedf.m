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
    return
end
srate=srates(1);
numSamples=edf_info.NumDataRecords*edf_info.NumSamples(1);

time_vect = linspace(0, numSamples-1, numSamples).'/srate;

%---get data
fprintf('reading data...\n')
edf_data=edfread(edf_file);
edf_data=cell2mat(edf_data{:,chanInds})';

% edf_info=edfinfo(edf_file);
% edf_data=edfread(edf_file);
%
% RecTime = seconds(edf_data.('Record Time')); %edf-specific
%
% if iscell(edf_data{:,1})
%
%     signal = cat(1,edf_data{:,1}{:});       % Concatenate The data
%     srate = numel(edf_data{1,1}{1})/mean(diff(RecTime));                                  % Sampling Intervals (Samples/Second)
%
% else
%     signal=edf_data{:,1};
%     srate = numel(edf_data{1,1})/mean(diff(RecTime));                                  % Sampling Intervals (Samples/Second)
%
% end
%
% time_vect = linspace(0, numel(signal)-1, numel(signal)).'/srate;
%
% %%quick version
% srates=edf_info.NumSamples./seconds(edf_info.DataRecordDuration)
%
% useChans={'EEG L Cleaned','EEG R Cleaned'};
% edf_selected=edfread(edf_file,'SelectedSignals',useChans)

%%


%data_edf.label={convertStringsToChars(edf_info.SignalLabels)}; %ZMax has non-unique signallabels...
data_edf.label=useChans(:);
data_edf.time={time_vect'};
data_edf.trial={edf_data};
data_edf.fsample=srate;
data_edf.sampleinfo=[1 numSamples];

%construct header from available data
data_edf.hdr=ft_fetch_header(data_edf);
data_edf.hdr.orig=edf_info;

