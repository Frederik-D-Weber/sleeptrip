% append either .edf or brainvision .eeg files, they can be compressed in a zip file 
% all files in a folder (including the subfolders) will be parsed for
% appending them. The order of merging is by the ascending order of
% filenames. Only works with SleepTrip.

% Copyright (C) 2019-, Frederik D. Weber
%
% This file is part of SleepTrip, see http://www.sleeptrip.org
% for the documentation and details.
%
%    SleepTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    SleepTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    SleepTrip is a branch of FieldTrip, see http://www.fieldtriptoolbox.org
%    and adds funtionality to analyse sleep and polysomnographic data. 
%    SleepTrip is under the same license conditions as FieldTrip.
%
%    You should have received a copy of the GNU General Public License
%    along with SleepTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$
clear

FileNamePostfixString = '_all_appended';
%OutputDataformat = 'edf_0.1uV_Ycuttoff';

%mainSpiSOPPath = 'D:\spisop_toolbox_beta2.3';
%mainSpiSOPPath = uigetdir('','choose spisop path, e.g. D:\spisop_toolbox_beta2.3');

%fieldtripPath = [mainSpiSOPPath '\fieldtrip_fw2016'];
fieldtripPath = uigetdir(pwd,'choose sleeptrip path, e.g. D:\sleeptrip');
external_codePath = [fieldtripPath filesep 'external' filesep 'enhanced_rdir'];

%folderPathWithFiles = 'B:\merge\merge';
folderPathWithFiles = uigetdir(fieldtripPath,'choose the input folder, e.g. C:\...\input');
%pathOutputFolder = 'B:\merge';
pathOutputFolder = uigetdir(fieldtripPath,'choose the output folder, e.g. C:\...\output');

if ~isdir(pathOutputFolder)
    mkdir(pathOutputFolder);
end

%addpath(mainSpiSOPPath)
addpath(fieldtripPath)
addpath(external_codePath)

st_defaults();

cd(folderPathWithFiles)

list_input = {'*.zip','*.edf','*.eeg'};
[selection,ok] = listdlg('PromptString','Input Data Format','SelectionMode','single','ListString',list_input);
if ~ok
    exit('next time select an input data format to run further')
end

data_format_ending = list_input{selection};

list_output = {'edf_0.1uV_Ycuttoff','brainvision_eeg_int16','brainvision_eeg_int32' 'brainvision_eeg_float32' 'edf_autoscale', 'edf_0.01uV_Ycuttoff'};
[selection,ok] = listdlg('PromptString','Output Data Format','SelectionMode','single','ListString',list_output);
if ~ok
    exit('next time select an output data format to run further')
end
data_format_output = list_output{selection};

edfFiles = rdir([folderPathWithFiles filesep '\**\' '*' data_format_ending]);


w = waitbar(0,'reading and merging data...Please wait!');

edfFilePathAndNames = {};

for k=1:length(edfFiles)
    
    edfFilePathAndNames{k} = edfFiles(k).name;
    
end

edfFilePathAndNames = sort(edfFilePathAndNames);



prompt = cellfun(@(s) [ 'before ' s],edfFilePathAndNames,'UniformOutput',false);
title = 'Add(+) or remove(-) seconds before datasets';
dims = [1 35];
definput = cellstr(num2str(zeros(numel(edfFilePathAndNames),1)))';
seconds_add_datasets_before = inputdlg(prompt,title,dims,definput);

prompt = cellfun(@(s) [ 'after ' s],edfFilePathAndNames,'UniformOutput',false);
title = 'Add(+) or remove(-) seconds after datasets';
dims = [1 35];
definput = cellstr(num2str(zeros(numel(edfFilePathAndNames),1)))';
seconds_add_datasets_after = inputdlg(prompt,title,dims,definput);

data = [];

oriDataFilename = 'appended';
oriDataExt = 'edf';
waitbar(0.02);
for k=1:length(edfFilePathAndNames)
    
    [fpathstr,fname,fext] = fileparts(edfFilePathAndNames{k});
    
    oriDataFilename = [fname];
    oriDataExt = [fext];
    
    cfg = [];
    cfg.continuous = 'yes'; %overwrite the trial uncontinuous data structure
    cfg.dataset = edfFilePathAndNames{k};
    cfg.channel = 'all';
    fprintf('dataset %i %s: read data\n',k,edfFilePathAndNames{k});
    data_temp = st_preprocessing(cfg);
    
    seconds_to_add = str2num(seconds_add_datasets_before{k});
    if seconds_to_add ~= 0;
        cfg = [];
        if seconds_to_add > 0
            samples_to_add = round(seconds_to_add*data_temp.fsample);
            for iTrTr = 1:numel(data_temp.trial)
                cfg.padtype = 'zero';
                data_temp.trial{iTrTr} = ft_preproc_padding(data_temp.trial{iTrTr}, cfg.padtype, samples_to_add, 0);
            end
            data_temp.time{1} = (0:(size(data_temp.trial{1},2)-1))/data_temp.fsample;
            data_temp.sampleinfo = [1 numel(data_temp.time{1})];
            
            % samples_to_precede = round(seconds_to_precede*data_temp.fsample);
            % data_temp.trial{1} = [zeros(size(data_temp.trial{1},1),samples_to_precede) data_temp.trial{1} ];
        elseif (seconds_to_add < 0) && (abs(seconds_to_add)*data_temp.fsample+1 < size(data_temp.trial{1},2))
            
            samples_to_add = round(seconds_to_add*data_temp.fsample);
            
            cfg.begsample = -samples_to_add+1;
            cfg.endsample = numel(data_temp.time{1});
            
            data_temp = ft_redefinetrial(cfg,data_temp);
            cfg = [];
            cfg.offset = samples_to_add+1;
            data_temp = ft_redefinetrial(cfg,data_temp);
            data_temp.sampleinfo = [1 numel(data_temp.time{1})];
        else
            error (['tried to cut too much a the beginning of file number' num2str(k) ', the cutting is longer than this file']);
        end
    end
    
    
    seconds_to_add = str2num(seconds_add_datasets_after{k});
    if seconds_to_add ~= 0;
        cfg = [];
        if seconds_to_add > 0
            samples_to_add = round(seconds_to_add*data_temp.fsample);
            for iTrTr = 1:numel(data_temp.trial)
                cfg.padtype = 'zero';
                data_temp.trial{iTrTr} = ft_preproc_padding(data_temp.trial{iTrTr}, cfg.padtype, 0, samples_to_add);
            end
            data_temp.time{1} = (0:(size(data_temp.trial{1},2)-1))/data_temp.fsample;
            data_temp.sampleinfo = [1 numel(data_temp.time{1})];
            
            % samples_to_precede = round(seconds_to_precede*data_temp.fsample);
            % data_temp.trial{1} = [zeros(size(data_temp.trial{1},1),samples_to_precede) data_temp.trial{1} ];
        elseif (seconds_to_add < 0) && (abs(seconds_to_add)*data_temp.fsample+1 < size(data_temp.trial{1},2))
            
            samples_to_add = round(seconds_to_add*data_temp.fsample);
            
            cfg.begsample = 1;
            cfg.endsample = numel(data_temp.time{1})+samples_to_add;
            
            data_temp = ft_redefinetrial(cfg,data_temp);
            %cfg = [];
            %cfg.offset = samples_to_add+1;
            %data_temp = ft_redefinetrial(cfg,data_temp);
            data_temp.sampleinfo = [1 numel(data_temp.time{1})];
        else
            error (['tried to cut too much a the beginning of file number' num2str(k) ', the cutting is longer than this file']);
        end
    end
    
    
    %     if seconds_to_add ~= 0;
    %
    %         samples_to_add = round(seconds_to_add*data_temp.fsample);
    %         data_temp.trial{1} = [data_temp.trial{1} zeros(size(data_temp.trial{1},1),samples_to_add)];
    %     end
    if k==1
        data = data_temp;
    else
        data.trial{1} = [data.trial{1} data_temp.trial{1}];
        data.time{1} = 0:(1/data.fsample):(size(data.trial{1},2)/data.fsample);
        if length(data.time{1}) == (1+length(data.trial{1}))
            data.time{1} = data.time{1}(1:(end-1));
        end
        data.sampleinfo = [1 size(data.trial{1},2)];
    end
    data_temp = [];
    waitbar(k/length(edfFilePathAndNames),w,[num2str(k) ' data set appended...']);
end


prompt = {'resample at sampling rate (Hz)\nPlease do NOT save in sampling rates that divide 2^n\n(e.g. 64 ,128 , 256, 512, ...)'};
title = 'Update sampling rate?';
dims = [1 35];
definput = cellstr(num2str([data.fsample]))';
updated_samplerate = inputdlg(prompt,title,dims,definput);
updated_samplerate = str2num(updated_samplerate{1});

if updated_samplerate ~= data.fsample
    cfg = [];
    cfg.resamplefs = updated_samplerate;%frequency at which the data will be resampled (default = 256 Hz)
    cfg.detrend = 'no';
    data = ft_resampledata(cfg,data);
end

hdr = data.hdr;
hdr.nSamples = size(data.trial{1},2);
hdr.nSamplesPre = 0;
hdr.nTrials = 1;
hdr.Fs = data.fsample;
hdr.nChans = length(data.label);
hdr.nTrials = 1;

tempOutputDataformat = data_format_output;
switch data_format_output
    case 'brainvision_eeg_int16'
        tempOutputDataformat = 'brainvision_eeg';
        hdr.brainvision_outformat = 'int16';%float32 int16 int32;
    case 'brainvision_eeg_int32'
        tempOutputDataformat = 'brainvision_eeg';
        hdr.brainvision_outformat = 'int32';%float32 int16 int32;
    case 'brainvision_eeg_float32'
        tempOutputDataformat = 'brainvision_eeg';
        hdr.brainvision_outformat = 'float32';%float32 int16 int32;
    case 'edf_autoscale'
        tempOutputDataformat = 'edf';
        hdr.edf_doautoscale = true;
    case 'edf_0.1uV_Ycuttoff'
        tempOutputDataformat = 'edf';
        hdr.edf_doautoscale = false;
        hdr.edf_accuracy = 0.1;
        hdr.edf_docutoff = true;
    case 'edf_0.01uV_Ycuttoff'
        tempOutputDataformat = 'edf';
        hdr.edf_doautoscale = false;
        hdr.edf_accuracy = 0.01;
        hdr.edf_docutoff = true;
    case 'edf_1uV_Ycuttoff'
        tempOutputDataformat = 'edf';
        hdr.edf_doautoscale = false;
        hdr.edf_accuracy = 1;
        hdr.edf_docutoff = true;
end
if strcmp(oriDataExt, '.zip')
    data_file_name = [pathOutputFolder filesep oriDataFilename FileNamePostfixString oriDataExt];
elseif strcmp(tempOutputDataformat, 'edf')
    data_file_name = [pathOutputFolder filesep oriDataFilename FileNamePostfixString '.edf'];
else
    data_file_name = [pathOutputFolder filesep oriDataFilename FileNamePostfixString];
end
ft_write_data([data_file_name], data.trial{:},'dataformat',tempOutputDataformat,'header',hdr);

close(w)

