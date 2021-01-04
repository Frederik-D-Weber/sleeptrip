function [cfg filepaths] = st_write_continuous_data(cfg, data)

% ST_WRITE_CONTINUOUS_DATA writes data out in various standard sleep eeg formats
%
% Use as
%   [cfg filepaths] = st_write_continuous_data(cfg, data)
%
% Required configuration parameters are:
%   cfg.filename  = string with the filename (or full path) to write to, if
%                   not existent folder will be created
%
% Optional configuration parameters are:
%   cfg.format  = string, either 
%                 'edf' for edf with fixed 0.1 uV/units accuracy and fixed
%                     cuttoffs to fit 16bit range
%                 'brainvision_int16' (= default)
%                 'brainvision_int32'
%                 'brainvision_float32'... for brainvision .eeg and .vhdr files
%                     using a 16bit format and int16 instead of the typical float32 
%                     saves half the disc space. Note that 16bit typically
%                     are typically enough for high-pass filtered data.
%   cfg.compress = string,if to compress the files, either 'no' or 'zip' (default = 'zip')
%                  zip the files take up only ~70% of space, so this is
%                  default.
%   cfg.posmarker = string, either 'yes' or 'no' (default = 'yes') 
%                   if a positivity marker of 402 samples (201 samples rising from
%                   0 to 100 units then 201 samples falling from 100 to 0
%                   again) a the beginning of the recording.
%                   this is to check that the data is later correctly
%                   interpreted
%
% See also FT_PREPROCESSING, FT_APPLY_MONTAGE

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

ttic = tic;
mtic = memtic;
functionname = getfunctionname();
fprintf([functionname ' function started\n']);

% set defaults
%cfg.channel  = ft_getopt(cfg, 'channel', 'all', 1);
cfg.format  = ft_getopt(cfg, 'format', 'brainvision_int16');
cfg.compress  = ft_getopt(cfg, 'compress', 'zip');
cfg.posmarker = ft_getopt(cfg, 'posmarker', 'yes');

if ~isfield(cfg, 'filename')
 ft_error('please specify cfg.filename')
end     

fprintf([functionname ' function initialized\n']);


if strcmp(cfg.posmarker,'yes')
    data.trial{1}(:,1:402) = repmat([((0:(1/200):1)*100) ((1:-(1/200):0)*100)],size(data.trial{1},1),1);
end

%first create a header for the exported file
hdr = [];
hdr.nSamples = size(data.trial{1},2);
hdr.nSamplesPre = 0;
hdr.nTrials = 1;
hdr.Fs = data.fsample;
hdr.nChans = length(data.label);
hdr.label = data.label;
hdr.nTrials = 1;

% chose a file name (without extension)
[pathstr, name, ext] = fileparts(cfg.filename);

file_name = [pathstr name];


filepaths = '';

switch cfg.format
    case 'edf'
        % export as edf
        data_format_output = 'edf';
        hdr.edf_doautoscale = false; % if the data should be autoscaled, then hdr.edf_accuracy and hdr.edf_docutoff are not relevant
        hdr.edf_accuracy = 0.1; % in uV, e.g. 1 0.5, 0.1, 0.01 ...
        hdr.edf_docutoff = true; % cut the data at the limits
        file_extension = '.edf';
        % this is the final full file path
        data_export_filepath = [file_name file_extension];
        ft_write_data(data_export_filepath, data.trial{:},'dataformat',data_format_output,'header',hdr);
        
        filepaths = data_export_filepath;
        
    case {'brainvision_int16', 'brainvision_int32', 'brainvision_float32'}
        
        % export as brainvision int16 format
        data_format_output = 'brainvision_eeg';
        switch cfg.format
            case 'brainvision_int16'
                hdr.brainvision_outformat = 'int16';%float32 int16 int32;
            case 'brainvision_int32'
                hdr.brainvision_outformat = 'int32';%float32 int16 int32;
            case 'brainvision_float32'
                hdr.brainvision_outformat = 'float32';%float32 int16 int32;
        end
        
        file_extension = ''; %brainvision exports do not have an extension
        % this is the final full file path
        data_export_filepath = [file_name file_extension];
        ft_write_data(data_export_filepath, data.trial{:},'dataformat',data_format_output,'header',hdr);
        
        filepaths = {[data_export_filepath '.eeg'], [data_export_filepath '.vhdr']};
        
    otherwise
        ft_error('cfg.format = %s is unknown', cfg.format)
end

if strcmp(cfg.compress, 'zip')
    %zip the files take up only ~70% of space
    fprintf('compressing the data in .zip and deleting uncompressed original.')
    zipfilepaths = filepaths;
    zipfilepath = [data_export_filepath '.zip'];
    zip(zipfilepath,zipfilepaths);
    if iscell(filepaths)
        delete(filepaths{:})
    else
        delete(filepaths)
    end
    filepaths = zipfilepath;
end


fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end
