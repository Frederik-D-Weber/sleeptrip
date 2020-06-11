% Copyright (C) 2020 Frederik D. Weber
%
% This file contains part of SleepTrip, see http://www.sleeptrip.org
% for the documentation and details and download at 
% https://github.com/Frederik-D-Weber/sleeptrip
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


%% pipeline for a sleeptrip analysis of 
% Forcato et al. 2020, Communications Biology

%%% PLEASE MAKE SURE:
%%% 1. YOU HAVE Matlab2013b or later
%%% 2. YOU HAVE ALL REQUIRED TOOLBOXES (typically the case)
%%% 3. YOUR SLEEPTRIP FOLDER DIRECTLY HAS THE FILES IN THEM, 
%%%    e.g. while extracting a zip from github maybe you need to refer to 
%%%    D:/sleeptrip-master/sleeptrip-master instead of D:/sleeptrip-master
%%% 4. NO FieldTrip NOR related toolboxes like EEGLAB are loaded or in the 
%%%    PATH of your MATLAB
%%% 5. You can downloaded example datasets here:
%%%    https://drive.google.com/open?id=160i_FG3zSZrBf9mr2zxZQX-u30OnLDMf
%%% 6. ...and the files from the example datasets are in the 
%%%    same folder as this script or somewhere in the matlab path 
   


%% initialization

%%% this asures the figures are later in vector format and editable 
set(0, 'DefaultFigureRenderer', 'painters');

%%% add sleeptrip to the path 
% e.g. addpath('D:/sleeptrip')
pathToSleepTrip = uigetdir('','choose spisop path, e.g. D:\sleeptrip');
addpath(pathToSleepTrip)

%%% disable some toolboxes if necessary and check if 
%%% 'signal_toolbox', 'signal_blocks' are available 
%%% because they are helpful to have.
% toggleToolbox('names')
% toggleToolbox('dsp','off')
% toggleToolbox('all','query')
% license('inuse')

%%% load all the defaults of SleepTrip and FieldTrip
st_defaults

[FileName,PathName,FilterIndex] = uigetfile({'*.edf';'*.zip';'*.*'},'Select zmax file (either one of the .edf or .zip');

[dummypathstr,FileNameWithoutExt,FileNameExt] = fileparts(FileName);

cfg = [];
cfg.dataset = [PathName FileName];
data = st_preprocessing(cfg);

% cfg = [];
% cfg.latency = [0 100];
% data = ft_selectdata(cfg, data);

%also read in the header for more info
hdr = ft_read_header(cfg.dataset);

recording_start_date = datetime(hdr.orig.T0(1),hdr.orig.T0(2),hdr.orig.T0(3),hdr.orig.T0(4),hdr.orig.T0(5),hdr.orig.T0(6),0);

%convert data to a table
datatable = cat(2,...
    table(datestr((recording_start_date + second(data.time{1}))','yyyy-mm-dd_HH:MM:SS.FFF'),'VariableNames',{'datetime'}),...
    array2table((1:numel(data.time{1}))','VariableNames',{'sample'}),...
    array2table(data.time{1}','VariableNames',{'seconds_since_recording_start'}),...
    array2table(data.trial{1}','VariableNames',genvarname(data.label)));

% write the file to csv
filepath_csv = [PathName FileNameWithoutExt '.csv'];

% [file,path] = uiputfile('*.csv','CSV file',filepath_csv);
% filepath_csv = [path file];

written_filepath = st_write_table(datatable,filepath_csv,',');

% messeag that this file is done
msgbox(['Finished! CSV file was written to ' written_filepath]);

