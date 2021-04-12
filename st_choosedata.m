function cfg = st_choosedata(cfg)

% ST_CHOOSEDATA choses datafile, headerfile, montagefile, datagrammerfile, scoringfile
%
% Use as
%   cfg = st_choosedata(cfg)
%
% Optional configuration parameters are:
%   cfg.datafile   = if to search for a data file either 'ask' or 'no' or a
%                    file path to the file. If the file exists it will not be asked for the
%                    file, if not it will ask for this file. 
%                    (default = 'ask')
%   cfg.headerfile = if to search for a header file either 'ask' or 'no' or a
%                    file path to the file. If the file exists it will not be asked for the
%                    file, if not it will ask for this file. 
%                    (default = 'no')
%   cfg.montagefile = if to search for a montage file either 'ask' or 'no' or a
%                    file path to the file. If the file exists it will not be asked for the
%                    file, if not it will ask for this file. 
%                    (default = 'no')
%   cfg.datagrammerfile = if to search for a datagrammer file either 'ask' or 'no'or a
%                    file path to the file. If the file exists it will not be asked for the
%                    file, if not it will ask for this file.
%                    (default = 'no')
%   cfg.scoringfile = if to search for a scoring file either 'ask' or 'no'or a
%                    file path to the file. If the file exists it will not be asked for the
%                    file, if not it will ask for this file.
%                    (default = 'no')
%
% See also ST_READ_SCORING, FT_PREPROCESSING, FT_APPLY_MONTAGE,
% ST_READ_MONTAGE, ST_APPLY_GRAMMER, ST_DATAGRAMMER, ST_SCORINGCONVERT

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

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
st_defaults
ft_preamble init
ft_preamble debug
% ft_preamble loadvar data
% ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% set defaults
cfg.datafile  = ft_getopt(cfg, 'datafile', 'ask', 1);
cfg.headerfile  = ft_getopt(cfg, 'headerfile', 'no', 1);
cfg.montagefile  = ft_getopt(cfg, 'montagefile', 'no', 1);
cfg.datagrammerfile  = ft_getopt(cfg, 'datagrammerfile', 'no', 1);
cfg.scoringfile  = ft_getopt(cfg, 'scoringfile', 'no', 1);

fprintf([functionname ' function initialized\n']);

skipNonData = false;
tmpfilepath = '';
[askdatafile hasdatafilefolder] = askfilepath(cfg.datafile);


if askdatafile
    answer_read = questdlg('Would you like to read in a data file?', ...
        'Read in scoring?', ...
        'Yes','No','No');
    if istrue(answer_read)
        if hasdatafilefolder
            tmpfilepath = cfg.datafile;
        end
        
        [dataset_file_name dataset_file_path dataset_file_filterindex] = uigetfile(...
            {'*.eeg;*.dat;*.edf;*.zip','Import formats (*.eeg,*.dat,*.edf)';...
            '*.eeg;*.dat','Brainvision files (*.eeg,*.dat)';...
            '*.edf','EDF files (*.edf)';...
            % '*.m', 'program files (*.m)';...
            % '*.fig','Figures (*.fig)';...
            % '*.mat','MAT-files (*.mat)';...
            '*.*',  'All other dataformats (*.*)'},...
            'Import Dataset',tmpfilepath);
        
        if isequal(dataset_file_name,0)
            cfg.datafile = [];
            skipNonData = true;
        else
            cfg.datafile = [dataset_file_path dataset_file_name];
            [temppath, tempname, tempext] = fileparts(cfg.datafile);
            tmpfilepath = temppath;
            if strcmp(tempext,'.eeg')
                cfg.headerfile = [temppath filesep tempname '.vhdr'];
            elseif strcmp(tempext,'.edf')
                cfg.headerfile = cfg.datafile;
            end
        end
    else
        cfg.datafile = [];
        skipNonData = true;
    end
else
    if strcmp(cfg.datafile,'no')
        cfg.datafile = [];
    end
end

if skipNonData
    if strcmp(cfg.headerfile,'no') || strcmp(cfg.headerfile,'ask')
        cfg.headerfile = [];
    end
    if strcmp(cfg.montagefile,'no') || strcmp(cfg.montagefile,'ask')
        cfg.montagefile = [];
    end   
    if strcmp(cfg.datagrammerfile,'no') || strcmp(cfg.datagrammerfile,'ask')
        cfg.datagrammerfile = [];
    end   
    if strcmp(cfg.scoringfile,'no') || strcmp(cfg.scoringfile,'ask')
        cfg.scoringfile = [];
    end
 
else
    
[askheaderfile hasheaderfilefolder] = askfilepath(cfg.headerfile);

if askheaderfile
    answer_read = questdlg('Would you like to read in a header file?', ...
        'Read in scoring?', ...
        'Yes','No','No');
    if istrue(answer_read)
        if hasheaderfilefolder
            tmpfilepath = cfg.headerfile;
        end
        [datasetheader_file_name datasetheader_file_path datasetheader_file_filterindex] = uigetfile(...
            {'*.vhdr;*.edf','Import formats (*.vhdr,*.edf)';...
            '*.vhdr','Brainvision files (*.vhdr)';...
            '*.edf','EDF files (*.edf)';...
            % '*.m', 'program files (*.m)';...
            % '*.fig','Figures (*.fig)';...
            % '*.mat','MAT-files (*.mat)';...
            '*.*',  'All Files (*.*)'},...
            'Import data header',tmpfilepath);
        if isequal(datasetheader_file_name,0)
            cfg.headerfile = [];
        else
            cfg.headerfile = [datasetheader_file_path datasetheader_file_name];
            [temppath, tempname, tempext] = fileparts(cfg.headerfile);
            tmpfilepath = temppath;
        end
    else
        cfg.headerfile = [];
    end
else
    if strcmp(cfg.headerfile,'no')
        cfg.headerfile = [];
    end
end


[askmontagefile hasmontagefilefolder] = askfilepath(cfg.montagefile);

if askmontagefile
    answer_read = questdlg('Would you like to read in a montage file?', ...
        'Read in scoring?', ...
        'Yes','No','No');
    
    if istrue(answer_read)
        if hasmontagefilefolder
            tmpfilepath = cfg.montagefile;
        end
        
        [montage_file_name montage_file_path montage_file_filterindex] = uigetfile(...
            {'*.txt;*.tsv;*.csv','Import formats (*.txt,*.tsv,*.csv)';...
            '*.txt','Text - Tab delimited (*.txt)';...
            '*.tsv','Text - Tab delimited (*.tsv)';...
            '*.csv','Text - Tab delimited (*.csv)';...
            % '*.m', 'program files (*.m)';...
            % '*.fig','Figures (*.fig)';...
            '*.*',  'All Files (*.*)'},...
            'Import montage',tmpfilepath);
        
        if isequal(montage_file_name,0)
            cfg.montagefile = [];
        else
            cfg.montagefile = [montage_file_path montage_file_name];
            [temppath, tempname, tempext] = fileparts(cfg.montagefile);
            tmpfilepath = temppath;
        end
    else
        cfg.montagefile = [];
    end
else
    if strcmp(cfg.montagefile,'no')
        cfg.montagefile = [];
    end
end

[askdatagrammerfile hasdatagrammerfilefolder] = askfilepath(cfg.datagrammerfile);

if askdatagrammerfile
    answer_read = questdlg('Would you like to read in a data grammer file?', ...
        'Read in scoring?', ...
        'Yes','No','No');
    
    if istrue(answer_read)
        if hasdatagrammerfilefolder
            tmpfilepath = cfg.datagrammerfile;
        end
        
        [datagrammer_file_name datagrammer_file_path datagrammer_file_filterindex] = uigetfile(...
            {'*.txt;*.tsv;*.csv','Import formats (*.txt,*.tsv,*.csv)';...
            '*.txt','Text - Tab delimited (*.txt)';...
            '*.tsv','Text - Tab delimited (*.tsv)';...
            '*.csv','Text - Tab delimited (*.csv)';...
            % '*.m', 'program files (*.m)';...
            % '*.fig','Figures (*.fig)';...
            '*.*',  'All Files (*.*)'},...
            'Import data grammer',tmpfilepath);
        
        if isequal(datagrammer_file_name,0)
            cfg.datagrammerfile = [];
        else
            cfg.datagrammerfile = [datagrammer_file_path datagrammer_file_name];
            [temppath, tempname, tempext] = fileparts(cfg.datagrammerfile);
            tmpfilepath = temppath;
        end
    else
        cfg.datagrammerfile = [];
    end
else
    if strcmp(cfg.datagrammerfile,'no')
        cfg.datagrammerfile = [];
    end
end

[askscoringfile hasscoringfilefolder] = askfilepath(cfg.scoringfile);

if askscoringfile
    answer_read = questdlg('Would you like to read in a scoring file?', ...
        'Read in scoring?', ...
        'Yes','No','No');
    
    if istrue(answer_read)
        if hasscoringfilefolder
            tmpfilepath = cfg.scoringfile;
        end
        
        [hyp_file_name hyp_file_path hyp_file_filterindex] = uigetfile(...
            {'*.txt;*.csv;*.tsv','Import formats (*.txt,*.csv,*.tsv)';...
            '*.txt','Text - Tab delimited (*.txt)';...
            '*.tsv','Text - Tab delimited (*.tsv)';...
            '*.csv','Comma Separated Values (*.csv)';...
            % '*.m', 'program files (*.m)';...
            % '*.fig','Figures (*.fig)';...
            '*.mat','MAT-files from SleepTrip scoring export(*.mat)';...
            '*.*',  'All Files (*.*)'},...
            'Import scoring',tmpfilepath);
        
        list_formats = {'SpiSOP/Schlafaus/sleepin','Zmax','Somnomedics English', 'FASST','SleepTrip'};
        list_formats_st = {'spisop','zmax','somnomedics_english','fasst','sleeptrip',};
        [indx_file_formats, selected_file_format] = listdlg('ListString',list_formats,'SelectionMode','single','PromptString',{'Select a file format.',''},'InitialValue',1,'Name','File format');
        
        cfg.scoringformat = list_formats_st{indx_file_formats};
        
        list_standards = {'AASM','Rechtschaffen&Kales'};
        list_standards_st = {'aasm','rk'};
        [indx_file_standards, selected_file_standard] = listdlg('ListString',list_standards,'SelectionMode','single','PromptString',{'Select a scoring standard.',['Scoring will be interpreted as this standard'],''},'InitialValue',1,'Name','Scoring standard');
        
        cfg.scoringstandard = list_standards_st{indx_file_standards};
        if isequal(hyp_file_name,0)
            cfg.scoringfile = [];
        else
            cfg.scoringfile = [hyp_file_path hyp_file_name];
            [temppath, tempname, tempext] = fileparts(cfg.scoringfile);
            tmpfilepath = temppath;
        end
    else
        cfg.scoringfile = [];
    end
else
    if strcmp(cfg.scoringfile,'no')
        cfg.scoringfile = [];
    end
end



end

% do the general cleanup and bookkeeping at the end of the function
ft_postamble debug
ft_postamble trackconfig
% ft_postamble previous data
% data = datanew
% ft_postamble provenance data
% ft_postamble history    data
% ft_postamble savevar    data



fprintf([functionname ' function finished\n']);
toc(ttic)
memtoc(mtic)
end

function [askfile hasfilefolder] = askfilepath(tocheckfile)
askfile = false;
hasfilefolder = false;

if isempty(tocheckfile)
    askfile = true;
else
    if strcmp(tocheckfile,'ask') %% might be a path
        askfile = true;
    else
        if ~strcmp(tocheckfile,'no')
            ex = exist(tocheckfile);
            if (ex==7)
                hasfilefolder = true;
                askfile = true;
            elseif (ex==2)
                askfile = false;
            else
                askfile = true;
            end
        end
    end
end
end
