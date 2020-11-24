function [dat, hdr] = read_axograph_data(filname, read_data)
% read in the header of Axograph files and if read_data is true also the
% data

% Copyright (C)   2019-,   Frederik D. Weber
% with code snipets from   Marco J Russo & Jeffrey Gill
% from https://github.com/inquietus/importaxo & https://github.com/CWRUChielLab/importaxographx
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

% .axgx file format:
%
% Byte		Format		Content
% 1-4		<char>		File type identifier (e.g. 'axgx');
% 5-8		<int32>		File version number (e.g. '6')
% 9-12		<int32>		Number of columns
% 13-16		<int32>		Number of elements in column
% 17-20		<int32>		Column type (see below)
% 21-24		<int32>		Number of bytes in column title (2*characters)
% 25-39		<char>		Column title, e.g. 'Time (s)' (each char preceded by null byte)
% 40-47		<double>	Sample interval (8-byte double)
% 48-55		<double>	Sample interval (8-byte double) Redundant?
% 56-59		<int32>		Number of elements in column
% 60-61		<int32>		Column type (see below)
% 62-63		<int32>		Number of bytes in column title (2*characters)
% 64-...	<char>		Column title
% ...		<datatype>	Data

fid = fopen(filname,'r','ieee-le.l64');

hdr.fileName   = filname;
hdr.fileType   = fread(fid, 4, '*char')';
hdr.fileVer    = fread(fid,1,'int32');

if hdr.fileVer - swapbytes(hdr.fileVer) >= 0
    %File is big-endian, close and reopen
    fclose(fid);
    fid = fopen(filname,'r','ieee-be.l64');
    hdr.fileType = fread(fid,4,'*char');
    hdr.fileVer = fread(fid,1,'int32');
end

hdr.NDatCol    = fread(fid,1,'int32');
hdr.nPoints    = fread(fid,1,'int32');
hdr.XColType   = fread(fid,1,'int32');
hdr.XTitleLen  = fread(fid,1,'int32');
hdr.XTitle	   = readtitle(fid,hdr.XTitleLen);
hdr.SampInt    = fread(fid,1,'double');
hdr.SampInt2   = fread(fid,1,'double'); %Redundant in data file (2nd channel?)

hdr.fsample    = hdr.SampInt;
hdr.nChannel   = hdr.NDatCol-1;
hdr.time       = double(1:hdr.nPoints).*hdr.SampInt;

dat = [];
if read_data
    dat = double(zeros(hdr.NDatCol-1,hdr.nPoints));
end

for iYCol = 1:(hdr.NDatCol-1)
    hdr.YCol(iYCol,1).nPoints   = fread(fid,1,'int32');
    hdr.YCol(iYCol,1).YColType  = fread(fid,1,'int32');
    hdr.YCol(iYCol,1).YTitleLen = fread(fid,1,'int32');
    hdr.YCol(iYCol,1).YTitle    = readtitle(fid, hdr.YCol(iYCol).YTitleLen);
    
    %%%%%%%%%%%%%%% REMOVE NULL CHARACTERS %%%%%%%%%%%%%%%
    % It seems that the types are not set correctly for  %
    % the version of AxoGraph X that I have because      %
    % titlelen is double what it should be and title has %
    % null characters inserted before each real          %
    % character. This code corrects these errors.        %
    hdr.YCol(iYCol).YTitleLen = hdr.YCol(iYCol).YTitleLen / 2;
    %hdr.YCol(iYCol).YTitle = hdr.YCol(iYCol).YTitle(2:2:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if read_data
        switch hdr.YCol(iYCol).YColType
            case 4  % column type is short
                dat(iYCol,:) = double(fread(fid,hdr.YCol(iYCol).nPoints, 'int16'));
            case 5  % column type is long
                dat(iYCol,:) = double(fread(fid,hdr.YCol(iYCol).nPoints, 'int32'));
            case 6  % column type is float
                dat(iYCol,:) = double(fread(fid,hdr.YCol(iYCol).nPoints, 'float32'));
            case 7  % column type is double
                dat(iYCol,:) = fread(fid,hdr.YCol(iYCol).nPoints, 'double');
            case 9  % 'series'
                data0 = fread(fid,2, 'double');
                %data(:, iYCol) = (1:1:hdr.YCol(iYCol).nPoints)*data0(2)+ data0(1);
                %%%%%%%%%%%%%%%%%%%%%%%% FIX SERIES START %%%%%%%%%%%%%%%%%%%%%%%%%
                % data0(1) is the starting value, data0(2) is the interval, so    %
                % the code above is off by one. This code corrects the error.     %
                dat(iYCol,:) = (0:1:hdr.YCol(iYCol).nPoints-1)*data0(2)+ data0(1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 10 % 'scaled short'
                scale = fread(fid,1, 'double');
                offset = fread(fid,1, 'double');
                data0 = fread(fid,hdr.YCol(iYCol).nPoints, 'int16');
                dat(iYCol,:)  = double(data0)*scale + offset;
            otherwise
        end
    else
        switch hdr.YCol(iYCol).YColType
            case 4  % column type is short
                fread(fid,hdr.YCol(iYCol).nPoints, 'int16');
            case 5  % column type is long
                fread(fid,hdr.YCol(iYCol).nPoints, 'int32');
            case 6  % column type is float
                fread(fid,hdr.YCol(iYCol).nPoints, 'float32');
            case 7  % column type is double
                fread(fid,hdr.YCol(iYCol).nPoints, 'double');
            case 9  % 'series'
                fread(fid,2, 'double');
            case 10 % 'scaled short'
                fread(fid,1, 'double');
                fread(fid,1, 'double');
                fread(fid,hdr.YCol(iYCol).nPoints, 'int16');
            otherwise
        end
    end
end

hdr.channels = {};
for iCh = 1:numel(hdr.YCol)
    hdr.channels{iCh} = hdr.YCol(iCh).YTitle;
end
hdr.channels = hdr.channels';

%Acquisition settings - from additional ~1700 bytes
ft = fread(fid,'int16=>char')';
infoStart = regexp(ft,'--- Acquisition Settings ---','once');
infoEnd = regexp(ft,'Inter-Pulse Interval Table Entries\s[\d+]','end','once');
hdr.Info = ft(infoStart:infoEnd);

hdr.Protocol = regexp(hdr.Info,'Protocol\s:\s(\w+)\s','tokens','once');

%Get acquisition channels
hdr.additionalChannelInfo = regexp(hdr.Info,'Acquisition channels : "(.+)"','tokens','once');

%Get number of episodes
hdr.NumEpisodes = regexp(hdr.Info,'Episodes\s(\d+)\s','tokens','once');
if ~isempty(hdr.NumEpisodes)
    hdr.NumEpisodes = str2num(hdr.NumEpisodes{:});
end
%Get number of pulses
hdr.NumPulses = regexp(hdr.Info,'Pulses\s(\d+)\s','tokens','once');
if ~isempty(hdr.NumPulses)
    hdr.NumPulses = str2num(hdr.NumPulses{:});
end
%Get pulse tables for protocol reconstruction
%NOTE: Exact structure and meaning of protocol table unclear
expProt = '(\s-?\d?\.?\d+){10,16}';
hdr.ProtocolTable = regexp(ft,expProt,'tokens');
if ~isempty(hdr.ProtocolTable)
    hdr.ProtocolTable = cellfun(@(x) str2num(char(x)),hdr.ProtocolTable,'UniformOutput',false);
end

fclose(fid);

end

function title = readtitle(fid,titlelen)
str = fread(fid,titlelen,'*char')';
title = str(~mod(1:titlelen,2));
end
