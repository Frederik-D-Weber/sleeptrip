function [dat] = read_brainvision_eeg(filename, hdr, begsample, endsample, chanindx)

% READ_BRAINVISION_EEG reads raw data from an EEG file
% and returns it as a Nchans x Nsamples matrix
%
% Use as
%   dat = read_brainvision_eeg(filename, hdr, begsample, endsample)
% where the header should be first read using read_brainvision_vhdr
%
% See also READ_BRAINVISION_VHDR, READ_BRAINVISION_VMRK

% Copyright (C) 2003-2011, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<5
    % read all channels
    chanindx = [];
end

if isequal(chanindx(:)', 1:hdr.NumberOfChannels);
    % read all channels
    chanindx = [];
end

% FIXME it would be nice to also implement the efficient reading of the
% selected channels for the other file formats but at the moment only the
% implementation of the binary multiplexed and vectorized formats is smart enough.

if strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'multiplexed') && any(strcmpi(hdr.BinaryFormat, {'int_16', 'int_32', 'ieee_float_32'}))
    
    switch lower(hdr.BinaryFormat)
        case 'int_16'
            sampletype = 'int16';
            samplesize = 2;
        case 'int_32'
            sampletype = 'int32';
            samplesize = 4;
        case 'ieee_float_32'
            sampletype = 'float32';
            samplesize = 4;
    end % case
    
    fid = fopen_or_error(filename, 'rb', 'ieee-le');
    
    if isempty(chanindx)
        % read all the channels
        fseek(fid, hdr.NumberOfChannels*samplesize*(begsample-1), 'cof');
        dat = fread(fid, [hdr.NumberOfChannels, (endsample-begsample+1)], [sampletype '=>float32']);
    else
        numsamples = (endsample-begsample+1);
        % read only the selected channels
        dat = zeros(length(chanindx), numsamples,'single');
        for chan = length(chanindx):-1:1
            fseek(fid, hdr.NumberOfChannels*samplesize*(begsample-1) + (chanindx(chan)-1)*samplesize, 'bof');
            dat(chan,:) = fread(fid, [1, numsamples], [sampletype '=>float32'], (hdr.NumberOfChannels-1)*samplesize);
        end
        
        
        
    end
    
    fclose(fid);
    
    % compute real microvolts using the calibration factor (resolution)
    % calib = diag(hdr.resolution(chanindx));
    % % using a sparse multiplication speeds it up
    % dat = full(sparse(calib) * dat);
    if isempty(chanindx)
        chanindxres = 1:size(dat,1);
    else
        chanindxres = chanindx;
    end
    if ~all(hdr.resolution(chanindxres) == 1) && ~all(isnan(hdr.resolution(chanindxres)))
        calib = reshape(hdr.resolution(chanindxres),[],1);
        %this is slower, and sparse multiplication is not supported with
        %single input:
        % dat = diag(calib) * dat;
        %this is faster
        for k = 1:size(dat,1)
            dat(k,:) = calib(k) .* dat(k,:);
        end
        
    end
    chanindxres = [];

    % don't do the channel selection again at the end of the function
    chanindx = [];
elseif strcmpi(hdr.DataFormat, 'binary') && strcmpi(hdr.DataOrientation, 'vectorized') && strcmpi(hdr.BinaryFormat, 'ieee_float_32')
    switch lower(hdr.BinaryFormat)
        case 'int_16'
            sampletype = 'int16';
            samplesize = 2;
        case 'int_32'
            sampletype = 'int32';
            samplesize = 4;
        case 'ieee_float_32'
            sampletype = 'float32';
            samplesize = 4;
    end % case
    
    fid = fopen_or_error(filename, 'rb', 'ieee-le');
    fseek(fid, 0, 'eof');
    hdr.nSamples = ftell(fid)/(samplesize*hdr.NumberOfChannels);
    
    numsamples = (endsample-begsample+1);
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % read all channels first and do NOT set chanindx = [];
    %   fseek(fid, 0, 'bof');
    %   dat = zeros(hdr.NumberOfChannels,numsamples,'single');
    %   for chan=1:hdr.NumberOfChannels
    %     fseek(fid, (begsample-1)*samplesize, 'cof');                 % skip the first N samples
    %     dat(chan,:) = fread(fid, numsamples, [sampletype '=>float32']);     % read these samples
    %     fseek(fid, (hdr.nSamples-endsample)*samplesize, 'cof');      % skip the last M samples
    %   end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %   % read only the selected channels
    %   dat = zeros(length(chanindx), numsamples,'single');
    %   for chan=1:length(chanindx)
    %    	fseek(fid, (begsample - 1) * hdr.NumberOfChannels * samplesize + (chanindx(chan) - 1) * samplesize, 'bof');
    %     dat(chan, :) = fread(fid, [1, numsamples], [sampletype '=>float32'], (hdr.NumberOfChannels - 1) * samplesize);
    %   end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% this is the fastest version and sparse in memory
    fseek(fid, 0, 'bof');
    if isempty(chanindx) && (endsample == hdr.nSamples) && (begsample == 1)% Read entire file
        dat = fread(fid, [hdr.nSamples, hdr.NumberOfChannels], [sampletype '=>float32']).';
    else
        if isempty(chanindx) % Read all channels
            chanindx = 1:hdr.NumberOfChannels;
        end
        
        dat = zeros(length(chanindx), numsamples,'single');
        ich = 0;
        for chan=1:hdr.NumberOfChannels
            fseek(fid, (begsample-1)*samplesize, 'cof');                 % skip the first N samples
            if ismember(chan,chanindx)
                ich = ich + 1;
                dat(ich,:) = fread(fid, numsamples, [sampletype '=>float32']);     % read these samples
            else
                fseek(fid, numsamples, 'cof');
            end
            fseek(fid, (hdr.nSamples-endsample)*samplesize, 'cof');      % skip the last M samples
        end
        
    end
    
    
    
    fclose(fid);
    
    % compute real microvolts using the calibration factor (resolution)
    % calib = diag(hdr.resolution(chanindx));
    % % using a sparse multiplication speeds it up
    % dat = full(sparse(calib) * dat);
    if isempty(chanindx)
        chanindxres = 1:size(dat,1);
    else
        chanindxres = chanindx;
    end
    if ~all(hdr.resolution(chanindxres) == 1) && ~all(isnan(hdr.resolution(chanindxres)))
        calib = reshape(hdr.resolution(chanindxres),[],1);
        %this is slower, and sparse multiplication is not supported with
        %single input:
        % dat = diag(calib) * dat;
        %this is faster
        for k = 1:size(dat,1)
            dat(k,:) = calib(k) .* dat(k,:);
        end
        
    end
    chanindxres = [];
    
    % don't do the channel selection again at the end of the function
    chanindx = [];
    
elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'multiplexed')
    fid = fopen_or_error(filename, 'rt');
    
    % skip lines if hdr.skipLines is set and not zero
    if isfield(hdr,'skipLines') && hdr.skipLines > 0
        for line=1:hdr.skipLines
            str = fgets(fid);
        end
    end
    
    for line=1:(begsample-1)
        % read first lines and discard the data in them
        str = fgets(fid);
    end
    dat = zeros(endsample-begsample+1, hdr.NumberOfChannels,'single');
    for line=1:(endsample-begsample+1)
        str = fgets(fid);         % read a single line with Nchan samples
        str(str==',') = '.';      % replace comma with point
        dat(line,:) = str2double(str);
    end
    fclose(fid);
    % transpose the data
    dat = dat';
    
elseif strcmpi(hdr.DataFormat, 'ascii') && strcmpi(hdr.DataOrientation, 'vectorized')
    % this is a very inefficient fileformat to read data from, since it requires to
    % read in all the samples of each channel and then select only the samples of interest
    fid = fopen_or_error(filename, 'rt');
    dat = zeros(hdr.NumberOfChannels, endsample-begsample+1);
    skipColumns = 0;
    for chan=1:hdr.NumberOfChannels
        % this is very slow, so better give some feedback to indicate that something is happening
        fprintf('reading channel %d from ascii file to get data from sample %d to %d\n', chan, begsample, endsample);
        
        % check whether columns have to be skipped
        if isfield(hdr,'skipColumns'); skipColumns = hdr.skipColumns; end
        
        str = fgets(fid);             % read all samples of a single channel
        str(str==',') = '.';          % replace comma with point
        
        if ~isempty(regexp(str(1:10), '[a-zA-Z]', 'once'))
            % the line starts with letters, not numbers: probably it is a channel label
            % find the first number and remove the preceding part
            sel   = regexp(str(1:10), ' [-0-9]');   % find the first number, or actually the last space before the first number
            label = str(1:(sel));                   % this includes the space
            str   = str((sel+1):end);               % keep only the numbers
            
            % as heading columns are already removed, set skipColumns to zero
            skipColumns = 0;
        end
        
        % convert the string into numbers and copy the desired samples over
        % into the data matrix
        tmp = str2double(str);
        dat(chan,:) = tmp(skipColumns+begsample:skipColumns+endsample);
    end
    fclose(fid);
    
else
    ft_error('unsupported sub-fileformat');
end

if ~isempty(chanindx)
    % for the the multiplexed and vectorized binary formats the channel selection was
    % already done in the code above, for the other formats the selection
    % still has to be done here
    dat = dat(chanindx,:);
end

% Convert to double, not necessary, as this is typically done higher order in ft_preprocessing
% ...for MATLAB < R14
%if str2double(version('-release')) < 14
%if ~isa(dat,'double')
%    dat = double(dat);
%end
end

