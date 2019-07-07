function filt = fw_filter_with_correction(B,A,dat,dir,hd,type,isFT)

% FILTER_WITH_CORRECTION applies a to the data and corrects
% edge-artifacts for one-pass filtering.
%
% Use as
%   [filt] = filter_with_correction(B,A,dat,dir);
% where
%   B,A        filter coefficients
%   dat        data matrix (Nchans X Ntime)
%   dir        optional filter direction, can be
%                'onepass'         forward filter only
%                'onepass-reverse' reverse filter only, i.e. backward in time
%                'twopass'         zero-phase forward and reverse filter (default)
%                'twopass-reverse' zero-phase reverse and forward filter
%                'twopass-average' average of the twopass and the twopass-reverse
%
% Note that a one- or two-pass filter has consequences for the
% strength of the filter, i.e. a two-pass filter with the same filter
% order will attenuate the signal twice as strong.

% Copyright (c) 2010, Stefan Klanke
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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
% $Id: filter_with_correction.m 7513 2013-02-20 16:54:20Z roboos $

if ~isstable(hd)
    error('filter unstable, calculated filter coefficients have poles on or outside the unit circle and will not be stable. Try a higher cutoff frequency or a different type/order of filter.');
end


poles = roots(A);
if any(abs(poles) >= 1)
  error('Calculated filter coefficients have poles on or outside the unit circle and will not be stable. Try a higher cutoff frequency or a different type/order of filter.');
end

dcGain = sum(B)/sum(A);

[d,N] = size(dat);

filt = [];
if isFT
    switch dir
        case 'onepass'
            offset = dat(:,1);
            dat = dat - repmat(offset,1,N);
            filt = filter(B, A, dat')' + repmat(dcGain*offset, 1, N);
        case 'onepass-reverse'
            offset = dat(:,end);
            dat  = fliplr(dat) - repmat(offset,1,N);
            filt = filter(B, A, dat')';
            filt = fliplr(filt) + repmat(dcGain*offset, 1, N);
        case 'twopass'
            % filtfilt does the correction for us
            filt = ft_fw_filtfilt(B, A, dat')';
        case 'twopass-reverse'
            % filtfilt does the correction for us
            filt = fliplr(ft_fw_filtfilt(B, A, fliplr(dat)')');
        case 'twopass-average'
            % take the average from the twopass and the twopass-reverse
            filt1 = ft_fw_filtfilt(B, A, dat')';
            filt2 = fliplr(ft_fw_filtfilt(B, A, fliplr(dat)')');
            filt  = (filt1 + filt2)/2;
        otherwise
            error('unsupported filter direction "%s"', dir);
    end
else
    
    switch dir
          case 'onepass'
            offset = dat(:,1);
            dat = dat - repmat(offset,1,N);
            filt = filter(B, A, dat')' + repmat(dcGain*offset, 1, N);
          case 'onepass-reverse'
            offset = dat(:,end);
            dat  = fliplr(dat) - repmat(offset,1,N);
            filt = filter(B, A, dat')';
            filt = fliplr(filt) + repmat(dcGain*offset, 1, N);
        case 'twopass'
            % filtfilt does the correction for us
            if strcmp(type,'IIRdesigned')
                filt = filtfilt(hd.sosMatrix,hd.ScaleValues, dat')';
            elseif strcmp(type,'FIRdesigned')
                filt = filtfilt(hd.Numerator,1, dat')';
            end
        case 'twopass-reverse'
            % filtfilt does the correction for us
            if strcmp(type,'IIRdesigned')
                filt = fliplr(filtfilt(hd.sosMatrix,hd.ScaleValues, fliplr(dat)')');
            elseif strcmp(type,'FIRdesigned')
                filt = fliplr(filtfilt(hd.Numerator,1, fliplr(dat)')');
            end
        case 'twopass-average'
            % take the average from the twopass and the twopass-reverse
            if strcmp(type,'IIRdesigned')
                filt1 = filtfilt(hd.sosMatrix,hd.ScaleValues, dat')';
                filt2 = fliplr(filtfilt(hd.sosMatrix,hd.ScaleValues, fliplr(dat)')');
                filt  = (filt1 + filt2)/2;
            elseif strcmp(type,'FIRdesigned')
                filt1 = filtfilt(hd.Numerator,1, dat')';
                filt2 = fliplr(filtfilt(hd.Numerator,1, fliplr(dat)')');
                filt  = (filt1 + filt2)/2;
            end
            
        otherwise
            error('unsupported filter direction "%s"', dir);
    end
    
end
end
