function [input] = st_apply_datagrammer(input, datagrammer, varargin)
% ST_APPLY_DATAGRAMMER applies a data grammer 
% with filtering and data transfer functions to data. 
% It is a wrapper function of ST_DATAGRAMMER
%
% Use as
%   [data]    = ft_apply_datagrammer(data, datagrammer,  ...)
%
% A datagrammer is specified as a structure with the fields
%   datagrammer.channel = Nx1 cell-array
%   datagrammer.grammer = Nx1 cell-array
%
% As an example, a datagrammer with different filterings for different 
% channels could look like this
%   datagrammer.channel  = {'EMG',   'F*'}
%   datagrammer.grammer  = {'hp 10', 'hp 0.5 lp 35'}
%   see ST_DATAGRAMMER for details on what datagrammer.grammer can be
%
% Additional options should be specified in key-value pairs and can be
%   'feedback'      = string, see FT_PROGRESS (default = 'text')
%   'warning'       = boolean, whether to show warnings (default = true)
%   'showcallinfo'  = string, 'yes' or 'no' (default = 'no')
%
% See also ST_DATAGRAMMER, ST_PREPROCESSING

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

if iscell(input) && iscell(input)
  % this represents combined EEG, ECoG and/or MEG
  for i=1:numel(input)
    input{i} = st_apply_datagrammer(input{i}, datagrammer, varargin{:});
  end
  return
end

% get optional input arguments
%keepunused    = ft_getopt(varargin, 'keepunused',  'no');
%inverse       = ft_getopt(varargin, 'inverse',     'no');
feedback      = ft_getopt(varargin, 'feedback',    'text');
showwarning   = ft_getopt(varargin, 'warning',     true);
%bname         = ft_getopt(varargin, 'balancename', '');
showcallinfo  = ft_getopt(varargin, 'showcallinfo', 'no');

if istrue(showwarning)
  warningfun = @warning;
else
  warningfun = @nowarning;
end

% check the consistency of the datagrammer
if numel(datagrammer.channel)~=numel(datagrammer.grammer)
  ft_error('the number of channels in the datagrammer is inconsistent with the number of grammers');
end

if isfield(input, 'chantype') && ~isequal(input.chantype, datagrammer.chantypeold)
  ft_error('inconsistent chantype in data and montage');
elseif isfield(input, 'chantypenew') && ~isequal(input.chantypenew, datagrammer.chantypeold)
  ft_error('inconsistent chantype in data and montage');
end

if isfield(input, 'channel') && isfield(input, 'grammer')
  inputtype = 'datagrammer';
elseif isfield(input, 'tra')
  inputtype = 'sens';
elseif ft_datatype(input, 'raw')
  inputtype = 'raw';
elseif ft_datatype(input, 'timelock')
  inputtype = 'timelock';
elseif isfield(input, 'fourierspctrm')
  inputtype = 'freq';
else
  inputtype = 'unknown';
end

switch inputtype
 % case 'datagrammer'
    
  case 'raw'
    % apply the montage to the raw data that was preprocessed using FieldTrip
    data = input;
    clear input
    for iChannels = 1:numel(datagrammer.channel)
        cfg_dg = [];
        cfg_dg.channel = datagrammer.channel(iChannels);
        cfg_dg.grammer = datagrammer.grammer(iChannels);
        data = st_datagrammer(cfg_dg, data);
    end

    % rename the output variable
    input = data;
    clear data

 % case 'timelock'
    
 % case 'freq'
    
  otherwise
    ft_error('unsupported input format');
end % switch inputtype



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % HELPER FUNCTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = indx2logical(x, n)
% y = false(1,n);
% y(x) = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nowarning(varargin)
return
