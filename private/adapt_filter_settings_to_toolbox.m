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

baseFirOrder100Hz = 442;
firfilterorder = (fsample/100)*baseFirOrder100Hz;
if mod(fix(firfilterorder),2) %% uneven
    firfilterorder = fix(firfilterorder) + 1;
else %% even
    firfilterorder = fix(firfilterorder);
end

if use_hp && (strcmp(core_cfg.hpfilttype,'IIRdesigned') || strcmp(core_cfg.hpfilttype,'FIRdesigned'))
    if ~license('checkout','signal_blocks') || strcmp(use_signal_blocks,'no')
        ft_warning('you are missing or not using the ''signal_blocks'' toolbox license or it is currently in use by someone else.\n... will try using ''signal_toolbox'' instead and ''firls'' or ''but'' design.')
        if strcmp(core_cfg.hpfilttype,'IIRdesigned')
            core_cfg.hpfilttype = 'but';
            if ~strcmp(UseFixedFilterOrder_hp,'yes')
                ft_warning('using FieldTrip default filter order for but which is probably 6.')
                core_cfg.hpfiltord = [];
            end
        elseif strcmp(core_cfg.hpfilttype,'FIRdesigned')
            if ~license('checkout','signal_toolbox') || strcmp(use_signal_toolbox,'no')
                ft_warning('you are missing or not using the ''signal_toolbox'' toolbox license or it is currently in use by someone else.\n ...using ''fir'' filter instead')
                core_cfg.hpfilttype = 'fir';
            else
                core_cfg.hpfilttype = 'firlsdesign';
                %core_cfg.StopToPassTransitionWidth_bp = StopToPassTransitionWidth_bp;% frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
                %core_cfg.PassToStopTransitionWidth_bp = PassToStopTransitionWidth_bp;% frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25
                %core_cfg.StopToPassTransitionWidth_hp = StopToPassTransitionWidth_hp;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2
                core_cfg.StopToPassTransitionWidth_hp = StopToPassTransitionWidth_hp_temp;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.1
                %core_cfg.PassToStopTransitionWidth_lp = PassToStopTransitionWidth_lp;
            end
            if ~strcmp(UseFixedFilterOrder_hp,'yes')
                core_cfg.hpfiltord = firfilterorder;
            end
        end
    end
end

if use_lp && (strcmp(core_cfg.lpfilttype,'IIRdesigned') || strcmp(core_cfg.lpfilttype,'FIRdesigned'))
    if ~license('checkout','signal_blocks') || strcmp(use_signal_blocks,'no')
        ft_warning('you are missing or not using the ''signal_blocks'' toolbox license or it is currently in use by someone else.\n... will try using ''signal_toolbox'' instead and ''firls'' or ''but'' design.')
        if strcmp(core_cfg.lpfilttype,'IIRdesigned')
            core_cfg.lpfilttype = 'but';
            if ~strcmp(UseFixedFilterOrder_lp,'yes')
                ft_warning('using FieldTrip default filter order for but which is probably 6.')
                core_cfg.lpfiltord = [];
            end
        elseif strcmp(core_cfg.lpfilttype,'FIRdesigned')
            if ~license('checkout','signal_toolbox') || strcmp(use_signal_toolbox,'no')
                ft_warning('you are missing or not using the ''signal_toolbox'' toolbox license or it is currently in use by someone else.\n ...using ''fir'' filter instead')
                core_cfg.lpfilttype = 'fir';
            else
                core_cfg.lpfilttype = 'firlsdesign';
                %core_cfg.StopToPassTransitionWidth_bp = StopToPassTransitionWidth_bp;% frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
                %core_cfg.PassToStopTransitionWidth_bp = PassToStopTransitionWidth_bp;% frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25
                %core_cfg.StopToPassTransitionWidth_hp = StopToPassTransitionWidth_hp;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2
                %core_cfg.StopToPassTransitionWidth_hp = StopToPassTransitionWidth_hp_temp;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.1
                core_cfg.PassToStopTransitionWidth_lp = PassToStopTransitionWidth_lp;
            end
            if ~strcmp(UseFixedFilterOrder_lp,'yes')
                core_cfg.lpfiltord = firfilterorder;
            end
        end
    end
end

if use_bp && (strcmp(core_cfg.bpfilttype,'IIRdesigned') || strcmp(core_cfg.bpfilttype,'FIRdesigned'))
    if ~license('checkout','signal_blocks') || strcmp(use_signal_blocks,'no')
        ft_warning('you are missing or not using the ''signal_blocks'' toolbox license or it is currently in use by someone else.\n... will try using ''signal_toolbox'' instead and ''firls'' or ''but'' design.')
        if strcmp(core_cfg.bpfilttype,'IIRdesigned')
            core_cfg.bpfilttype = 'but';
            if ~strcmp(UseFixedFilterOrder_bp,'yes')
                ft_warning('using FieldTrip default filter order for but which is probably 6.')
                core_cfg.bpfiltord = [];
            end
        elseif strcmp(core_cfg.bpfilttype,'FIRdesigned')
            if ~license('checkout','signal_toolbox') || strcmp(use_signal_toolbox,'no')
                ft_warning('you are missing or not using the ''signal_toolbox'' toolbox license or it is currently in use by someone else.\n ...using ''fir'' filter instead')
                core_cfg.bpfilttype = 'fir';
            else
                core_cfg.bpfilttype = 'firlsdesign';
                core_cfg.StopToPassTransitionWidth_bp = StopToPassTransitionWidth_bp;% frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
                core_cfg.PassToStopTransitionWidth_bp = PassToStopTransitionWidth_bp;% frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25
                %core_cfg.StopToPassTransitionWidth_hp = StopToPassTransitionWidth_hp;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2
                %core_cfg.StopToPassTransitionWidth_hp = StopToPassTransitionWidth_hp_temp;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.1
                %core_cfg.PassToStopTransitionWidth_lp = PassToStopTransitionWidth_lp;
            end
            if ~strcmp(UseFixedFilterOrder_bp,'yes')
                core_cfg.bpfiltord = firfilterorder;
            end
        end
    end
end

