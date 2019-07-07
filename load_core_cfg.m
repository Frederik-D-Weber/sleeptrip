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

use_signal_blocks = 'no'; %either 'no', 'yes'
use_signal_toolbox = 'yes'; %either 'no', 'yes'

core_cfg = [];

core_cfg.feedback = 'no';% if there is interupting feedback from the FieldTrip toolbox functions though not totally avoidable either no or text or textbar or gui default no
core_cfg.precision = 'double';% either single or double default double
core_cfg.dftfilter = 'no';% line noise removal using discrete fourier transform either no or yes default yes 
core_cfg.dftfreq = [50 100 150];% line noise frequencies in Hz for DFT filter default 50 100 150
core_cfg.bpfilttype = 'FIRdesigned';% digital filter type for band pass filtering. Also but for fieldtrip Butterworth IIR filter or fir for fieldtrip FIR filter using Matlab fir1 function is supported. Note that in case of fir the Fpass*** frequencies are then the -3dB cutt-off frequencies default FIRdesigned
core_cfg.bpfiltdir = 'twopass';% filter direction either twopass or twopass-reverse or twopass-average default twopass
core_cfg.bpinstabilityfix = 'no';% deal with filter instability either no for only detect and give error or reduce for reduce the filter order or split for split the filter in two lower-order filters to apply sequentially default no
UseFixedFilterOrder_bp = 'no';% ignore attenuation values (i.e. Apass AstopLeft AstopRight) and use a fixed filter order for band pass filtering either yes or no default no
FilterOrder_bp = 400;% in case of UseFixedFilterOrder_bp is set to yes then this order is used for band pass filtering 
UseTwoPassAttenuationCorrection_bp = 'no';% if in case of a two-pass fitlering the attenuation values should be halfed before each pass. This results in closer attenuation than requested. either yes or no default no
core_cfg.lpfilttype = 'FIRdesigned';% digital filter type for low pass filtering. Also but for fieldtrip Butterworth IIR filter or fir for fieldtrip FIR filter using Matlab fir1 function is supported. Note that in case of fir the Fpass*** frequencies are then the -3dB cutt-off frequencies default FIRdesigned
core_cfg.lpfiltdir = 'twopass';% filter direction either twopass or twopass-reverse or twopass-average default twopass
core_cfg.lpinstabilityfix = 'no';% deal with filter instability either no for only detect and give error or reduce for reduce the filter order or split for split the filter in two lower-order filters to apply sequentially default no
UseFixedFilterOrder_lp = 'no';% ignore attenuation values (i.e. Apass AstopLeft AstopRight) and use a fixed filter order for low pass filtering either yes or no default no
FilterOrder_lp = 400;% in case of UseFixedFilterOrder_lp is set to yes then this order is used for low pass filtering 
UseTwoPassAttenuationCorrection_lp = 'no';% if in case of a two-pass fitlering the attenuation values should be halfed before each pass. This results in closer attenuation than requested. either yes or no default no
core_cfg.hpfilttype = 'IIRdesigned';% digital filter type either IIRdesigned for Butterworth IIR filter or FIRdesigned for FIR filter using Matlab filter (one-pass) or filtfilt (two-pass) function. note. Also but for fieldtrip Butterworth IIR filter or fir for fieldtrip FIR filter using Matlab fir1 function is supported. Note that in case of fir the Fpass*** frequencies are then the -3dB cutt-off frequencies default IIRdesigned
core_cfg.hpfiltdir = 'twopass';% filter direction either twopass or twopass-reverse or twopass-average default twopass
core_cfg.hpinstabilityfix = 'no';% deal with filter instability either no for only detect and give error or reduce for reduce the filter order or split for split the filter in two lower-order filters to apply sequentially default no
UseFixedFilterOrder_hp = 'yes';% ignore attenuation values (i.e. Apass AstopLeft AstopRight) and use a fixed filter order for low pass filtering either yes or no default yes
FilterOrder_hp = 4;% in case of UseFixedFilterOrder_hp is set to yes then this order is used for high pass filtering default 4
UseTwoPassAttenuationCorrection_hp = 'no';% if in case of a two-pass fitlering the attenuation values should be halfed before each pass. This results in closer attenuation than requested. either yes or no default no
Apass = 0.001;% Attenuation of bandpass and lowpass and highpass (ripples) as one parameter for all filters in db default 0.001
AstopLeft = 100;% Attenuation of left stop band (<FstopLeft) ripples for band pass filter and high pass filter in db default 100
AstopRight = 100;% Attenuation of right stop band (>FstopRight) ripples for band pass filter and low pass filter in db default 100
StopToPassTransitionWidth_bp = 1.25;% frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 1.25
PassToStopTransitionWidth_bp = 1.25;% frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25
StopToPassTransitionWidth_hp = 0.2;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.2
StopToPassTransitionWidth_hp_predownsample = 0.1;% for high pass filter frequency in Hz that is substacted to left pass frequency (e.g. FpassLeft) to get left stop frequency (e.g. FstopLeft) default 0.1
PassToStopTransitionWidth_lp = 1.25;% for low pass filter frequency in Hz that is added to right pass frequency (e.g. FpassRight) to get right stop frequency (e.g. FstopRight) default 1.25
UseFTfiltfilt = 'no';% if to use the fieldtrip internal open source filtfilt function (i.e. yes) or the one from the matlab signal toolbox (i.e. no).
MaximizeFilterOrderIfFixedFilterOrderIsUsed = 'no';% default no

core_cfg.use_ft_filtfilt = strcmp(UseFTfiltfilt,'yes');

useTwoPassFiltering_bp = 'no';

useTwoPassFiltering_lp = 'no';

useTwoPassFiltering_hp = 'no';

if ~isempty(strfind(core_cfg.bpfiltdir,'two'))
    useTwoPassFiltering_bp = 'yes';
end

if ~isempty(strfind(core_cfg.lpfiltdir,'two'))
   useTwoPassFiltering_lp = 'yes';
end

if ~isempty(strfind(core_cfg.hpfiltdir,'two'))
    useTwoPassFiltering_hp = 'yes';
end

