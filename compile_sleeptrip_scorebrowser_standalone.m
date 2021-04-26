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


%% minimal example pipeline for a sleeptrip analysis

%%% PLEASE MAKE SURE:
%%% 1. YOU HAVE Matlab2013b or later
%%% 2. YOU HAVE ALL REQUIRED TOOLBOXES (typically the case)
%%% 3. YOUR SLEEPTRIP FOLDER DIRECTLY HAS THE FILES IN THEM, 
%%%    e.g. while extracting a zip from github maybe you need to refer to 
%%%    D:/sleeptrip-master/sleeptrip-master instead of D:/sleeptrip-master
%%% 4. NO FieldTrip NOR related toolboxes like EEGLAB are loaded or in the 
%%%    PATH of your MATLAB
%%% 5. You have downloaded the example datasets here:
%%%    https://drive.google.com/open?id=160i_FG3zSZrBf9mr2zxZQX-u30OnLDMf
%%% 6. ...and the files from the example datasets are in the 
%%%    same folder as this script or somewhere in the matlab path
   


%% initialization

%%% add sleeptrip to the path
%addpath('D:/sleeptrip')
%cd('D:/sleeptrip')
pathToSleepTrip = uigetdir('','choose spisop path, e.g. D:\sleeptrip');
addpath(pathToSleepTrip)
cd(pathToSleepTrip)


%%% mabye disable some toolboxes if necessary and check if 
%%% 'signal_toolbox', 'signal_blocks' are available 
%%% because they are helpful to have.
% toggleToolbox('names')
% toggleToolbox('dsp','off')
% toggleToolbox('all','query')
% license('inuse')

%%% clear and load all the defaults of SleepTrip and FieldTrip
clear ft_defaults;
clear st_defaults;
st_defaults;

% do not use my personal defaults, but rather FieldTrip standard defaults
global st_default
st_default = [];
global ft_default
ft_default = [];

%check if matlab compiler is installed, licencded
%license('checkout','Compiler')
%ans =
%     1
%chec connection to the lic sever is ok:
%!mcc
mcc -mv sleeptrip_scorebrowser_standalone.m

delete mccExcludedFiles.log
delete requiredMCRProducts.txt
delete readme.txt