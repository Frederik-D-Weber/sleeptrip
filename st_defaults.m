function st_defaults

% FT_DEFAULTS (ending with "s") sets some general settings in the global variable
% ft_default (without the "s") and takes care of the required path settings. You can
% call this function in your startup.m script. This function is also called at the
% begin of all FieldTrip functions.
%
% The global configuration defaults are stored in the global "ft_default" structure.
% The ft_checkconfig function that is called by many FieldTrip functions will merge
% these global configuration defaults with the cfg ctructure that you pass to
% the FieldTrip function that you are calling.
%
% The global options and their default values are
%   ft_default.colorscheme    = string, can be 'dark', 'restlees', 'bright' (default = 'dark')
%
% If you want to overrule these default settings, you can add something like this in your startup.m script
%   st_defaults
%   global st_default
%   st_default.option1 = value1
%   st_default.option2 = value2

%------- toolbox selection (FieldTrip)-------- 
%by default, Fieldtrip uses its own (slow) versions of (at least) 3 toolboxes:
%attempt to use Matlab's built-in toolboxes (Fieldtrip will check if these are available)

%first set ft_default (variable) as global
global ft_default

%set defaults for below tookboxes by adding fields to ft_default
ft_default.toolbox.signal='matlab'; %compat or matlab
ft_default.toolbox.stats='matlab';
ft_default.toolbox.images='matlab';

%ft_defaults (function) now considers ft_default (variable) and will setup paths accordingly
ft_defaults

% global ft_defaults
% ft_default.trackusage = 'no'

global st_default

if ~isfield(st_default, 'colorscheme'),       st_default.colorscheme    = 'dark';      end


try
    % this contains licensed Z3score auto scoring functions
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB', 1, 1);
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB/utilities', 1, 1);
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB/utilities/encoder', 1, 1);
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB/utilities/jsonlab', 1, 1);

    % enhanced rdir
    ft_hastoolbox('external/enhanced_rdir', 1, 1);

    % natsort
    ft_hastoolbox('external/natsortfiles', 1, 1);

    % REMs detector by Marek Adamczyk
    ft_hastoolbox('external/REMs_detector/marek_adamczyk/REMdetector', 1, 1);

end

%add cleaning functions
try; ft_hastoolbox('cleaning',1,1); end
try; ft_hastoolbox('cleaning/channels',1,1); end

% avoid using the dsp toolbox from 2021-06-03 on
% toggleToolbox('dsp','off')

%addpath(fullfile(fileparts(which('ft_defaults')), 'sleep'));
