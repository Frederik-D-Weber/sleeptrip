function st_defaults

ft_defaults

% global ft_defaults
% ft_default.trackusage = 'no'

try
    % this contains licensed Z3score auto scoring functions
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB', 1, 1);
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB/utilities', 1, 1);
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB/utilities/encoder', 1, 1);
    ft_hastoolbox('external/automatic_sleep_scoring/z3score/z3score-api/cfslib-MATLAB/utilities/jsonlab', 1, 1);
    
    % enhanced rdir
    ft_hastoolbox('external/enhanced_rdir', 1, 1);

    
    % REMs detector by Marek Adamczyk
    ft_hastoolbox('external/REMs_detector/marek_adamczyk/REMdetector', 1, 1);

end

%addpath(fullfile(fileparts(which('ft_defaults')), 'sleep'));
