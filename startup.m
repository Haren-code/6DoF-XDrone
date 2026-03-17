%% ============================================================
%  STARTUP SCRIPT FOR 6DOF PROJECT
%  Initializes paths and checks required dependencies
% =============================================================

fprintf('\n--- Initializing 6DOF Simulation Environment ---\n')

%% 1. Locate project root
project_root = fileparts(mfilename('fullpath'));

%% 2. Add all subfolders to MATLAB path
addpath(genpath(project_root));

%% 3. Remove unwanted folders from path (optional but recommended)
exclude_folders = {
    '.git'
    'data'
    'results'
    'docs'
};

for i = 1:length(exclude_folders)
    
    folder = fullfile(project_root, exclude_folders{i});
    
    if isfolder(folder)
        rmpath(genpath(folder));
    end
    
end

fprintf('Project folders added to path.\n')

function project_startup()

persistent initialized
if ~isempty(initialized)
    return
end
initialized = true;

fprintf('\n--- Initializing 6DOF Project ---\n')

%% Add project folders
project_root = fileparts(mfilename('fullpath'));
addpath(genpath(project_root));

fprintf('Project paths added.\n')

%% Required toolboxes
required_toolboxes = {
    'Aerospace_Blockset'
    'Aerospace_Toolbox'
    'Control_Toolbox'
    'Curve_Fitting_Toolbox'
    'Robust_Control_Toolbox'
    'Simulink'
    'Statistics_Toolbox'
};

missing = {};

for i = 1:length(required_toolboxes)

    if ~license('test', required_toolboxes{i})
        missing{end+1} = required_toolboxes{i}; %#ok<AGROW>
    end

end

if ~isempty(missing)

    fprintf('\nMissing required toolboxes:\n')
    
    for i = 1:length(missing)
        fprintf(' - %s\n', missing{i})
    end

    error('Install the missing toolboxes before running the simulation.')

else
    fprintf('All required toolboxes available.\n')
end

fprintf('Initialization complete.\n\n')

end

%% 5. MATLAB version check (optional)

min_version = '9.10'; % example: R2021a

if verLessThan('matlab', min_version)
    
    error(['MATLAB version too old. Minimum required: ', min_version])
    
end

fprintf('MATLAB version OK.\n')

%% 6. Clean workspace defaults (optional)

format compact
format short g

fprintf('Environment ready.\n')
fprintf('------------------------------------------------\n\n')