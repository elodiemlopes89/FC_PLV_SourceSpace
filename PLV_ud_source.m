function out_seg = PLV_ud_source(PatID, N, Stim, FR)
%
% PLV_ud_source - This function computes the phase locking value (PLV) 
% between EEG sources in different brain regions for a given patient.
% It performs preprocessing, segmentation, source analysis, and computes 
% undirected (PLVu) and directed (PLVd) network measures. Additionally, 
% it outputs global connectivity metrics for the given EEG segments.
%
% Estimation of Brain Functional Connectivity (FC) Networks in the source space
% Fieldtrip is used to construct the Head Model
%
% Headmodel: brainstorm_ICBM152_mni
%
% Elodie M. Lopes (elodie.m.lopes@inesctec.pt)
% Doctoral Program in Biomedical Engineering (FEUP)
% Supervisor: João P. Cunha (INESC TEC, Porto, Portugal)
% 2024
%
% Syntax:
%   out_seg = PLV_ud_source(PatID, N, Stim, FR)
%
% Inputs:
%   - PatID: A string representing the patient ID (e.g., '01').
%   - N: A string or number representing the segment number (e.g., '1').
%   - Stim: A character ('a' or 'b') indicating whether the data is 
%     from after ('a') or before ('b') stimulation.
%   - FR: A vector of two values, [fmin, fmax], specifying the frequency range 
%     of interest for filtering (e.g., [1, 30] Hz).
%
% Outputs:
%   - out_seg: A structure containing the following fields:
%       - eeg: EEG data for the segment.
%       - sources: The source-reconstructed EEG signals.
%       - net_u: Undirected network connectivity matrix.
%       - net_d: Directed network connectivity matrix.
%       - GC: Global connectivity metrics (Clustering coefficient, Node degree, etc.)
%       - ROIs: Labels for the regions of interest.
%

%% (I) Add required MATLAB packages
% Add paths to necessary packages (ensure they are available for your system).
% PACKAGES:
% Fieldtrip - 20191119 (@Free Software Foundation, Inc, 1991)
% BrainNet Viewer (https://www.mathworks.com/matlabcentral/fileexchange/68881-brainnet-viewer)
% Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/)
% EEG Preprocessing custom-made package (https://github.com/elodiemlopes89/Electrophysiological-Data-Preprocessing-Visualization)

% Directory setup
dir = pwd;  % Get the current working directory.
dir_Preprocessing = [dir, '/Packages/Preprocessing_EEG'];
dir_FT = [dir, '/Packages/fieldtrip-20191119'];
dir_FT_fileio = [dir, '/Packages/fieldtrip-20191119/fileio'];
dir_FT_utilities = [dir, '/Packages/fieldtrip-20191119/utilities'];
dir_BCT = [dir, '/Packages/2019_03_03_BCT'];
dir_BNV = [dir, '/Packages/BrainNetViewer_20191031'];

% Add paths to each toolbox
addpath(dir_Preprocessing); 
addpath(dir_FT); 
addpath(dir_FT_fileio); 
addpath(dir_FT_utilities); 
addpath(dir_BCT); 
addpath(dir_BNV);

%% (II) Define data directories
dir_data = [dir, '/Data'];  % Directory containing EEG data files.

%% (2) Preprocessing EEG data

% Load EEG data for the specified patient and segment
filename = ['PD', PatID, Stim, '_seg', N, '.mat'];
load([dir_data, '/', filename]);  % Load the data file
cd(dir);  % Change back to the current directory

% Set up output segment information
out_seg.info = ['PD', PatID, Stim, '_seg', N];

% Initial variables from loaded data
labels0 = channels;  % Channel labels
freq = sr;  % Sampling rate
eeg = data;  % EEG data matrix
clear data;  % Clear original data to save memory

% Remove irrelevant channels
labels2remove = {'Event'; 'F11'; 'F12'; 'P11'; 'P12'};  % Channels to remove
AM = 1;  % Average montage flag (set to 1)
fs = 256;  % Sampling frequency (256 Hz)

% Set frequency range for filtering
f_min = FR(1);
f_max = FR(2);

% Apply preprocessing to EEG data
[new_data, labels] = preprocessing(eeg, labels0, freq, labels2remove, AM, fs, f_min, f_max);

%% (3) Electrode Info (Specific to the Patient)
% Load electrode information (standard 1005 electrode layout)
grad = ft_read_sens('standard_1005.elc', 'senstype', 'eeg');

% Define electrode positions for the patient
elec.label = labels';
indice1 = [4, 19, 41, 63, 84, 21, 43, 65, 29, 34, 45, 6, 27, 49, 71, 86, 25, 47, 69, 39, 12, 67, 8, 16, 30, 38];
A = grad.chanpos(indice1, :);  % Specific electrode positions for this patient

% Define additional reference positions for F11, F12, P11, P12
r = 87.54;  % Radius for positioning
a = pi / 180;  % Degree to radians conversion
F11 = [r * cos(-130 * a) * cos(-40 * a), r * sin(-130 * a) * cos(-40 * a), r * cos(-130 * a)];
F12 = [r * cos(130 * a) * cos(40 * a), r * sin(130 * a) * cos(40 * a), r * cos(130 * a)];
P11 = [r * cos(-130 * a) * cos(40 * a), r * sin(-130 * a) * cos(40 * a), r * cos(-130 * a)];
P12 = [r * cos(130 * a) * cos(-40 * a), r * sin(130 * a) * cos(-40 * a), r * cos(130 * a)];

% Define positions for B and C regions (electrodes)
indice2 = [51, 50];  % Electrode indices for B
B = grad.chanpos(indice2, :);
indice3 = [74, 82, 31, 37];  % Electrode indices for C
C = grad.chanpos(indice3, :);

% Combine all electrode positions into a matrix M
M = [A; B; C]; 

% Assign electrode positions to output structure
elec.chanops = M;
elec.elecpos = M;
elec.chantype = repmat({'eeg'}, size(M, 1), 1);
elec.chanunit = repmat({'V'}, size(M, 1), 1);
elec.unit = 'mm';

out_seg.channels = labels';  % Store channel labels

%% (4) Segmentation of EEG Data

% Split data into segments of 20 seconds
N_samples = length(new_data(1, :));  % Total number of samples
N_20s = 5120;  % 20 seconds worth of data (256 Hz * 20)
N_seg = 4;  % Number of segments to analyze

% Divide the data into segments of 20 seconds each
seg{1, 1} = new_data(:, 1:5120);
for i = 2:N_seg
    seg{1, i} = new_data(:, (i - 1) * N_20s + 1 : i * N_20s);
end

% Initialize containers for storing segment-related data
seg_pat = cell(1, N_seg);
source_pat = cell(1, N_seg);
net_u_pat = cell(1, N_seg);
net_d_pat = cell(1, N_seg);
GC_pat = cell(1, N_seg);

%% (5) Main Processing Loop

fprintf('Beginning...'), clock = tic;  % Start timing the processing

% Loop through each segment and perform source analysis
for i_seg = 1:N_seg
    % A. Convert the segment to FieldTrip structure
    elecLabel = labels;
    fsample = fs;
    data = data2fieldtrip(seg{1, i_seg}, elecLabel, fsample, 20); %data2fieldtrip.m function
    
    seg_pat{i_seg} = data.trial{1, 1};  % Store segment data
    
    % B. Prepare Head Model (Template)
    load('brainstorm_ICBM152_mni', 'headmodel');  % Load template head model
    mesh = struct;  % Create mesh structure for BEM model
    shells = {'head', 'skull', 'brain'};  % List of BEM model shells

    % Loop through shells and prepare mesh for BEM
    for i = 1:length(shells)
        mesh(i).pos = headmodel.(shells{i}).pos;
        mesh(i).tri = headmodel.(shells{i}).tri;
        mesh(i).unit = 'mm';
        mesh(i).coordsys = 'mni';
    end

    % Set up configuration for BEM model
    cfg = struct;
    cfg.tissue = {'scalp', 'skull', 'brain'};
    cfg.method = 'bemcp';
    headmodel_ft = ft_prepare_headmodel(cfg, mesh);

    % C. Prepare Electrode Model
    cfg = struct;
    cfg.channel = labels';
    cfg.method = 'project';
    cfg.headshape = headmodel_ft.bnd(1);
    elec2 = ft_electroderealign(cfg, elec);

    % D. Prepare Source Grid
    grid = struct;
    grid.pos = headmodel.cortex.pos;
    grid.inside = true(length(grid.pos), 1);

    % E. Prepare Leadfield
    cfg = struct;  % Initialize configuration
    cfg.headmodel = headmodel_ft;  % BEM head
