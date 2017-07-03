% flt2GND() - Create a Mass Univariate ERP Toolbox GND struct from a set of
%             .flt files created by Phil Holcomb's ERP software 
%
% EXAMPLE USAGE:
%  >> GND = flt2GND(infiles, 'bsln', -100, 'sampling_rate', 200, 'use_bins', 1:6)
%
% REQUIRED INPUTS:
%  infiles        - Either a cell array of strings specifying .flt file
%                   names or a string giving the fullpath of a text file
%                   with an .flt file name on each line. In either case,
%                   .flt file names should be *filenames only*. The filepath
%                   for their location is given via the filepath input
%  bsln           - The number of milliseconds before 0 in each epoch. This
%                   is specified in the "Presample" field of your .scp file
%  sampling_rate  - Sampling rate of the date in Hz
%  use_bins       - Which bins to include in the GND struct
%
% OPTIONAL INPUTS:
%  filepath       - Full filepath for .flt files. {default: current working
%                   directory}
%  bin_desc       - cell array with condition names for bins 
%                   {default: no descriptions}
%  avgdump_exe    - filepath and filename for AVGDUMP.EXE
%                   {default: 'C:\BIN\AVGDUMP.EXE'}
%  chanlocs_file  - full path for chanlocs file
%                   {default: no chanlocs information}
%  save_GND       - 'yes' or 'no' or filename. If 'yes', a gui will appear prompting
%                   the user to give a name and location to save the GND
%                   variable. If a filename, it will save automatically.
%                   {default: 'yes'}
%  plot_gui       - 'yes' or 'no'. If 'yes', a gui for exploring the data
%                   in the newly created GND variable will pop up when the 
%                   function finishes. {default: 'yes'}
%  n_electrodes   - number of electrodes {default: 32}
%  n_samples      - number of sample points in each epoch {default: 256}
%  exp_name       - Name of the study to add to the GND struct. 
%                   {default: ''}
%  verblevel      - An integer specifiying the amount of information you want
%                   this function to provide about what it is doing during runtime.
%                    Options are:
%                      0 - quiet, only show errors, warnings, and EEGLAB reports
%                      1 - stuff anyone should probably know
%                      2 - stuff you should know the first time you start working
%                          with a data set {default value}
%                      3 - stuff that might help you debug (show all
%                          reports)
%
% OUTPUT:
%  GND            - Mass Univariate Toolbox GND struct
%
% Notes:
% -The GND fields odelay, cals, and condesc are specific to Kutaslab data.
% Other labs should be able to ignore them.
%
% AUTHOR: Eric Fields
% VERSION DATE: 3 July 2017
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed. 
%This is a beta version of this software. It needs additional testing 
%and should not be considered error free. 

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%% REVISION HISTORY %%%%%%%%%%%%%%%%%%%%%%
%4/6/17   - First version
%4/7/17   - Several minor updates and improvements
%4/8/17   - Added ability to specify different number of channels and
%           samples; some code refactoring; a few errors corrected
%4/18/17  - Added checks and informative messages for incorrect filepath and
%           incorrect number of sampling points or electrodes; added check
%           and error message for non-existent .flt file or incorrect
%           filepath; fixed bug where gui_erp would cause error with no
%           chanlocs file
%6/1/17   - Fixed bug when using bins that don't start at 1
%6/13/17  - GND.indiv_bin_ct now filled with -1
%6/14/17  - bsln_wind field now spelled correctly
%6/27/17  - Can now auto-save
%7/3/17   - Filnames can now have spaces

function GND = flt2GND(infiles, varargin)

    warning('You are using a beta version of flt2GND. It needs further testing and should not be considered error free.');
    
    global VERBLEVEL
    
    %flt2GND relies on the Windows program AVGDUMP.EXE
    if ~ispc()
        error('Because it relies on AVGDUMP.EXE, flt2GND only works on Windows computers.')
    end
    
    %% Parse input
    
    %If infiles is a text file, import information to cell array
    if ischar(infiles) && exist(infiles, 'file')==2
        %Get cell array of flt files from text file
        f_in = fopen(infiles);
        subs = textscan(f_in, '%s', 1000);
        subs = subs{1};
        fclose(f_in);
    elseif ~iscell(infiles)
        error('The infiles input does not appear to be valid.');
    end

    %Assign keyword-value pair inputs to variables
    for i = 1:length(varargin)
        input = varargin{i};
        if ischar(input)
            switch input
                case 'filepath'
                    filepath = varargin{i+1};
                case 'bsln'
                    bsln = varargin{i+1};
                case 'sampling_rate'
                    srate = varargin{i+1};
                case 'use_bins'
                    use_bins = varargin{i+1};
                case 'bin_desc'
                    bin_desc = varargin{i+1};
                case 'avgdump_exe'
                    avgdump_exe = varargin{i+1};
                case 'exp_name'
                    exp_name = varargin{i+1};
                case 'chanlocs_file'
                    chanlocs_file = varargin{i+1};
                case 'save_GND'
                    save_GND = varargin{i+1};
                case 'plot_gui'
                    plot_gui = varargin{i+1};
                case 'n_electrodes'
                    n_electrodes = varargin{i+1};
                case 'n_samples'
                    n_samples = varargin{i+1};
                case 'verblevel'
                    VERBLEVEL = varargin{i+1};
            end
        end
    end
    
    %Check for required inputs
    if ~exist('bsln', 'var')
        error('''bsln'' input is required')
    end
    if ~exist('srate', 'var')
        error('''sampling_rate'' input is required')
    end
    if ~exist('use_bins', 'var')
        error('''use_bins'' input is required')
    end
    
    %Set defaults for missing inputs
    if ~exist('filepath', 'var')
        filepath = pwd();
    end
    if ~exist('avgdump_exe', 'var')
        avgdump_exe = 'C:\BIN\AVGDUMP.EXE';
    end
    if ~exist('exp_name', 'var')
        exp_name = '';
    end
    if ~exist('chanlocs_file', 'var')
        chanlocs_file = false;
    end
    if ~exist('save_GND', 'var')
        save_GND = 'yes';
    end
    if ~exist('plot_gui', 'var')
        plot_gui = 'yes';
    end
    if ~exist('n_electrodes', 'var')
        n_electrodes = 32;
    end
    if ~exist('n_samples', 'var')
        n_samples = 256;
    end
    if ~exist('bin_desc', 'var')
        bin_desc = repmat({''}, 1, length(use_bins));
    end
    if isempty(VERBLEVEL)
        VERBLEVEL = 2;
    end
    
    %Formatting and errors
    if length(use_bins) ~= length(bin_desc)
        error('When using the ''bin_desc'' option, you must provide the same number of descriptions as bins you are importing.')
    end
    if ~exist(avgdump_exe, 'file')
        error('%s does not exist. Please specify a correct filepath and filename for AVGDUMP', avgdump_exe)
    end
    if chanlocs_file
        if ~exist('readlocs', 'file')
            error('To read chanlocs information, you must have EEGLAB installed and running.')
        end
        if ~exist(chanlocs_file, 'file')
            error('%s does not exist. Please specify a valid filepath and filename for the chanlocs file.')
        end
    else
        watchit('No chanlocs file specified. Tests that rely on channel location will not be possible.')
        if ~strcmpi('plot_gui', 'no') && ~strcmpi('plot_gui', 'n')
            watchit('Visualization of GND requires a chanlocs file. No visualization will be created.')
            plot_gui = 'no';
        end
    end
    if ~strcmpi(plot_gui, 'no') && ~strcmpi(plot_gui, 'n') && ~exist('icadefs', 'file')
        watchit('Visualization of GND requires EEGLAB to be installed and running. No visualization will be created.')
        plot_gui = 'no';
    end
    for i = 1:length(subs)
        if ~exist(fullfile(filepath, subs{i}), 'file')
            error('%s does not exist. Please check the ''infiles'' or ''filepath'' input.', fullfile(filepath, subs{i}))
        end
    end
    if ~strcmpi(save_GND, 'yes') && ~strcmpi(save_GND, 'y')&& ~strcmpi(save_GND, 'no') && ~strcmpi(save_GND, 'n') && ~strcmpi(save_GND(end-3:end), '.GND')
        error('''save_GND'' input must be ''yes'', ''no'', or a valid GND filename')
    end

    
    %% Make GND
    
    %Some useful numbers
    n_subs = length(subs);
    n_bins = length(use_bins);

    %Make GND struct
    GND = struct;
    GND.exp_desc         = exp_name;
    GND.filename         = '';
    GND.filepath         = '';
    GND.saved            = 'no';
    GND.grands           = NaN(n_electrodes, n_samples, n_bins);
    GND.grands_stder     = NaN(n_electrodes, n_samples, n_bins);
    GND.grands_t         = NaN(n_electrodes, n_samples, n_bins);
    GND.sub_ct           = ones(1, n_bins) * n_subs;
    GND.chanlocs         = NaN;
    GND.bin_info         = struct('bindesc', bin_desc, 'condcode', repmat({1}, 1, n_bins));
    GND.condesc          = {'Experiment (not cal pulses)'};
    GND.time_pts         = NaN(1, n_samples);
    GND.bsln_wind        = [-abs(bsln) 1000/srate];
    GND.odelay           = [];
    GND.srate            = srate;
    GND.indiv_fnames     = cell(1, n_subs);
    GND.indiv_subnames   = cell(1, n_subs);
    GND.indiv_traits     = [];
    GND.indiv_bin_ct     = -ones(n_subs, n_bins);
    GND.indiv_bin_raw_ct = NaN;
    GND.indiv_erps       = NaN(n_electrodes, n_samples, n_bins, n_subs);
    GND.indiv_art_ics    = num2cell(NaN(1, n_subs));
    GND.cals             = [];
    GND.history          = {};
    GND.t_tests          = [];

    %Calculate times (in ms) corresponding to each sampling point for
    %time_pts field of GND
    incr = 1000/srate;
    end_point = (n_samples - abs(bsln)/incr - 1) * incr;
    GND.time_pts = -abs(bsln):incr:end_point;

    %Get ERP data from .flt files for indiv_erps field of GND
    old_dir = cd(filepath);
    for s = 1:n_subs
        flt_file = subs{s};
        if VERBLEVEL
            fprintf('Importing data from %s\n', flt_file);
        end
        GND.indiv_fnames{s}   = fullfile(filepath, flt_file);
        GND.indiv_subnames{s} = flt_file(1:end-4);
        for b = 1:n_bins
            flt_bin = use_bins(b);
            try
                [~, avgdump_output] = system(sprintf('"%s" "%s" %d -t -fm %d', avgdump_exe, flt_file, flt_bin, n_electrodes));
                GND.indiv_erps(:, :, b, s) = str2num(avgdump_output)'; %#ok<ST2NM>
            catch ME
            %If the data returned by AVGDUMP is not the expected size, give
            %an informative error message:
                    if strcmpi(ME.identifier, 'MATLAB:subsassigndimmismatch')
                        if size(str2num(avgdump_output),1) ~= n_samples %#ok<ST2NM>
                            error('Your data has %d sampling points, not %d. Please revise ''n_samples'' input.', ...
                                  size(str2num(avgdump_output),1), n_samples);  %#ok<ST2NM>
                        elseif size(str2num(avgdump_output),2) ~= n_elecrodes %#ok<ST2NM>
                            error('Your data has %d electrodes, not %d. Please revise ''n_electrodes'' input.', ...
                                  size(str2num(avgdump_output),2), n_electrodes);  %#ok<ST2NM>
                        else
                            rethrow(ME)
                        end
                    else
                        rethrow(ME)
                    end
            end
        end
    end
    cd(old_dir);
    GND.indiv_erps = GND.indiv_erps/100; %AVGDUMP returns data in hundreds of microvolts
    if abs(mean(mean(mean(mean(GND.indiv_erps(:, 1:sum(GND.time_pts<0), :, :)))))) > .01
        watchit('The average of the pre-zero period is greater than .01 microvolts. Check that data is baselined correctly.')
    end

    %Calculate grands, grands_stder, and grands_t from ERP data
    GND.grands       = mean(GND.indiv_erps, 4);
    GND.grands_stder = std(GND.indiv_erps, 0, 4) / sqrt(n_subs);
    GND.grands_t     = GND.grands ./ GND.grands_stder;

    %Add chanlocs information
    if chanlocs_file
        GND.chanlocs = readlocs(chanlocs_file);
        GND.chanlocs = rmfield(GND.chanlocs, 'sph_theta_besa');
        GND.chanlocs = rmfield(GND.chanlocs, 'sph_phi_besa');
        [GND.chanlocs.ref] = deal('');
        [GND.chanlocs.type] = deal('');
        [GND.chanlocs.urchan] = deal([]);
        [GND.chanlocs.sph_radius] = deal(NaN);
        GND.chanlocs = orderfields(GND.chanlocs, {'Y', 'X', 'Z', 'labels', 'sph_theta', 'sph_phi', 'sph_radius', 'theta', 'radius', 'ref', 'type', 'urchan'});
    end
    
    %Save GND
    if ~strcmpi(save_GND, 'no') && ~strcmpi(save_GND, 'n')
        if strcmpi(save_GND(end-3:end), '.GND')
            if isempty(fileparts(save_GND))
                GNDpath = pwd;
                GNDname = save_GND;
            else
                GNDpath = fileparts(save_GND);
                [~, GNDname] = fileparts(save_GND);
                GNDname = [GNDname '.GND'];
            end
            GND = save_matmk(GND, GNDname, GNDpath, 1);
        else
            GND = save_matmk(GND, 'gui');
        end
    end
    
    %Visually inspect GND
    if ~strcmpi(plot_gui, 'no') && ~strcmpi(plot_gui, 'n')
        gui_erp(GND);
    end
    
end
