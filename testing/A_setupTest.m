%FMUT unit testing
%General setup
%
%AUTHOR: Eric Fields
%VERSION DATE: 12 June 2020

clear all; close all; path(pathdef); %#ok<CLALL>

%Add EEGLAB and MUT to path
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%Remove FMUT directory within EEGLAB plugins folder
rmpath(fileparts(which('FclustGND')));

%Add FMUT folder to search path
FMUT_dir = fileparts(fileparts(which('A_setupTest')));
addpath(FMUT_dir);

%Clear everything before testing
close all;
clear all;  %#ok<CLALL>

global test_xls_output;
test_xls_output = true;

%Clear spreadsheet outputs folder
if test_xls_output
    FMUT_dir = fileparts(fileparts(mfilename('fullpath')));
    cd(fullfile(FMUT_dir, 'testing', 'outputs'));
    delete *.xlsx
    cd(fullfile(FMUT_dir, 'testing'));
end

%Make sure we are using the right FMUT functions
assert(strcmp(fileparts(which('FclustGND')), FMUT_dir));
