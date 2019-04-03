%FMUT unit testing
%General setup
%AUTHOR: Eric Fields
%VERSION DATE: 3 April 2019

%Add EEGLAB and MUT to path
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%Remove FMUT directory within EEGLAB plugins folder
rmpath(fullfile(fileparts(which('eeglab')), 'plugins', 'FMUT_0.4.1'));

%Add FMUT folder to search path
FMUT_dir = fileparts(fileparts(mfilename('fullpath')));
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
