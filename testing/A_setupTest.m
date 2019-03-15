%FMUT unit testing
%General setup
%AUTHOR: Eric Fields
%VERSION DATE: 15 March 2019

%Add FMUT folder to search path
FMUT_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(FMUT_dir);
addpath('C:\Users\ecfne\Documents\MATLAB\dmgroppe-Mass_Univariate_ERP_Toolbox-d1e60d4');

%Add EEGLAB and MUT to path
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

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
