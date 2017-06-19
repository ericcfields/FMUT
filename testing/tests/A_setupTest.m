%FMUT unit testing
%General setup
%AUTHOR: Eric Fields
%VERSION DATE: 15 June 2017

if ispc()

    addpath('R:\Public\GK_lab\Eric\FMUT_development\FMUT_functions');
    addpath('R:\Public\GK_lab\Eric\FMUT_development\FMUT_functions\in_progress');

    %Clear spreadsheet outputs folder
    cd('R:\Public\GK_lab\Eric\FMUT_development\testing\tests\outputs')
    delete *.xlsx
    cd('R:\Public\GK_lab\Eric\FMUT_development\testing\tests')

    %Add EEGLAB and MUT to path
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    %Clear everything before testing
    clear all;  %#ok<CLSCR>
    close all;
    
elseif ismac()
    
    addpath('/Users/efield02/Desktop/ncl_eeglab-erplab-pipeline');
    addpath('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT_functions');
    addpath('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT_functions/in_progress');

    %Clear spreadsheet outputs folder
    cd('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/testing/tests/outputs')
    delete *.xlsx
    cd('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/testing/tests')

    %Add EEGLAB and MUT to path
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    %Clear everything before testing
    close all;
    clearvars;
    
end
