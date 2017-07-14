%FMUT unit testing
%General setup
%AUTHOR: Eric Fields
%VERSION DATE: 14 July 2017

global test_xls_output;
test_xls_output = true;

if ispc()

    addpath('R:\Public\GK_lab\Eric\FMUT_development\FMUT');
    %addpath('R:\Public\GK_lab\Eric\FMUT_development\FMUT\dev');

    %Clear spreadsheet outputs folder
    if test_xls_output
        cd('R:\Public\GK_lab\Eric\FMUT_development\FMUT\testing\outputs')
        delete *.xlsx
        cd('R:\Public\GK_lab\Eric\FMUT_development\FMUT\testing\')
    end

    %Add EEGLAB and MUT to path
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    %Clear everything before testing
    clear all;  %#ok<CLSCR>
    close all;
    
elseif ismac()
    
    addpath('/Users/efield02/Desktop/ncl_eeglab-erplab-pipeline');
    addpath('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT');
    %addpath('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT/dev');

    %Clear spreadsheet outputs folder
    if test_xls_output
        cd('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT/testing/outputs')
        delete *.xlsx
        cd('/Volumes/as-rsch-ncl1$/Public/GK_lab/Eric/FMUT_development/FMUT/testing/')
    end

    %Add EEGLAB and MUT to path
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

    %Clear everything before testing
    close all;
    clearvars;
    
end
