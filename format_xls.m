%Add useful formatting to FMUT output spreadsheets
%
%REQUIRED INPUTS
% spreadsheet - name of spreadsheet file to format
%
%VERSION DATE: 16 June 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

%%%%%%%%%%%%%%%%%%%  REVISION LOG   %%%%%%%%%%%%%%%%%%%
% 5/17/17   - First version
% 5/23/17   - Updated to use compiled version of formatting code
% 6/14/17   - Fixed problem with spreadsheet names including spaces
% 6/16/17   - Updated for compatability with Mac

function format_xls(spreadsheet)
    func_dir = fileparts(which('format_xls'));
    try
        py_addpath(func_dir)
        if ~isdir(fileparts(spreadsheet))
            spreadsheet = [pwd filesep spreadsheet];
        end
        py.fmut.format_xls(spreadsheet)
    catch
        if ispc()
            system(sprintf('"%s\\py_fmut.exe" --format_xls "%s"', func_dir, spreadsheet));
        else
            try
                system(sprintf('"%s/py_fmut" --format_xls "%s"', func_dir, spreadsheet));
            catch
                watchit(sprintf('Unable to format spreadsheet.\nSee FMUT documentation for issues with spreadsheet formatting on non-Windows systems and a possible workaround.'))
            end
        end
    end
end
