%Add useful formatting to FMUT output spreadsheets
%
%EXAMPLE USAGE
% >> format_xls('results.xlsx')
%
%REQUIRED INPUTS
% spreadsheet - name of spreadsheet file to format
%
%VERSION DATE: 17 November 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function format_xls(spreadsheet)
    
    %Check is spreadsheet exists and return gracefully if it doesn't
    if ~exist(spreadsheet, 'file')
        watchit(sprintf('Spreadsheet formatting error:\n%s does not exist.', spreadsheet));
        return;
    end
    
    %Find FMUT directory
    func_dir = fileparts(which('format_xls'));
    
    %Try using .py script
    try
        py_addpath(func_dir)
        if ~isdir(fileparts(spreadsheet))
            spreadsheet = [pwd filesep spreadsheet];
        end
        py.fmut.format_xls(spreadsheet)
        
    %If .py version doesn't work, use pyinstaller version
    catch
        if ispc()
            try
                system(sprintf('"%s\\py_fmut.exe" --format_xls "%s"', fullfile(func_dir, 'py_fmut_win'), spreadsheet));
            catch
                watchit('Unable to format spreadsheet.');
            end
        else
            try
                system(sprintf('"%s/py_fmut" --format_xls "%s"', func_dir, spreadsheet));
            catch
                watchit(sprintf('Unable to format spreadsheet.\nSee FMUT documentation for issues with spreadsheet formatting on non-Windows systems and a possible workaround.'))
            end
        end
    end
    
end
