%Add useful formatting to FMUT output spreadsheets
%
%EXAMPLE USAGE
% >> format_xls('results.xlsx')
%
%REQUIRED INPUTS
% spreadsheet - name of spreadsheet file to format
%
%VERSION DATE: 29 March 2019
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function format_xls(spreadsheet, alpha)

    if nargin < 2
        alpha = 0.05;
    end

    %If full path is not provided assume spreadsheet is the current working
    %directory
    if ~isdir(fileparts(spreadsheet)) %#ok<ISDIR>
        spreadsheet = [pwd filesep spreadsheet];
    end

    %Check is spreadsheet exists and return gracefully if it doesn't
    if ~exist(spreadsheet, 'file')
        watchit(sprintf('Spreadsheet formatting error:\n%s does not exist.', spreadsheet));
        return;
    end
    
    %Find FMUT directory
    func_dir = fileparts(which('format_xls'));
    
    %Try using .py script
    try
        py_addpath(func_dir);
        py.fmut.format_xls(spreadsheet, alpha);
        
    %If .py version doesn't work, use pyinstaller version
    catch
        try
            if ispc()
                system(sprintf('"%s\\py_fmut.exe" --format_xls "%s"', fullfile(func_dir, 'py_fmut_win'), spreadsheet));
            else
                system(sprintf('"%s/py_fmut" --format_xls "%s"', fullfile(func_dir, 'py_fmut_mac'), spreadsheet));
            end    
        catch
            watchit(sprintf('Unable to format spreadsheet.'))
        end
        
    end
        
end
