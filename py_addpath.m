%Add directory to import search path for the instance of 
%the Python interpreter currently controlled by MATLAB
%
%EXAMPLE USAGE
% >> py_addpath('C:\Documents\ERPResults')
%
%REQUIRED INPUTS
% directory      - Directory to add the Python import search path
% MATLAB_too     - If true (or 1), directory will also be added to the
%                  MATLAB path. {default: false}
%
%OPTIONAL OUTPUT
% new_py_path    - a cell array of the directories on the updated
%                  Python path; to get this output without updating the 
%                  Python path, use an empty string as the input:
%                  py_path = py_addpath('')
%
%VERSION DATE: 23 June 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function new_py_path = py_addpath(directory, MATLAB_too)
    
    %check input
    if ~ischar(directory)
        error('Input must be a string')
    elseif ~exist(directory, 'dir') && ~isempty(directory)
        error('%s is not a valid directory', directory)
    end
    
    %add directory to Python path if not already present
    if ~any(strcmp(get_py_path(), directory))
        py_path = py.sys.path;
        py_path.insert(int64(1), directory);
    end
    
    %add directory to MATLAB path if requested
    if nargin>1 && MATLAB_too
        addpath(directory);
    end
    
    %optionally return ammended path.sys as cell array
    if nargout
        new_py_path = get_py_path();
    end
    
end

function current_py_path = get_py_path()
%Function to return the current python search path as a cell array of strings
    current_py_path = cellfun(@char, cell(py.sys.path), 'UniformOutput', 0)';
end
