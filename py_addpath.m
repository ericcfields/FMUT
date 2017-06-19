function new_py_path = py_addpath(directory, MATLAB_too)
%py_addpath(directory)
%
%Add directory to import search path for the instance of 
%the Python interpreter currently controlled by MATLAB
%
%   This function adds the given directory to sys.path for the current session
%   of the Python interpreter called by MATLAB. 
%   Once the directory is added, functions can be used from any module within that
%   directory just as they are from the Python standard library within MATLAB.
%   
%   Optional second argument also adds directory to MATLAB path (1 or true) 
%   or not (0 or false). The default is false. Example: 
%       py_addpath(directory,1) %this updates the MATLAB path too
%
%   Optional return value is a cell array of the directories on the updated
%   Python path. Example: 
%       new_py_path = py_addpath(directory)
%
%   To get this output without updating the Python path, use an empty string 
%   as the input:
%       py_path = py_addpath('')
%   
%   AUTHOR: Eric Fields
%   VERSION: 1.1.0
%   VERSION DATE: 27 April 2017
    
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
