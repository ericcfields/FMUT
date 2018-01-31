%Add Apache POI library to MATLAB's static Java class path
%
%EXAMPLE USAGE
% >> add_poi_path
%
%INPUT
% segment - if '-static', add POI library .jar files to the static path; if
%           '-dynamic', add to the dynamic path. {default: '-static'}
%
%VERSION DATE: 20 November 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function add_poi_path(segment)
    
    %Default to updating static path
    if ~nargin
        segment = '-static';
    end
    
    %Get full path for POI .jar files not on static Java class path
    [~, missing_poi_files] = get_poi_paths();
    
    %Nothing to do if all files are on static path
    if isempty(missing_poi_files)
        return;
    end
    
    if strcmpi(segment, '-static')

        %Full path for javaclasspath file
        static_file = [prefdir filesep 'javaclasspath.txt'];
        %Append missing POI files to javaclasspath.txt file
        f_out = fopen(static_file, 'a');
        for i = 1:length(missing_poi_files)        
            fprintf(f_out, '%s\n', missing_poi_files{i});
        end
        fclose(f_out);
    
        %POI files will not actually be added to the Java class path until
        %MATLAB is restarted
        msgbox('Close and restart MATLAB to finish adding the Apache POI library to the Java class path.')
        
    elseif strcmpi(segment, '-dynamic')
        
        %Add any file no on the static or dynamic path to the dynamic path
        for i = 1:length(missing_poi_files)
            jar_file = missing_poi_files{i};
            if ~ismember(jar_file, javaclasspath('-dynamic'))
                javaaddpath(jar_file);
            end
        end
        
    end
    
end
    