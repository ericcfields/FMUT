%Return the full paths of the .jar files of the Apache POI library
%
%EXAMPLE USAGE
% >> get_poi_paths
%
%OUTPUT
% poi_files          - Full paths for all .jar files in the Apache POI
%                      library
% missing_poi_files  - Full paths for .jar files in the Apache POI library
%                      that are not on the static Java class path
%
%VERSION DATE: 20 November 2017
%AUTHOR: Eric Fields
%
%NOTE: This function is provided "as is" and any express or implied warranties 
%are disclaimed.

%Copyright (c) 2017, Eric Fields
%All rights reserved.
%This code is free and open source software made available under the 3-clause BSD license.

function [poi_files, missing_poi_files] = get_poi_paths()
    
    %Find directory with POI library
	poi_dir = fullfile(fileparts(which('get_poi_paths')), 'poi_library');
    
    %Find all .jar files within directory
    poi_files = dir(poi_dir);
    poi_files = {poi_files.name};
    poi_files = poi_files(cell2mat(cellfun(@(x) ~isempty(strfind(x, '.jar')), poi_files, 'UniformOutput', false))); %#ok<STREMP>
    poi_files = cellfun(@(x) [poi_dir filesep x], poi_files, 'UniformOutput', false);
    
    %Find which .jar files (if any) are not on the static Java class path
    missing_poi_files = poi_files(~ismember(poi_files, javaclasspath('-static')));

end