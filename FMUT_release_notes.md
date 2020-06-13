# 0.5.1

## Installation instructions

1. Download FMUT_0.5.1.zip.
2. Unzip it.
3. Put the folder of files someplace sensible (e.g., in your "MATLAB" directory).
4. [Add the folder](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) to the set of "paths" MATLAB searches when trying to answer a function call.
5. For Mac and Linux users: Run `>> add_poi_path` at the MATLAB command line, then close and restart MATLAB.

The FMUT documentation can be found at: https://github.com/ericcfields/FMUT/wiki

## Release Notes

### Bug fixes
* Error is now produced for unequal sample sizes with split plot designs.

# 0.5.0

## Installation instructions

1. Download FMUT_0.5.0.zip.
2. Unzip it.
3. Put the folder of files someplace sensible (e.g., in your "MATLAB" directory).
4. [Add the folder](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) to the set of "paths" MATLAB searches when trying to answer a function call.
5. For Mac and Linux users: Run `>> add_poi_path` at the MATLAB command line, then close and restart MATLAB.

The FMUT documentation can be found at: https://github.com/ericcfields/FMUT/wiki

## Release Notes

### New features
* `FfdrGND` (and relevant sub-functions) can now include a sphericity correction as long as there are no more than two factors with more than two levels.

### Bug fixes
* `get_mean_amplitude` now removes commas from group and bin names in csv output to avoid formatting issues
* Updated Python code for spreadsheet formatting so that highlighting will reflect the alpha level of the test (rather than always assuming alpha=0.05). Due to problems with pyinstaller, the frozen versions have not been updated.
* fixed bug in `fmut.py` due to changes in openpyxl's cell.column format 
* Changed deprecated remove_sheets method in fmut.py
* An error is now generated when a chan_hood matrix is not symmetrical

### Refactoring
(These changes are unlikely to affect most users, but could affect your code if you have called some sub-functions directly.)  

* Basic permutation ANOVA functions now output F\_obs explicitly rather than implicitly (as the first permutation in F\_dist)
* Code applying the Fmax and cluster corrections to ANOVA output have now been moved to separate functions: `Fclust_corr` and `Fmax_corr`. This was primarily to allow for more efficient simulation work (i.e., calculate the permutation ANOVA once and then applying each correction to the output).


# 0.4.0-beta

## Installation instructions

1. Download FMUT_0.4.0-beta.zip.
2. Unzip it.
3. Put the folder of files someplace sensible (e.g., in your "MATLAB" directory).
4. [Add the folder](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) to the set of "paths" MATLAB searches when trying to answer a function call.
5. For Mac and Linux users: Run `>> add_poi_path` at the MATLAB command line, then close and restart MATLAB.

The FMUT documentation can be found at: https://github.com/ericcfields/FMUT/wiki

## Release Notes

There are no known backwards incompatible changes from the 0.3.4 release. It is recommended that all users upgrade.

### New features
* New function `get_mean_amplitude`: Extract the mean amplitude for each condition across the time points/electrodes where a mass univariate test revealed a significant effect. This is useful for testing whether an ERP effect correlates with other measures (e.g., behavioral measures, individual differences).


### Bug fixes
* Fixed formatting on cluster_id tab of spreadsheet output to handle more than 20 clusters
* Updated Python code for spreadsheet formatting to account for deprecated openpyxl functions and and upcoming Python syntax changes for escape characters

# 0.3.4-beta

## Installation instructions

1. Download FMUT_0.3.4-beta.zip.
2. Unzip it.
3. Put the folder of files someplace sensible (e.g., in your "MATLAB" directory).
4. [Add the folder](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) to the set of "paths" MATLAB searches when trying to answer a function call.
5. For Mac and Linux users: run `>> add_poi_path` at the MATLAB command line. Close and restart MATLAB.

The FMUT documentation can be found at: https://github.com/ericcfields/FMUT/wiki

## Release Notes

### Minor changes
* Removed warning about three-factor designs.
* Specifying a factor with one level now produces an informative error message.

### Bug fixes
* Fixed a bug leading to an error message for black and white raster plots.

# 0.3.3-beta

## Installation instructions

1. Download FMUT_0.3.3-beta.zip
2. Unzip it
3. Put the folder of files someplace sensible (e.g., in your "MATLAB" directory)
4. [Add the folder](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) to the set of "paths" MATLAB searches when trying to answer a function call
5. **New to this release:** For Mac and Linux users, run `>> add_poi_path` at the MATLAB command line. Close and restart MATLAB.

The FMUT documentation can be found at: https://github.com/ericcfields/FMUT/wiki

## Release Notes

This is a bug fix release. There are no known backwards incompatible changes from the previous 0.3.x releases. It is recommended that all users upgrade.

### Bug fixes
* Time points no longer display in scientific notation in spreadsheet output.
* `ttest2xls` now reports correct number of subjects on the test_summary page for between subjects designs. Previously the number was one less than the true number of subjects.
* Spreadsheet outputting no longer skips some clusters on the cluster summary page.
* For Mac and Linux users the `Ftest2xls` and `ttest2xls` functions call MATLAB's `javaaddpath` function, which runs `clear all`. This can clear global variables and cause other problems. Two changes in this release help address this:
	* FMUT functions were updated to accommodate the behavior of `javaaddpath` without producing an error
	* A full solution is provided by the new function `add_poi_path`, which adds the POI libraries to the static path, avoiding the need for `javaaddpath`. See modification to installation instructions above.

### Minor and internal changes
* Changes to pyinstaller packed versions of fmut.py (spreadsheet formatting):
  * Now built from Python 3.6
  * Now use one-folder option rather than one-file, which runs faster

# 0.3.2-beta

This is a bug fix release. There are no known backwards incompatible changes from the previous 0.3.x releases. It is recommended that all users upgrade.

### Minor changes
* Updated function documentation
* Updated `py_addpath` function to handle relative paths better
* Can now use other corrections with FDR functions (undocumented)

### Bug fixes
* When the time window input specifies a time outside the epoch, an informative error message is now produces. Previously the analyses was run on a shortened time window.

# 0.3.1-beta

## Release Notes

This is a bug fix release. It is recommended that all users upgrade.

### Minor changes
* The beta warning displayed for the main functions has been made less severe, reflecting extensive additional testing of the software.
* `perm_rbANOVA` and `perm_spANOVA` now have input option to determine whether data is reduced before analysis. This is an internal change intended for testing purposes.

### Bug fixes
* Tests with between-subjects factors now report the correct number of subjects in spreadsheet output
* Spreadsheet formatting now removes all blank sheets (rather than just Sheet1)

# 0.3.0-beta

## Release notes

### Major changes
* FMUT now supports ANOVA designs with a between subjects factor. This is implemented via the new functions `FmaxGRP`, `FclustGRP`, and `FfdrGRP` which take a Mass Univariate Toolbox `GRP` variable as the input.

### Backwards compatibility
* Two new fields, `use_groups` and `group_n`, have been added to the `F_tests` struct. Attempts to add results to a GND with results missing this field will result in an error.
* Inputs to some subfunctions have changed: `calc_Fclust`, `calc_Fmax`, and `calc_param_ANOVA`  now take `cond_subs` as the second argument to specify the between subjects structure. `reduce_data` and `get_int_res` now require `cond_subs` and/or `dims` inputs.


# 0.2.0-beta

## Release notes

### Major changes
* The `'int_method'` input to permutation-based functions has been eliminated. It was ambiguous what this input meant for more complex designs with multiple interactions where an exact test is possible for some of the effects and not others and for effects where a combination of data reduction, restricted permutation, and/or permutation residuals is ideal. All effects are now reduced to the maximum extent possible. See the [FMUT Documentation](https://github.com/ericcfields/FMUT/wiki/Mass-univariate-statistics-and-corrections#how-fmut-calculates-effects-in-factorial-designs) for more details.
* Related changes in the code base mean that FMUT functions can now handle a wider array of designs beyond three-way designs.

### Backwards compatibility
* Function calls with the `'int_method'` input will now result in an error
* The structure of the `F_tests` field of the `GND` now has a new field, `exact_test`, which contains a Boolean for each effect indicating whether it was an exact test. This change in the structure will generate an error if you try to add new results to a `GND` with results in the old format.
* The `test_results` struct returned by several sub-functions also now has an `exact_test` field



# 0.1.1-beta

## Release notes

This is a bug fix release. It is recommended that all users upgrade.

### Bug fixes
* flt2GND can now handle spaces in AVGDUMP filename and path or flt filename
* flt2GND now saves correctly when ‘yes’ is specified for the ‘save_GND’ input
* Ftest2xls now correctly handles long effect names

# 0.1.0-beta


## Release notes

First beta release.