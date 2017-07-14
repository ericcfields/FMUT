%Run all test outside of MATLAB's unit testing framework
%
%AUTHOR: Eric Fields
%VERSION DATE: 14 July 2017

global test_xls_output
user_resp = input('Test spreadsheet output?(y/n) ', 's');
if strcmpi(user_resp, 'n')
    test_xls_output = false;
else
    test_xls_output = true;
end

tic
A_setupTest
B_FmaxGNDTest
C_FclustGNDTest
D_FfdrGNDTest
E_VariousTest
F_FmaxGRPTest

fprintf('\n\n')
fprintf('All tests done.!\n')
toc
fprintf('\n\n')