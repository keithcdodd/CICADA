function test_option1_sha256_fixture
% Read-only deterministic SHA-256 test for the companion helper.
root=fileparts(fileparts(mfilename('fullpath')));
helpers=fullfile(root,'helper_functions');
p=path;
restore=onCleanup(@()path(p)); %#ok<NASGU>
addpath(helpers,'-begin');
fixture=fullfile(root,'tests','sha256_fixture.txt');
actual=cicada_task_glm_sha256(fixture);
expected='2d412275ac8c0ce1e624db06fd3cbbcbd13c400f48aae3900ffd7238aed674ce';
assert(strcmp(actual,expected),'SHA256 fixture mismatch');
fprintf('SHA256 fixture check PASSED: %s\n',actual);
end
