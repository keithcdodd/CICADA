function test_reml_cached_vs_reference()
% Self-contained synthetic numerical regression test. No image files,
% CICADA runs, FSL calls, or output writes. MATLAB runtime required.
folder=fileparts(mfilename('fullpath'));
root=fileparts(folder);
oldPath=path;
oldRng=rng;
restorePath=onCleanup(@()path(oldPath)); %#ok<NASGU>
restoreRng=onCleanup(@()rng(oldRng)); %#ok<NASGU>
addpath(root);
addpath(fullfile(root,'helper_functions'),'-begin');
addpath(fullfile(folder,'reference'),'-begin');
rng(20261101,'twister');

N=72;
t=(0:N-1)';
X0=[ones(N,1), sin(t/11), cos(t/13), ...
    (t-mean(t))/N, randn(N,1)];
C0=[0 1 -1 0 0; 0 0 1 -1 0];
Y=X0*[1000 1001 998 1002; .8 0 -.3 .4; ...
    .1 .4 .2 -.7; .2 .0 .1 -.2; .4 -.1 .6 .2] + randn(N,4);

for model={'AR1','AR2'}
    opts=struct('model',model{1},'maxBatchVoxels',2);
    cache=cicada_task_glm_reml_prepare(X0,C0,opts);
    p1=cicada_task_glm_reml_fit_batch(Y(:,1:2),cache);
    p2=cicada_task_glm_reml_fit_batch(Y(:,3:4),cache);
    ref=cicada_task_glm_reml_core(Y,X0,C0, ...
        struct('model',model{1},'maxVoxels',4));
    score=[p1.remlProfileScore p2.remlProfileScore];
    cope=[p1.cope p2.cope];
    vcope=[p1.varcopePlugIn p2.varcopePlugIn];
    r1=[p1.r1Hat p2.r1Hat];
    k2=[p1.kappa2Hat p2.kappa2Hat];
    assert(max(abs(score-ref.remlProfileScore))<1e-8);
    assert(max(abs(cope(:)-ref.cope(:)))<1e-8);
    assert(max(abs(vcope(:)-ref.varcopePlugIn(:)) ...
        ./ref.varcopePlugIn(:))<1e-9);
    assert(isequal(r1,ref.r1Hat) && isequal(k2,ref.kappa2Hat));
    fprintf('%s: matched reference, two batches, two contrasts\n',model{1});
end
fprintf('CICADA TASK GLM SELF-TEST PASSED (synthetic numerical test only)\n');
end
