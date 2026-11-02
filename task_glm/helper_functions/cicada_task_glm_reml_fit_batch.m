function R=cicada_task_glm_reml_fit_batch(Y,cache)
% CICADA_TASK_GLM_REML_FIT_BATCH  Fit one original-BOLD voxel batch.
%
% cache = cicada_task_glm_reml_prepare(X0,C0,struct('model','AR2'));
% R = cicada_task_glm_reml_fit_batch(Y64,cache);
%
% Y is N x V ORIGINAL BOLD (not CICADA cleaned or Group output).
% Uses scan-level covariance and coefficient operators from PREPARE;
% computes voxelwise REML grid minima, GLS betas, COPE, and provisional
% plug-in VARCOPE. No I/O, outputs to disk, or automatic model selection.
%
% For a fixed original-time covariance A, GLS residual sum of squares is
% (U'Y)' (U'AU)^(-1) (U'Y). Thus the REML quadratic form cached at the
% selected parameter is also the GLS RSS and gives sigma2 = q/(N-P).
% This identity is why full per-voxel rewhitening is unnecessary here.

assert(isstruct(cache) && isscalar(cache) && ...
    isfield(cache,'status') && ...
    strcmp(cache.status,'PREPARED_FIXED_ICA_CONDITIONAL'), ...
    'Expected a cache from cicada_task_glm_reml_prepare');
assert(isnumeric(Y) && ismatrix(Y), ...
    'Y must be a numeric N x V matrix');
Y=double(Y);
[N,V]=size(Y);
assert(N==cache.N && V>=1 && V<=cache.maxBatchVoxels && ...
    all(isfinite(Y(:))), ...
    'Y time dimension, batch size or nonfinite values invalid');
D=cache.residualDimension;
Z=cache.U'*Y;
assert(all(sum(Z.^2,1)>0), ...
    'Zero restricted residual energy for at least one voxel');

bestScore=inf(1,V);
bestIndex=zeros(1,V);
bestQuadratic=nan(1,V);
for g=1:cache.nGridCandidates
    W=cache.cholRestricted(:,:,g)\Z;
    q=sum(W.^2,1);
    assert(all(isfinite(q) & q>0), ...
        'Invalid restricted quadratic at candidate %d',g);
    score=cache.logDetRestricted(g)+D*log(q/D);
    take=score<bestScore;
    bestScore(take)=score(take);
    bestIndex(take)=g;
    bestQuadratic(take)=q(take);
end
assert(all(bestIndex>0) && all(isfinite(bestScore)) && ...
    all(isfinite(bestQuadratic)),'REML profile fitting failed');

P=cache.P;
K=cache.nContrasts;
beta=zeros(P,V);
cope=zeros(K,V);
varcope=zeros(K,V);
sigma2=bestQuadratic/D;
for g=unique(bestIndex)
    cols=find(bestIndex==g);
    beta(:,cols)=cache.betaOperator(:,:,g)*Y(:,cols);
    cope(:,cols)=cache.C0*beta(:,cols);
    varcope(:,cols)=cache.contrastGeometry(:,g)*sigma2(cols);
end
assert(all(isfinite(beta(:))) && all(isfinite(cope(:))) && ...
    all(isfinite(varcope(:)) & varcope(:)>0), ...
    'Invalid GLS beta/COPE/VARCOPE');

rHat=cache.params(bestIndex,1)';
kHat=cache.params(bestIndex,2)';
R=struct();
R.status='PROVISIONAL_FIXED_ICA_CONDITIONAL';
R.model=cache.model;
R.N=cache.N;
R.nVoxels=V;
R.designRank=P;
R.residualDimension=D;
R.nContrasts=K;
R.nGridCandidates=cache.nGridCandidates;
R.r1Grid=cache.r1Grid;
R.k2Grid=cache.k2Grid;
R.r1Hat=rHat;
R.kappa2Hat=kHat;
R.remlProfileScore=bestScore;
R.bestGridIndex=bestIndex;
R.r1AtGridBoundary=(rHat==cache.r1Grid(1) | ...
    rHat==cache.r1Grid(end));
R.kappa2AtGridBoundary=false(1,V);
if strcmp(cache.model,'AR2')
    R.kappa2AtGridBoundary=(kHat==cache.k2Grid(1) | ...
        kHat==cache.k2Grid(end));
end
R.beta=beta;
R.cope=cope;
R.varcopePlugIn=varcope;
R.residualScale2=sigma2;
R.scaledWhitenedDesignCondition=cache.scaledCondition(bestIndex);
R.rcondWhitenedR=cache.rcondR(bestIndex);
R.notes=['Provisional fixed-ICA plug-in variances; no calibrated t, ' ...
    'p or z inference. No data written by this function.'];
end
