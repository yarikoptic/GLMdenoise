% Author: Kendrick Kay
% Source: private correspondence

% run GLMdenoise (initial call)
results = GLMdenoisedata(model,data,stimdur,tr,[],[],[],'GLMdenoisefigures');

% get the noise regressors
noisereg = cellfun(@(x) x(:,1:results.pcnum),results.pcregressors,'UniformOutput',0);

% get the initial HRF seed
hrfknobs = normalizemax(getcanonicalhrf(stimdur,tr)');

% initialize
hrfs = zeros(18,prod(xyzsize));
betas = zeros(prod(xyzsize),56);
R2 = zeros(xyzsize);

% do it
cache = [];
tic;
for zz=1:prod(xyzsize)

  % extract voxel time-series
  data0 = cellfun(@(x) subscript(squish(x,3),{zz ':'}),data,'UniformOutput',0);

  % call GLMdenoise
  [results0,cache] = GLMestimatemodel(model,data0,stimdur,tr,'optimize',hrfknobs,0, ...
                                      struct('extraregressors',{noisereg},'hrfthresh',-Inf, ...
                                             'suppressoutput',1),cache);

  % record results
  hrfs(:,zz) = results0.modelmd{1};
  betas(zz,:) = results0.modelmd{2};
  R2(zz) = results0.R2;

  % report time
  if mod(zz,100)==0
    etime = toc;
    ctime = (prod(xyzsize)-zz)/100*etime/60;
    fprintf('analysis for 100 voxels took %.2f seconds (estimated time until completion: %.1f minutes).\n',etime,ctime);
    tic;
  end

end
