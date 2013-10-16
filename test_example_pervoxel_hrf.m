% create some fake but feasible data 
%addpath('/home/data/famface/try-customhrf/GLMdenoise')
%addpath('/home/data/famface/try-customhrf/GLMdenoise/utilities')

nnuisances = 3;                         % number of random nuisance variables to be mixed in through the voxels
anuisances = 0.6;
anoise = 0.2;

srandom_variance = 0.2;

baseline = 100;
stimdur = 1;
tr = 2.0;
bdur = round(8/tr); %  block duration in volumes -- not quite a gamma, right? ;)
opt = struct();
% opt.modelpv = 1;
opt.maxpolydeg = 0;                     %  for some reason even if virtually noiseless -- looks like hrf was accounting for some filtering of the model on the edges... thought that it might have been due to polynomials... wrong
opt.seed = 0;

design = {};
% very cruel simple case
% design = {{[2, 19, 30]; [8,12,27]}; {[0, 12, 29]; [6, 32, 37]}};
nruns = 3;
tsize = 100;
xyzsize = [4, 6, 6];

% somewhat more realistic although still stupid case
%nruns = 6;
%tsize = 100;
%xyzsize = [15, 16, 10];

% very silly design, I know
for r=1:nruns
    % 2 EVs
    design{r} = {1:12:tsize*tr - bdur; 5:14:tsize*tr - bdur};
end

% seed RNG so we get the same (let's hope) result across runs. rand('seed', 0) isn't effective btw
RandStream.setDefaultStream(RandStream('mt19937ar','seed', 0));
data = {};
for r=1:nruns
    data{r} = random('norm', 0, 1, [xyzsize, tsize]) * anoise;
end

data_signal = {};
data_nuisances = {};

% add some very basic block signal at those onsets
for chunk=1:length(design)
    data_signal{chunk} = zeros(size(data{chunk}));
    data_nuisances{chunk} = zeros(size(data{chunk}));
    for ev=1:size(design{chunk}, 1)
        for i=1:length(design{chunk}{ev})
            vol = round(design{chunk}{ev}(i)/tr)+1;
            % for both conditions
            uvol = min(tsize, vol+bdur-1);
            nvols = length(vol:uvol);
            % let's have 3 informative voxels
            data_signal{chunk}(ev,2,3,vol:uvol) = data_signal{chunk}(ev,2,3,vol:uvol) + 1 + randn(1, 1, 1, nvols)*srandom_variance;
            % first two with separate EV in each, while 3rd one would have all mixed equally
            data_signal{chunk}(3,2,3,vol:uvol) = data_signal{chunk}(3,2,3,vol:uvol) + 1 + randn(1, 1, 1, nvols)*srandom_variance;
        end
    end
    % add some common noise (nuisances) components
    % to be mixed in with all the voxels with arbitrary weights
    nuisances = random('norm', 0, 1, nnuisances, tsize) * anuisances;
    nvoxels = prod(xyzsize);
    mixin = random('norm', 0, 1, nvoxels, nnuisances)*0.3;
    data_nuisances{chunk} = reshape(mixin*nuisances, [xyzsize, tsize]);
    data{chunk} = data{chunk} + data_nuisances{chunk} + data_signal{chunk} + baseline;
end
% just plotting first run as a sample
if 1
x = 1:tsize;
for i=1:3
    figure; plot(x, squeeze(data{1}(i,2,3,:))-baseline, 'k', ...
                 x, squeeze(data_nuisances{1}(i,2,3,:)), 'b-', ...
                 x, squeeze(data_signal{1}(i,2,3,:)), 'r');
    title(i)
end
end

%figure; plot(squeeze(data{1}(1,2,3,:)))
%figure; plot(squeeze(data{1}(2,2,3,:)))

% use design represented in type 2 -- arrays with onset volumes marked
% so we could check other than 'assume' types
% INTERESTING/TODO:  that results in quite different estimates, so we
% need to check if conversion is done correctly since here shouldn't matter...
% Actually estimated model (according to generated figures) looks bad!
% # of nuisance PCs though estimated nicely/correctly
design = convert_design3to2(design, data, tr);

% run the beast
hrfknobs = [];                           % no assumption
hrfmodel = 'fir'; hrfknobs = 10;
[results_denoise, denoiseddata] = GLMdenoisedata(design, data, stimdur, tr, hrfmodel, hrfknobs, opt, 'GLMdenoisefigures')

if hrfmodel == 'fir'
 % plot for our two mighty EVs
 for ev=1:2
  ev_hrf = squeeze(results_denoise.modelmd(ev,2,3,ev,:));
  figure; stem(ev_hrf);
  legend('HRF... should be block starting at 0')
  ev_rec = conv(squeeze(design{1}(:, ev)), ev_hrf);
  figure; plot(x, squeeze(design{1}(:, ev)), 'g', ...
               x, ev_rec(1:tsize), 'b', ...
               x, squeeze(data{1}(ev,2,3,:)) - baseline, 'r', ...
               x, squeeze(denoiseddata{1}(ev, 2, 3, :)) - baseline, 'm');
  legend('onsets', 'reconstructed (up to scale)', 'original', 'denoised')
  title(ev)
 end
end
