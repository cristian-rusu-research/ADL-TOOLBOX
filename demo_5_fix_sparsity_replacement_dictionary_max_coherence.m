%%% start fresh
% if there is a variable called CLEAN_UP, then check if we start fresh
if (~exist('CLEAN_UP', 'var')) || (CLEAN_UP ~= 0)
    close all
    clear
    clc
end

rng(1337);                  %%% random seed

picname = 'mandrill.png';   %%% picture to extract patches from
s1 = 8;                     %%% patchheight
s2 = s1;                    %%% patchwidth
d = s1*s2;   
N = 1000;                   %%% signals per iteration

%%%% load picture 
pic = imread(strcat('images/',picname));
pic = im2double(pic);
[d1,~,d3]=size(pic);
if d3 > 1
    pic = 0.2989*pic(:,:,1)+0.5870*pic(:,:,2)+0.1140*pic(:,:,3);
end
pic = imresize(pic,256/d1);   %%% get 256x256 or so 
[d1,d2] = size(pic);
Nmax = (d1-s1+1)*(d2-s2+1);
N = min(Nmax,N);
%%%% end of loading picture

%%% input data for dictionary learning
[Y, ~] = pic2patches(pic,s1,s2);

%%% set parameters for toolbox
params = [];
out =  []; in = [];
Y = Y(:, 1:5000);             %%% restrict dataset to 5000 signals
params.Y = Y;

% number of atoms
params.K = 'R';     %%% replacement mode
params.K_MIN = 80;
params.K_MAX = 80;

params.COHERENCE_MAX = 0.5; %%% set the maximum allowed coherence

params.S = 4;   % fix the sparsity level
params.L = 2;   % number of candidate atoms

params.func_preprocess = @remove_mean;
params.func_initialization = @random_initialization;
params.func_dictionary_update = @ksvd_batch;
params.func_generate_atoms = @generate_repcand_itkrm_batch;
params.func_sparse_approximation = @thresholding;
params.func_atoms_scores = @atom_scores_count;
params.func_prune_unused_atoms = @prune_least_used_atom;
params.func_statistics = {@frobenius_norm_squared_dictionary_size @coherence};
params.func_do_every_iteration = @plotdico;
params.func_stop_early = @alwaysno;

%%% call the toolbox
[out, in] = dla_toolbox(params);

disp(['Average time of the sparse approx. ' num2str(mean(out.times.sparse_approximation))]);
disp(['Average time of the dictionary update ' num2str(mean(out.times.dictionary_update))]);

figure; plot([out.statistics.frobenius_norm_squared]);
xlabel('Iteration');
ylabel('Objective function value');
figure; plot([out.statistics.dictionary_size]);
xlabel('Iteration');
ylabel('Size of dictionary');
figure; plot([out.statistics.coherence]);
xlabel('Iteration');
ylabel('Coherence');
