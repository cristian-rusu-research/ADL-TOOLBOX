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

% setup of learning parameters
params.D = d;
params.Y = Y;
params.ITER_MAXIMUM = 30;
params.N = Nmax;
params.DICTIONARY_SIZE_MODE = 'ADAPTIVE';
params.K = 1;
params.K_MIN = 3;
params.K_MAX = 4*params.D;
params.SPARSITY_MODE = 'ADAPTIVE';
params.SPARSITY_NOISE_RISK = 1/2;
params.S = 1;

% parameters for pruning coherent atoms 
params.COHERENCE_MAX = 0.7;
params.func_prune_unused_atoms = @default_prune_max_coherent_atoms;

% parameters for pruning unused atoms 
M = d*log(d);
params.THRESHOLD_REMOVE_UNUSED_ATOMS = 1/M; %%% error if removed
params.TAU = 2*log(2*params.N/M);         %%% error if removed
params.func_atoms_scores = @atom_scores_count;
params.func_prune_unused_atoms = @default_prune_unused_atoms;
params.func_atoms_scores = @atom_scores_weighted_count;
params.func_sparse_approximation = @adaptive_pursuit;
params.func_dictionary_update = @itkrm_batch;
params.func_generate_atoms = @generate_repcand_itkrm_batch;

params.func_preprocess = @remove_mean;
params.func_initialization = @random_initialization;

params.func_statistics = {@dictionary_size @coherence};
params.func_stop_early = @alwaysno;

[out, in] = dla_toolbox(params);

disp(['Average time of the sparse approx. ' num2str(mean(out.times.sparse_approximation))]);
disp(['Average time of the dictionary update ' num2str(mean(out.times.dictionary_update))]);


figure; plot([out.statistics.dictionary_size]);
xlabel('Iteration');
ylabel('Size of dictionary');
figure; plot([out.statistics.coherence]);
xlabel('Iteration');
ylabel('Coherence');


[~,p] = sort(in.ATOM_SCORES(end,:),'descend');
dico = [ones(d,1)/sqrt(d), in.DICT(:,p)];

figure; imagesc(showdico(dico)); axis off
