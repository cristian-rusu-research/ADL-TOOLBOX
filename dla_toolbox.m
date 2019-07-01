function [out, in] = dla_toolbox(in)
%ADL_TOOLBOX
%   This algorithm ...

if(nargin ~= 1)
    error('syntax: [out, in] = dla_toolbox(in)');
end

%% make sure we show warnings
warning('on', 'all');

%% check that all directories exist
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%
%%%%%%%%%%

%% we need this, internal functions
addpath('./private_functions');

%% initialization functions in case dataset has to be synthetic
addpath('./initialization');


%% generate a synthetic dataset
if (isfield(in, 'Y'))
    if ischar(in.Y) && strcmp(clean_string(in.Y), 'SYNTHETIC')

        if (~isfield(in, 'SYNTHETIC_DICTIONARY_TYPE'))
            in.SYNTHETIC_DICTIONARY_TYPE = 'RANDOM';
            warning('SYNTHETIC_DICTIONARY_TYPE set to ''RANDOM');
        end

        if (isfield(in, 'D'))
            in.SYNTHETIC_DATA_DIMENSION = in.D;
            warning(['SYNTHETIC_DATA_DIMENSION set to ' num2str(in.D)]);
        end
        
        if (~isfield(in, 'SYNTHETIC_DATA_DIMENSION'))
            in.SYNTHETIC_DATA_DIMENSION = 128;
            warning(['SYNTHETIC_DATA_DIMENSION set to ' num2str(in.SYNTHETIC_DATA_DIMENSION)]);
        end
        in.D = in.SYNTHETIC_DATA_DIMENSION;

        if (isfield(in, 'D') && isfield(in, 'SYNTHETIC_DATA_DIMENSION'))
            if (in.D ~= in.SYNTHETIC_DATA_DIMENSION)
                error('Synthetic data dimension does not match data dimension.');
            end
        end

        if (~isfield(in, 'SYNTHETIC_DICTIONARY_SIZE'))
            if isfield(in, 'K') && isnumeric(in.K)
                in.SYNTHETIC_DICTIONARY_SIZE = in.K;
            else
                in.SYNTHETIC_DICTIONARY_SIZE = round(1.5*in.SYNTHETIC_DATA_DIMENSION);
            end
            warning(['SYNTHETIC_DICTIONARY_SIZE set to ' num2str(in.SYNTHETIC_DICTIONARY_SIZE)]);
        end

        if (~isfield(in, 'SYNTHETIC_DATA_SPARSITY'))
            if isfield(in, 'S') && isnumeric(in.S)
                in.SYNTHETIC_DATA_SPARSITY = in.S;
            else
                in.SYNTHETIC_DATA_SPARSITY = 6;
            end
            warning(['SYNTHETIC_DATA_SPARSITY set to ' num2str(in.SYNTHETIC_DATA_SPARSITY)]);
        end

        if (~isfield(in, 'SYNTHETIC_DATA_SPARSITY_RANGE'))
            in.SYNTHETIC_DATA_SPARSITY_RANGE = 0;
            warning(['SYNTHETIC_DATA_SPARSITY_RANGE set to ' num2str(in.SYNTHETIC_DATA_SPARSITY_RANGE)]);
        end

        if (~isfield(in, 'SYNTHETIC_DATA_DECAY'))
            in.SYNTHETIC_DATA_DECAY = 0.1;
            warning(['SYNTHETIC_DATA_DECAY set to ' num2str(in.SYNTHETIC_DATA_DECAY)]);
        end

        if (~isfield(in, 'SYNTHETIC_DATA_NOISE_LEVEL'))
            in.SYNTHETIC_DATA_NOISE_LEVEL = 0.25/sqrt(in.D);
            warning(['SYNTHETIC_DATA_NOISE_LEVEL set to ' num2str(in.SYNTHETIC_DATA_NOISE_LEVEL)]);
        end

        if (~isfield(in, 'SYNTHETIC_DATA_OUTLIERS'))
            in.SYNTHETIC_DATA_OUTLIERS = 0.05;
            warning(['SYNTHETIC_DATA_OUTLIERS set to ' num2str(in.SYNTHETIC_DATA_OUTLIERS)]);
        end

        if (~isfield(in, 'SYNTHETIC_DATA_NORMALIZE'))
            in.SYNTHETIC_DATA_NORMALIZE = 0;
            warning(['SYNTHETIC_DATA_NORMALIZE set to ' num2str(in.SYNTHETIC_DATA_NORMALIZE)]);
        end

        if strcmp(in.SYNTHETIC_DICTIONARY_TYPE, 'RANDOM')
            in.SYNTHETIC_DICT = randn(in.SYNTHETIC_DATA_DIMENSION, in.SYNTHETIC_DICTIONARY_SIZE);
            in.SYNTHETIC_DICT = bsxfun(@rdivide, in.SYNTHETIC_DICT, sqrt(sum(in.SYNTHETIC_DICT.^2)));
        else
            if strcmp(in.SYNTHETIC_DICTIONARY_TYPE, 'DCT_DIRAC')
                in.SYNTHETIC_DICTIONARY_SIZE = min(in.SYNTHETIC_DICTIONARY_SIZE, 2*in.SYNTHETIC_DATA_DIMENSION);
                in.SYNTHETIC_DICT = [eye(in.SYNTHETIC_DATA_DIMENSION), dctmtx(in.SYNTHETIC_DATA_DIMENSION)'];
                p = randperm(2*in.SYNTHETIC_DATA_DIMENSION);
                in.SYNTHETIC_DICT = in.SYNTHETIC_DICT(:, p(1:in.SYNTHETIC_DICTIONARY_SIZE));
            end
        end

        if ~isfield(in, 'func_generate_synthetic_dataset')
            in.func_generate_synthetic_dataset = @default_generate_synthetic_dataset;
            warning('func_generate_synthetic_dataset set to default_generate_synthetic_dataset');
        end

        in = in.func_generate_synthetic_dataset(in);
    end
end


%% check if dataset is OK
if (~isfield(in, 'Y'))
    error('Dataset in.Y not there.');
end
if (any(any(isnan(in.Y))))
    error('Dataset in.Y contains NaN.');
end
if (any(any(isinf(in.Y))))
    error('Dataset in.Y contains Inf.');
end


%% remove zero columns
energy = sum(in.Y.^2); emptydatapoints = find (energy == 0);
in.Y(:, emptydatapoints) = []; energy(emptydatapoints) = [];
clear emptydatapoints; clear energy;


%% get size of dataset
[in.D, in.N] = size(in.Y);

if (in.N < 1)
    error('Dataset in.Y is empty.');
end


%% default parameters
% is the number of atoms set?
if (~isfield(in, 'K'))
    in.K = 2*in.D;
    warning(['Dictionary size set to ' num2str(in.K)]);
end

% check the size of the dataset vs the size of the dictionary
if (isnumeric(in.K))
    in.K = round(in.K);
    % make sure 1 <= K <= N
    in.K = max(1, in.K); in.K = min(in.K, in.N);
    
    if (in.N < in.K + 1)
        disp('Less training signals than atoms => trivial solution');
        in.DICT = bsxfun(@rdivide, in.Y, sqrt(energy));
        in.X = sparse(1:in.N, 1:in.N, sqrt(energy));
        return;
    end
end

% check the dictionary size mode
if (~isfield(in, 'DICTIONARY_SIZE_MODE'))
    in.DICTIONARY_SIZE_MODE = 'FIXED';
    if (ischar(in.K))
        if (strcmp(clean_string(in.K), 'ADAPTIVE')) || (strcmp(clean_string(in.K), 'ADAPT')) || (strcmp(clean_string(in.K), 'A'))
            in.DICTIONARY_SIZE_MODE = 'ADAPTIVE';
        end

        if (strcmp(clean_string(in.K), 'REPLACEMENT')) || (strcmp(clean_string(in.K), 'REPLACE')) || (strcmp(clean_string(in.K), 'R'))
            in.DICTIONARY_SIZE_MODE = 'REPLACEMENT';
        end
    end
else
    %% make sure DICTIONARY_SIZE_MODE is either A, R, or F
end

% is the sparsity set?
if (~isfield(in, 'S'))
    in.S = min(round(log2(in.D)), in.K);
    warning(['Sparsity set to ' num2str(in.S)]);
end

% make sure 1 <= S <= min(D, K)
if (isnumeric(in.S)) 
    in.S = round(in.S);
    in.S = max(1, in.S); 
    in.S = min(in.S, min(in.D, in.K));
    %in.SPARSITY_MODE = 'FIXED';
end

% check the sparsity mode
if (~isfield(in, 'SPARSITY_MODE'))
    in.SPARSITY_MODE = 'FIXED';
    if (ischar(in.S))
        if (strcmp(clean_string(in.S), 'ADAPTIVE') || strcmp(clean_string(in.S), 'A'))
            in.SPARSITY_MODE = 'ADAPTIVE';
        end        
    end
else
    if (strcmp(clean_string(in.SPARSITY_MODE), 'ADAPTIVE') || strcmp(clean_string(in.SPARSITY_MODE), 'A'))
        in.SPARSITY_MODE = 'ADAPTIVE';
    else
        in.SPARSITY_MODE = 'FIXED';
    end
end

% if the dictionary size mode is adaptive and you do not supply a maximum
% number of iterations, we consider that the algorithm continues until the
% size of the dictionary reaches (or passes) K_MIN
if strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') && (isfield(in, 'K_MIN'))
    if (isfield(in, 'ITER_MAXIMUM'))
        if (~isinf(in.ITER_MAXIMUM))
%             warning('Do not set both K_MIN and ITER_MAXIMUM! Running for as many iterations as necessary to reach K_MIN.');
%             in.ITER_MAXIMUM = 50;
            
            if (~isfield(in, 'ITER_POLISHING_STEPS'))
                in.ITER_POLISHING_STEPS = 200;
            end
        end
    end
    
    if (in.K_MIN == in.K_MAX)
%         in.ITER_MAXIMUM = 50;
        
        if (~isfield(in, 'ITER_POLISHING_STEPS'))
            in.ITER_POLISHING_STEPS = 200;
        end
    end
end

% is the maximum number of iterations set?
if (~isfield(in, 'ITER_MAXIMUM'))
    in.ITER_MAXIMUM = 100;
    warning(['Maximum number of iterations set to ' num2str(in.ITER_MAXIMUM)]);
else
    if (isnumeric(in.ITER_MAXIMUM))
        in.ITER_MAXIMUM = max(in.ITER_MAXIMUM, 10);
    else
        in.ITER_MAXIMUM = 100;
    end
end

%% default algoritm
addpath('./preprocess');
if (~isfield(in, 'func_preprocess'))
    in.func_preprocess = @remove_mean;
    warning('Pre-process function set to @remove_mean');
end

if (~isfield(in, 'func_initialization'))
    in.func_initialization = @random_initialization;
    warning('Initization function set to @random_initialization');
end

addpath('./dictionary_update');
if (~isfield(in, 'func_dictionary_update'))
    in.func_dictionary_update = @itkrm_batch;
    warning('Dictionary update function set to @itkrm_batch');
end

addpath('./sparse_approximation');
if (~isfield(in, 'func_sparse_approximation'))
    in.func_sparse_approximation = @thresholding;
    warning('Sparse approximation function set to @thresholding');
end

addpath('./objective_function');
if (~isfield(in, 'func_statistics'))
    in.func_statistics = @frobenius_norm_squared;
    warning('Statistics function set to @frobenius_norm_squared');
end

if (~isfield(in, 'func_stop_early'))
    in.func_stop_early = @alwaysno;
end

if (~isfield(in, 'func_do_every_iteration'))
    in.func_do_every_iteration = @donothing;
end

% initialization for adaptive dictionary mode
if (strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') || strcmp(in.DICTIONARY_SIZE_MODE, 'REPLACEMENT'))
    addpath('./atoms_prune');
    addpath('./train_replacement_atoms');
    
    if (~isfield(in, 'THRESHOLD_REMOVE_UNUSED_ATOMS'))
        in.THRESHOLD_REMOVE_UNUSED_ATOMS = 50;
        warning(['THRESHOLD_REMOVE_UNUSED_ATOMS set to ' num2str(in.THRESHOLD_REMOVE_UNUSED_ATOMS)]);
    end
    
    if (~isfield(in, 'TAU'))
        in.TAU = 2*log(2*in.N/50);
        warning(['TAU set to ' num2str(in.TAU)]);
    end

    if (strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE'))
        % start from a small dictionary
        if (~isfield(in, 'K')) || ((isfield(in, 'K') && ischar(in.K)))
            in.K = round(in.D);
        end

        if (~isfield(in, 'L'))
            in.L = round(log2(in.D));
        end

        if (~isfield(in, 'K_MIN'))
            in.K_MIN = in.K;
        end
        
        if (~isfield(in, 'K_MAX'))
            in.K_MAX = 4*in.D;
        end
    end
    
    if (strcmp(in.DICTIONARY_SIZE_MODE, 'REPLACEMENT'))
        
        if (~isfield(in, 'L'))
            in.L = round(log2(in.D));
        end

        if (~isfield(in, 'K_MIN'))
            in.K_MIN = 4*in.D;
        end
        
        if (~isfield(in, 'K_MAX'))
            in.K_MAX = 4*in.D;
        end
        
        if (in.K_MIN ~= in.K_MAX)
            error('In REPLACEMENT mode K_MIN = K_MAX');
        end
        
        in.K = in.K_MIN;
        
        if (~isfield(in, 'L'))
            warning('In REPLACEMENT mode L is set adaptively with each iteration! Ignoring current value.');
            in.L = 0;
        end
    end
    
    if (~isfield(in, 'func_check_if_remove_unused_atoms'))
        in.func_check_if_remove_unused_atoms = @rule_remove_unused_atoms;
        
        if (~isfield(in, 'func_atoms_scores'))
            in.func_atoms_scores = @atom_scores_count;
        end

        if (~isfield(in, 'func_prune_unused_atoms'))
            in.func_prune_unused_atoms = @default_prune_unused_atoms;
        end
    end
    
    if (isfield(in, 'COHERENCE_MAX'))
        in.COHERENCE_MAX = min(in.COHERENCE_MAX, 1);
        in.COHERENCE_MAX = max(in.COHERENCE_MAX, 0);
        
        if (~isfield(in, 'func_check_if_remove_coherent_atoms'))
            in.func_check_if_remove_coherent_atoms = @rule_remove_coherent_atoms;

            if (~isfield(in, 'func_prune_max_coherent_atoms'))
                in.func_prune_max_coherent_atoms = @default_prune_max_coherent_atoms;
            end
        end
    end
    
    if (~isfield(in, 'func_check_if_add_atoms'))
        in.func_check_if_add_atoms = @rule_add_atoms;
        
        if (~isfield(in, 'func_add_atoms'))
            in.func_add_atoms = @add_atoms;
        end
    end
    
end

% initialization for adaptive sparsity mode
if (strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))% && ~strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE'))
    in.S = 'ADAPTIVE';
    
    if (~isfield(in, 'S_MIN'))
        in.S_MIN = 1;%floor(log2(min(in.K, in.D))/2);
        warning(['S_MIN is set to ' num2str(in.S_MIN)]);
    end

    if (~isfield(in, 'S_MAX'))
        in.S_MAX = floor(in.D/log(in.D));
        warning(['S_MAX is set to ' num2str(in.S_MAX)]);
    end
end

%% get out of the way, toolbox starts
% preprocess the data
tic; in = in.func_preprocess(in); out.times.preprocess = toc;

% initialize the dictionary and sparse representations
tic; in = in.func_initialization(in); out.times.initialization = toc;

out.times.sparse_approximation = zeros(in.ITER_MAXIMUM, 1);
out.times.dictionary_update = zeros(in.ITER_MAXIMUM, 1);
out.times.adaptive_dictionary_size = zeros(in.ITER_MAXIMUM, 1);

if strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE')
    in.K_PER_ITERATION = zeros(in.ITER_MAXIMUM, 1);
end

if (strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') || strcmp(in.DICTIONARY_SIZE_MODE, 'REPLACEMENT'))
    in.ATOM_SCORES = [];
end

if (strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') && strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
%     error('Not done yet!');
end

%%%%%%% TO DO PUT THIS SOMEWHERE
if (~isfield(in, 'ATOM_SCORES_MEMORY'))
    in.ATOM_SCORES_MEMORY = round(log2(in.D));
end

% the iterative process
in.ITER_CURRENT = 0;
while (in.ITER_CURRENT < in.ITER_MAXIMUM)
    in.ITER_CURRENT = in.ITER_CURRENT + 1;
    
    % update the sparse representations
    tic; in = in.func_sparse_approximation(in); out.times.sparse_approximation(in.ITER_CURRENT) = toc;
    
    %%% suggestion: for minimal invasiveness
    %%% if mode = adaptive or replacement call generate_repcand here
    %%% where both the old dictionary with true coefficients 
    %%% and residuals are still available to generate L candidates
    if strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') || strcmp(in.DICTIONARY_SIZE_MODE, 'REPLACEMENT')
%         in = generate_repcand_itkrm_batch(in);
        in = in.func_generate_atoms(in);
    end
    
    % update the dictionary
    tic; in = in.func_dictionary_update(in); out.times.dictionary_update(in.ITER_CURRENT) = toc;
    
    % if dictionary size is adaptive
    tic;
    if (strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') || strcmp(in.DICTIONARY_SIZE_MODE, 'REPLACEMENT'))
        % update atom scores
        in = in.func_atoms_scores(in);
        
%         old_K = in.K;
        if (isfield(in, 'COHERENCE_MAX'))
            if (in.func_check_if_remove_coherent_atoms(in))
                % if maximum coherence is set, prune highly coherent atoms
                in = in.func_prune_max_coherent_atoms(in);
            end
        end
        
        if (in.func_check_if_remove_unused_atoms(in))
            % remove unused atoms 
            in = in.func_prune_unused_atoms(in);
        end
        
%         if (strcmp(in.DICTIONARY_SIZE_MODE, 'REPLACEMENT'))
%             in.L = old_K - in.K;
%         end
%         clear old_K;
        
        if (strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') && in.func_check_if_add_atoms(in)) ...
                    || (in.K < in.K_MIN)
            % update sparse representations before adding new atoms,
            % func_add_atoms needs to have update information
%             in = in.func_sparse_approximation(in);
            
            % add new, good, atoms          
            %%%% with the new generic function that adds atoms
            in = in.func_add_atoms(in);

            % update sparse representations after adding new atoms
%             in = in.func_sparse_approximation(in);
        end
        
        % if sparsity mode is adaptive, update S_MAX is case in.K has
        % changed
        if (strcmp(in.SPARSITY_MODE, 'ADAPTIVE'))
            in.S_MAX = floor(min(in.K, in.D)/2);
        end
    end
    out.times.adaptive_dictionary_size(in.ITER_CURRENT) = toc;
    
    % keep track of statistics
    aux = struct;
    if (length(in.func_statistics) > 1)
        for index = 1:length(in.func_statistics)
            aux = in.func_statistics{index}(in, aux);
        end
    else
        aux = in.func_statistics(in, aux);
    end
    
    if (~isempty(aux))
        out.statistics(in.ITER_CURRENT) = aux;
    end
    clear aux;
    
    in.statistics = out.statistics;
    
    % stop early because we converged or some other condition is true
    if (in.func_stop_early(in))
        out.statistics = out.statistics(1:in.ITER_CURRENT);
        break;
    end
    
    % maximum dictionary size reached, stop and polish the dictionary we
    % have for a fixed number of extra iterations no longer growing the
    % dictionary
    if strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE') && (in.K >= in.K_MAX)
        % resize dictionary and sparse representations
        if (in.K > in.K_MAX)
            in.K = in.K_MAX;
            in.DICT = in.DICT(:, 1:K_MAX);
            in.X = in.X(1:K_MAX, :);
        end
        
        % prescribed dictionary size reached, enter replacement mode
        if isinf(in.ITER_MAXIMUM)
            in.DICTIONARY_SIZE_MODE = 'REPLACEMENT';
%             in.ITER_MAXIMUM = in.ITER_CURRENT + in.ITER_POLISHING_STEPS;
        end
    end

    % after the first 1/2 iterations, also allow sparsity 0 (i.e., check for
    % pure noise data points)
    if (strcmp(in.SPARSITY_MODE, 'ADAPTIVE')) && (in.ITER_CURRENT >= round(1/2*in.ITER_MAXIMUM))
        in.S_MIN = 0;
    end
    
    % after the first 3/4 iterations, when sparsity mode is adaptive check
    % for the average sparisty. If that is too low, decide that the data
    % does not have a sparse representation in the dictionary
    if (strcmp(in.SPARSITY_MODE, 'ADAPTIVE')) && (in.ITER_CURRENT >= round(3/4*in.ITER_MAXIMUM))
        if (nnz(in.X)/in.N < 1)
            in.DICT = [];
            in.X = [];
            warning('Your data does not seem to be sparsely representable in a dictionary! Returning an empty solution.');
            break;
        end
    end
    
    % in case of adaptive dictionary size, check the evolution of the size
    if strcmp(in.DICTIONARY_SIZE_MODE, 'ADAPTIVE')
        in.K_PER_ITERATION(in.ITER_CURRENT) = in.K;
        
        % if in the past 100 iterations the size of the dictionary mostly
        % decreased (rather than increasing) then put a warning and
        % increase the number of candidates
        if (in.ITER_CURRENT >= 100)
            if (sum(diff(in.K_PER_ITERATION(end-99:end))) < 1)
                %%% TO RE-THINK THIS IN THE FUTURE
                % WHAT DO WE DO IF THE ITERATIONS ARE GOING ON BUT THE SIZE
                % OF THE DICTIONARY DOES NOT INCREASE
%                 in.L = max(round(1.1*in.L), in.L+1);
%                 warning(['The size K of your dictionary does not increase! Increasing L to ' num2str(in.K) '.']);
            end
        end
    end
    
    if (length(in.func_do_every_iteration) > 1)
        for index = 1:length(in.func_do_every_iteration)
            in = in.func_do_every_iteration{index}(in);
        end
    else
        in = in.func_do_every_iteration(in);
    end
end

%% for DEBUG purposes
STOP = 1;
end
