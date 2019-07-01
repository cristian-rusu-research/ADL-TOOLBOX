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

params.FIG_ERROR_TEXT = 'default alg';
params.func_do_every_iteration = @ploterror;

%%% call the toolbox
[out, in] = dla_toolbox(params);

disp(['Average time of the sparse approx. ' num2str(mean(out.times.sparse_approximation))]);
disp(['Average time of the dictionary update ' num2str(mean(out.times.dictionary_update))]);

figure; imagesc(showdico(in.DICT)); axis off
title('final dictionary atoms');
