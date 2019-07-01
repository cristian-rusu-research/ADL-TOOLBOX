close all
clear
clc

demos = {
            %'demo_0_blind_run.m', ...
            %'demo_1_ksvd_algorithm.m', ...
            %'demo_2_itkrm_algorithm.m', ...
            %'demo_3_hybrid.m', ...
            'demo_4_fix_sparsity_replacement_dictionary.m', ...
            'demo_5_fix_sparsity_replacement_dictionary_max_coherence.m', ...
            'demo_6_fix_sparsity_adapt_dictionary.m', ...
            'demo_7_adaptive_sparsity_and_dictionary.m', ...
            };

while(1)
%     order = randsample(1:length(demos), length(demos));
    order = 1:length(demos);
    
    CLEAN_UP = 0;
    
    for demo_i = order
        decodedValue = ['Running ' demos{order(demo_i)}];

        fig = figure;
        hPan = uipanel(fig,'Units','normalized');

        uicontrol(hPan, 'Style','text', 'HorizontalAlignment','center', ...
        'FontSize', 16, 'Units','normalized', 'Position',[0.2 0.4 0.6 0.2], ...
        'String',decodedValue);
    
        pause(3);
        
        close all;
        
        run(demos{order(demo_i)});
        
        fprintf('\nPress any key to continue with the next demo ...\n\n');
        pause();
        
        close all
    end
    
    STOP = 1;
end
