function generate_QQ_plot(estimation_path, data_path, output_path, ...
                                   fct_form, basis, knots, force_n, sim_indices, plot_pairs)
    % Generate QQ plots for specific simulations and i-j pairs after the fact
    %
    % INPUTS:
    %   estimation_path: Path to estimation results
    %   data_path:       Path to simulated data
    %   output_path:     Path where plots will be saved
    %   fct_form:        Function form (e.g., 'LogNormal', 'Pareto', 'Estimates')
    %   basis:           Basis function (e.g., 'LogNormal', 'Spline')
    %   knots:           Number of knots (typically 1)
    %   force_n:         Force n option (e.g., 'n0', 'n1')
    %   sim_indices:     Vector of simulation indices to plot (e.g., [1, 5, 10])
    %   plot_pairs:      Cell array of [ii, jj] pairs (e.g., {[1,2], [1,3], [2,3]})
    
    % Load the simulated data
    simfile = fullfile(data_path, "Simulated_data_" + fct_form + "_alt.mat");
    load(simfile, 'M', 'Q_lnx_all');
    
    % Set up p_grid (should match what was used in simulation)
    p_grid = (1:99)./100;
    
    fprintf('\nGenerating QQ plots...\n');
    fprintf('Function form: %s | Basis: %s | Knots: %d | Force n: %s\n', ...
            fct_form, basis, knots, force_n);
    fprintf('Simulations to plot: %s\n', mat2str(sim_indices));
    fprintf('i-j pairs to plot: ');
    for p = 1:numel(plot_pairs)
        fprintf('[%d,%d] ', plot_pairs{p}(1), plot_pairs{p}(2));
    end
    fprintf('\n\n');
    
    % Loop over simulations
    for s = 1:numel(sim_indices)
        sim_id = sim_indices(s);
        fprintf('Processing simulation %d (%d of %d)...\n', sim_id, s, numel(sim_indices));
        
        % Create data structure for this simulation
        Mi = M(:,:,sim_id);
        d = make_data_simulation(Mi, knots, 'cubic');
        
        % Run regQQ to get the regression results and stored data
        l = regQQ(d, Q_lnx_all(:, sim_id), p_grid, basis, knots, force_n);
        
        % Loop over i-j pairs
        for p = 1:numel(plot_pairs)
            pair = plot_pairs{p};
            ii = pair(1);
            jj = pair(2);
            
            fprintf('  Generating plot for pair [%d, %d]...\n', ii, jj);
            
            % Generate the plot
            figure_plotter_QQ_minimal(output_path, l, ii, jj, sim_id, ...
                                     fct_form, basis, knots);
        end
    end
    
    fprintf('\nAll plots saved to: %s\n', output_path);
    fprintf('Plot naming convention: QQ_scatter_[fct_form]_[basis]_k[knots]_sim[id]_i[ii]_j[jj].pdf\n');
end


%% Example usage:
% ----------------------------------------------------------------------

% % After running simulations, generate plots for specific cases
% estimation_path = 'path/to/estimation';
% data_path = 'path/to/data';
% output_path = 'path/to/plots';
% 
% fct_form = 'LogNormal';
% basis = 'LogNormal';
% knots = 1;
% force_n = 'n0';
% 
% % Choose which simulations and pairs to plot
% sim_indices = [1, 5, 10];           % Plot simulations 1, 5, and 10
% plot_pairs = {[1,2], [1,3], [2,3]}; % Plot pairs (1,2), (1,3), (2,3)
% 
% % Generate all plots
% generate_QQ_plots_posthoc(estimation_path, data_path, output_path, ...
%                           fct_form, basis, knots, force_n, ...
%                           sim_indices, plot_pairs);
% 
% % This will create PDF files like:
% %   QQ_scatter_LogNormal_LogNormal_k1_sim1_i1_j2.pdf
% %   QQ_scatter_LogNormal_LogNormal_k1_sim1_i1_j3.pdf
% %   QQ_scatter_LogNormal_LogNormal_k1_sim1_i2_j3.pdf
% %   QQ_scatter_LogNormal_LogNormal_k1_sim5_i1_j2.pdf
% %   ... etc


%% Quick test with a single simulation and pair:
% ----------------------------------------------------------------------

% % Load one simulation
% load('path/to/Simulated_data_LogNormal_alt.mat', 'M', 'Q_lnx_all');
% 
% % Set up
% Mi = M(:,:,1);  % First simulation
% d = make_data_simulation(Mi, 1, 'cubic');
% p_grid = (1:99)./100;
% 
% % Run regression
% l = regQQ(d, Q_lnx_all(:,1), p_grid, 'LogNormal', 1, 'n0');
% 
% % Generate plot for pair (1,2)
% figure_plotter_QQ_minimal('output/', l, 1, 2, 1, 'LogNormal', 'LogNormal', 1);
% 
% % Check output: output/QQ_scatter_LogNormal_LogNormal_k1_sim1_i1_j2.pdf