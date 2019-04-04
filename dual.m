%% Settings Dual %%



function dual(name, nrVideos, varargin)



    switch name
        
        case 'robustNMF_GL'
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                                % Plot?
            defaultParams{2} = [2];                                % Rank
            defaultParams{3} = [150];                               % Iterations
            defaultParams{4} = [0];                                % Sparsity of S
            defaultParams{5} = [0.05,1];                         % eps of S

            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity','epsilon'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            evaluationDual(@robustNMF_GL, params, paramsNames, nrVideos);
            
     
            
        otherwise
            disp('Method not availiable.');
    end
end


