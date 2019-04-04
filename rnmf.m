%% Settings RNMF %%



function ids = rnmf(name, nrVideos, varargin)

    
    % Global Parameter
    masksize =      'Adapt';
    method   = 'wholeVideo'; % movingWindow vs wholeVideo
    globalParams = readArgs({masksize,method}, varargin, {'masksize','method'});
    

    switch name
        
        case 'robustNMF'
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                                % Plot?
            defaultParams{2} = [2,3,5];                            % Rank
            defaultParams{3} = [40];                               % Iterations
            defaultParams{4} = [0, 0.05, 0.1];                     % Sparsity of S

            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
            
        case 'robustNMF_DS'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                    % Plot?
            defaultParams{2} = [2,5,10];                % Rank
            defaultParams{3} = [40];                   % Iterations
            defaultParams{4} = [0,0.1];                % Sparsity of S
            defaultParams{5} = [0, 0.5, 1];            % D_XY(S)
            defaultParams{6} = [0, 0.5];               % D_Z(S)


            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity','DxyS','DzS'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_DS2, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
        case 'robustNMF_NonlocalMeans'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                                 % Plot?
            defaultParams{2} = [2];                               % Rank
            defaultParams{3} = [40];                                % Iterations
            defaultParams{4} = [0];                            % Sparsity of S
            defaultParams{5} = [1];                               % Similarity
            defaultParams{6} = [2];                             % Similarity Type
            defaultParams{7} = [8];                                 % Splitsize

            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity','similarity','simtype','splitsize'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_NonlocalMeans, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            

            
        case 'robustNMF_L21'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                  % Plot?
            defaultParams{2} = [2,5,10];                    % Rank
            defaultParams{3} = [40];                               % Iterations
            defaultParams{4} = [0, 0.005, 0.01];% Sparsity of S



            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_L21, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
            
        case 'robustNMF_L1'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                  % Plot?
            defaultParams{2} = [2,3,5];              % Rank
            defaultParams{3} = [10];                 % Iterations
            defaultParams{4} = [0,1,2];              % Sparsity of S


            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_L1B, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
        case 'robustNMF_TVS'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                  % Plot?
            defaultParams{2} = [10];               % Rank
            defaultParams{3} = [40];                 % Iterations
            defaultParams{4} = [0.1];            % Sparsity of S
            defaultParams{5} = [0, 0.05, 0.1, 0.2];            % TV of S
            defaultParams{6} = [0];               % TV of S in time?


            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity','TvS','time'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_TVS, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
        case 'robustNMF_TVW'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                  % Plot?
            defaultParams{2} = [2];                  % Rank
            defaultParams{3} = [40];                 % Iterations
            defaultParams{4} = [0.1];           % Sparsity of S
            defaultParams{5} = [0, 10];            % TV of W

            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity','TvW'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_TVW2, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
        case 'robustNMF_excludeWHS'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                  % Plot?
            defaultParams{2} = [6];               % Rank
            defaultParams{3} = [40];                 % Iterations
            defaultParams{4} = [0];      % Sparsity of S
            defaultParams{5} = [0.5]; % Excluding

            % Create all Combinations of paramters
            paramsNames = {'plot','rank','iterations','sparsity','excluding'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            numberComb = evaluationRnmf(@robustNMF_excludeWHS, params, paramsNames, nrVideos, globalParams{2}, globalParams{1});
            ids = 1:numberComb;
            
        otherwise
            disp('Method not availiable.');
            ids = -1;
    end
end


