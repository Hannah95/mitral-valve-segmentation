%% Settings RNMF %%



function segment(name, onMethod, ids, varargin)

    % Global Parameter
    postProcessing = 'both'; % both, true, false
    cropped = 'false'; % true or false
    globalParams = readArgs({postProcessing,cropped}, varargin, {'postProcessing','cropped'});
    


    switch name
        
        case 'segmentCVPlus'
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                                % Plot?
            defaultParams{2} = [20];                               % Iterations
            defaultParams{3} = [0.5];                              % Thresholding S
            defaultParams{4} = [0.04];                             % TVxy(M)
            defaultParams{5} = [0];                                % TVz(M)
            defaultParams{6} = [1];                                % Thresholding W*H
            defaultParams{7} = [0.1];                                % Area
            defaultParams{8} = [0.01];                             % t




            % Create all Combinations of paramters
            paramsNames = {'plot','iterations','threshS','TVxyM','TVzM','threshWH','area','thresh'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            evaluationSegment(@segmentCVPlus, params, paramsNames,onMethod, ids, globalParams{1},globalParams{2});
            
        case 'segmentMBO'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = [0];                                % Plot?
            defaultParams{2} = [50];                               % Iterations
            defaultParams{3} = [1];                       % alpha -> Sim Matrix
            defaultParams{4} = [5];                                % patchsizeXY
            defaultParams{5} = [0];                                % patchsizeZ
            defaultParams{6} = [100];                              % nrSimilarities
            defaultParams{10} = [50];                              % nrEV
            defaultParams{7} = [0.1];                                % dataterm
            defaultParams{8} = [1,2,3];                                % simType
            defaultParams{9} = [8];                                % splitSize




            % Create all Combinations of paramters
            paramsNames = {'plot','iterations','alpha','patchsizeXY','patchsizeZ',...
                           'nrSimilarities','dataterm','simType','splitSize','nrEV'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            evaluationSegment(@segmentMBO, params, paramsNames,onMethod, ids, globalParams{1},globalParams{2});
            
        case 'thresholding'
            
            % Default Parameters
            defaultParams = {};
            defaultParams{1} = linspace(0,1,40);                                % threshold

            % Create all Combinations of paramters
            paramsNames = {'threshold'};
            
            % Read Args
            params = readArgs(defaultParams, varargin, paramsNames);

            % Run Evaluation for all Parameters
            evaluationSegment(@thresholding, params, paramsNames,onMethod, ids, globalParams{1},globalParams{2});
            

            
        otherwise
            disp('Method not availiable.');
    end
end


