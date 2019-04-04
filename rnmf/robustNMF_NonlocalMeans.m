function ret = robustNMF_NonlocalMeans(X3D,vars)

% Energy = 0.5 * norm(X-W*H-S,'fro').^2+ lambda*sum(sum(abs(S))) - my/2 * S(:)' * L * S(:);

%% Get Parameters

var_plot = vars(1); % Plot?
r = vars(2);        % Rank
iter = vars(3);     % No iteration
lambda = vars(4);   % Sparsity
my = vars(5);       % Smoothness
simType = vars(6);  % Type of similarity measurement
splitSize = vars(7);% Splitsize

%% Precalculate W , H and S


sz = size(X3D,3);
numberSplit = ceil(sz/splitSize);
SFinal = [];
E = [];

for split = 1:numberSplit
    fBeg = (split - 1) * splitSize + 1;
    fEnd = min(fBeg+splitSize-1,sz);
    
    X3DTemp = X3D(:,:,fBeg:fEnd);
    X = reshape(X3DTemp, size(X3DTemp,1) * size(X3DTemp,2), size(X3DTemp,3));
    
    r = min(vars(2),size(X,2));
    [W,H] = nnmf(X,r);
    S = max(X - W * H, 0);
    
    %% Calculate Graph Laplace
    if my > 0
        opt={};
        opt.alpha = 1;
        opt.patchsizeXY = 7;
        opt.patchsizeZ = 0;
        opt.nrSim = 90;
        opt.method = 'dense';
        opt.simType = simType;
        opt.ensureRank = 0;
        [sim,opt] = getSimilarityMatrixNonlocalMeans(X3DTemp, opt);
        
        if strcmp(opt.method,'dense')
            XX = sim.XX;
            XY = sim.XY;
            [eigVal,eigVec] = nystroem(XX,XY,opt.nrSim);
                eigVec(opt.permutation,:) = eigVec;
            
%             eigVal = eigVal(1:nrEV,1:nrEV);
%             eigVec = eigVec(:,1:nrEV);
            
            L = @(u) eigVec * (eigVal * (eigVec' * u));
        else
            L = getGraphLaplace(sim);
        end

    else
        L = @(u) u;
    end
    
    
    
    
    %% Calculate Graph Laplace
%     ide = speye(size(L));
%     lipschitz = normest(ide - my * L); % 0.67
%     tau = 2/lipschitz*0.5; %tau in ]0,2/L[
    tau = 0.001;
    
    %% Run Algorithm
    
    % Prox Function
    prox = @(tau_,lambda_,v_) min(max(v_ - tau_ * lambda_, 0), 1);
    
    % Create Figure and Energy
    if(var_plot)
        fig = figure('Name', 'Energy RNMF');
    end
    
    
    for i = 1:iter
        
        % Update S
        S = S(:);
        for j = 1:5
            D = X - W * H;
            DVecPlusS = -D(:) + S;
            LS = my*L(S);
            gradF = DVecPlusS - LS;
            S = prox(tau,lambda, S - tau * gradF);
        end
        S = reshape(S,size(X,1), size(X,2));
        
        
        % Update W & H
        S_X = S - X;
        W( W < 1e-10 ) = 1e-10;
        H( H < 1e-10 ) = 1e-10;
        W =(( abs(S_X * H') - (S_X * H') ) ./ (2 * (W * (H * H')))) .* W;
        H =(( abs(W' * S_X) - (W' * S_X) ) ./ (2 * (W' * W * H)) ) .* H;
        
        % Norm
        tmp = sqrt(sum(W.^2));
        if sum(tmp) > 0
            W = W ./ tmp;
            H = H .* tmp';
        end
        
        
        
        % Energy
        SLong = (S(:));
        SLS = SLong' * L(SLong);
        Energy = 0.5 * norm(X-W*H-S,'fro').^2+ lambda*sum(sum(abs(S))) - my/2 * SLS;
        E = [E, Energy];
        
        % Plot
        if var_plot == true
            
            set(0, 'currentfigure', fig);
            subplot(1,2,1)
            plot(E);
            
            
            subplot(1,2,2);
            S_1 = reshape(S(:,1), size(X3DTemp,1), size(X3DTemp,2));
            imagesc(S_1);colorbar;
            title('S');
            
            drawnow;
        end
        
    end
    
    
    SFinal = cat(3, SFinal, reshapeas(S,X3DTemp));

end

% Full W and H
X = reshape(X3D, size(X3D,1) * size(X3D,2), size(X3D,3));
[W,H] = nnmf(X,r);
    
% Return Values
ret = struct('W',W,'H',H,'S',SFinal,'Energy', E);

end

