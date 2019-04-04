function ret = robustNMF_TVW(X3D,vars)
% || X - W * H - S||^2 +    lambda1 * ||S||_1 + 
%                           lambda2 * TVxy(W) 



 
var_plot =  vars(1); % Plot
r =         vars(2); % Rank
iter =      vars(3); % Iterations
lambda1 =   vars(4); % Sparsity
lambda2 =   vars(5); % TVxy 

lambda3 = 1;
%% Variables


% Reshape Video
X = reshape(X3D,[size(X3D,1) * size(X3D,2), size(X3D,3)]);

% Initialize W,H,S,M
[W,H] = nnmf(X,r);
S = max(X - W * H,0);
Dxy = getDerivativeXY(size(X3D,1), size(X3D,2), 1);
DxyL = getDerivativeXY(size(X3D,1), size(X3D,2), r);


% Stepsize D
sigmaW = 0.8 / normest(DxyL);
tauW =  sigmaW;

% Energy and Plot
E = [];
if var_plot == true
    fig = figure('Name','RNMF with TV(W)');
end


% Energy
Energy = @(X,W,H,S)0.5 * norm(X - W * H - S,'fro')^2 + ...
                        lambda1 * sum(abs(S(:))) + ...
                        lambda2 * sum(abs(DxyL * W(:)));


%% Run Algorithm
for i = 1:iter
    
    % Update S
    S = X - W * H;
    ind = S > lambda1;
    S( ind ) = S( ind ) - lambda1;
    S( ~ind ) = 0;
    

    
    H = H';
    XS = X - S;
    
    % Update W
    for col = 1:r
        u = W(:,col);
        hj = H(:,col);
        
        pW1 = 0; 
        u_strich = u;
        
        HX = H;
        WX = W;
        HX(:,col) = 0;
        WX(:,col) = 0;
        
        Xj = 0;
        for it = 1:r
            Xj = Xj + WX(:,it) * HX(:,it)';
        end
        Xj = XS - Xj;
        
        for j = 1:15
            pW1 = min(max(pW1 + sigmaW * Dxy * u_strich, -lambda2), lambda2);
            
            u_old = u;
            u = (u - tauW * Dxy' * pW1 + lambda3 * tauW * Xj * hj) / (1 + lambda3 * tauW * hj' * hj);
            u = min(max(u,0),1);
            u_strich = 2 * u - u_old;
              
        end
        W(:,col) = u;
    end
    
    H = H';

    

    
    
    % Update H 
    eps = 0.0001;
    H = H';
    XS = X - S;
    
    for col = 1:r
        wj = W(:,col);
        hj = H(:,col);

        WW = (W' * W);
        XW = (XS' * W);

        H(:,col) = 1/(wj' * wj) * max(eps, XW(:,col) - H * WW(:,col) + hj * wj' * wj);        
    end
    
    H = H';
    

    

    
    
    % Energy
    EnergyNew = Energy(X,W,H,S);
        
    E = [E,EnergyNew];
        
    % Plot
    if var_plot == true
        set(0, 'currentfigure', fig)

        subplot(2,2,1)
        plot(E);
        title('Energy');
        
        subplot(2,2,2);
        W_1 = reshape(W(:,1),size(X3D,1), size(X3D,2));
        imagesc(W_1);colorbar;
        title('W_1');

        subplot(2,2,3);
        S_5 = reshape(S(:,5), size(X3D,1), size(X3D,2));
        imagesc(S_5);colorbar;
        title('S_5');

        
        drawnow;
    end
    
end



%% Return Values


% Create Return
ret = struct('W',W,'H',H,'S',reshapeas(S,X3D),'Energy', E);

end





