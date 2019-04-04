function ret = robustNMF_DS(X3D,vars)
%E = 0.5 * norm(X - W * H - S,'fro').^2 + lambda1 * norm(S,1) + lambda2 * norm(D_XY * S,2) + lambda3 * norm(D_Z * S,2);


% Reshape Video
X = reshape(X3D,[size(X3D,1) * size(X3D,2), size(X3D,3)]);


% Input Variables
var_plot = vars(1); % Plot
r = vars(2);        % Rank
iter = vars(3);     % Iterations
lambda1 = vars(4);  % Sparsity
lambda2 = vars(5);  % D_XY(S)
lambda3 = vars(6);  % D_Z(S)




%% Variables

% Precalculate W and H
[W,H] = nnmf(X,r);
S = max(X - W * H,0);


% Derivative
[~,Dx,Dy,D_Z] = getDerivativeXYZ(size(X3D,1),size(X3D,2),size(X3D,3));
D_XY = cat(1,Dx,Dy);

% Step Size
tau = 0.1;

% Energy
EnergyAll = [];

% Prox Function
prox = @(tau_,lambda_,v_) min(max(v_ - tau_ * lambda_, 0), 1);

% Plot
if var_plot == true
    fig = figure('Name','RNMF with TV(S)');
end

%% Algorithm
for i = 1:iter
    
    % Vectorize X - W * H
    R = X - W * H; R = R(:);
    
    % Update S
    u = S(:);
    
    for j = 1:5
        
        % Derivative F
        DF = -R + u + lambda2 * 2 * (D_XY') * D_XY * u + lambda3 * 2 * (D_Z') * D_Z * u;
        u = prox(tau,lambda1,u - tau * DF);
        
    end
    
    % Reshape S back
    S = reshape(u, size(S,1), size(S,2));

    % Update W and H
    S_X = S - X;
    W( W < 1e-10 ) = 1e-10;
    H( H < 1e-10 ) = 1e-10;
    W =(( abs(S_X * H') - (S_X * H') ) ./(2 * (W * (H * H')))) .* W;
    H =(( abs(W' * S_X) - (W' * S_X) )./ (2 * (W' * W * H)) ) .* H;
    
    % Norm
    tmp = sqrt(sum(W.^2));
    if sum(tmp) > 0
        W = W ./ tmp;
        H = H .* tmp';
    end
    
    
    % Calculate Energy
    Energy = 0.5 * norm(X-W*H-S,'fro').^2 + ...
            lambda1 * sum(abs((S(:)))) + ...
            lambda2 * norm((D_XY * S(:)),2) + ...
            lambda3 * norm((D_Z * S(:)),2);
        
    EnergyAll = [EnergyAll,Energy];
    
    % Plot
    if var_plot == true
        set(0, 'currentfigure', fig)

        subplot(1,2,1)
        plot(EnergyAll);
        title('Whole Energy');
        
        subplot(1,2,2);
        S_5 = reshape(S(:,5), size(X3D,1), size(X3D,2));
        imagesc(S_5);colorbar;
        title('S');
        
        drawnow;
    end
end


%% Return Values


% Create Return
ret = struct('W',W,'H',H,'S',reshapeas(S,X3D),'Energy', EnergyAll);

end


%% Functions

function D = getDerivativeXY(n,m) 
% n = dim(x), m = dim(y)

D1 = speye(n);
D1 = D1 + sparse(1:n,max(1,0:n-1),-1,n,n);
D1 = kron(speye(m),D1);

D2 = speye(m);
D2 = D2 + sparse(max(1,0:m-1),1:m,-1,m,m);
D2 = kron(D2',speye(n));

D = [D1 ; D2];
end

function D = getDerivativeZ(n,m,z)
% n = dim(x), m = dim(y), z = dim(z)
n = n * m;
m = z;

D2 = speye(m);
D2 = D2 + sparse(max(1,0:m-1),1:m,-1,m,m);
D2 = kron(D2',speye(n));

D = D2;
end



