function ret = robustNMF_TVS(X3D,vars)
%E = 0.5 * norm(X - W * H - S,'fro').^2 + lambda2 * norm(D * S,1) + lambda1 * norm(S,1)


% Reshape Video
X = reshape(X3D,[size(X3D,1) * size(X3D,2), size(X3D,3)]);


% Input Variables
var_plot = vars(1); % Plot
r = vars(2);        % Rank
iter = vars(3);     % Iterations
lambda1 = vars(4);  % Sparsity
lambda2 = vars(5);  % TV(S)
time = vars(6);     % Derivative in time or spacial

%% Variables

% Precalculate W and H
[W,H] = nnmf(X,r);
S = max(X - W * H,0);


% Create Primal Parameter
u = S(:);
u_strich = u;


% Create Dual Parameter
p1 = 0;
p2 = 0;
% K2 = getDerivative(size(X3D,1),size(X3D,2));
% K2 = kron(speye(size(S,2)),K2);
% K2 = repmat(K2,1,size(S,2));
K1 = 1;%speye(size(u,1));

[~,Dx,Dy,Dz] = getDerivativeXYZ(size(X3D,1),size(X3D,2),size(X3D,3));
K2 = cat(1,Dx,Dy);
if time == 1
    K2 = Dz;
end
if time == 2
    K2 = cat(1,K2,Dz);
end

% Step Size
sigma = 0.2;
tau = sigma;

% Energy
EnergyAll = [];

if var_plot == true
    fig = figure('Name','RNMF with TV(S)');
end

%% Algorithm
for i = 1:iter
    
    % Vectorize X - W * H
    vecXWH = X - W * H; vecXWH = vecXWH(:);
    
    % Update S
    u = S(:);
    
    for j = 1:5
        % Dual
        p1 = min(max(p1 + sigma * K1 * u_strich, -lambda1), lambda1);
        p2 = min(max(p2 + sigma * K2 * u_strich, -lambda2), lambda2);

        % Primal
        u_old = u;
        u = (u - tau * K1' * p1 - tau * K2' * p2 + tau * vecXWH) / (1 + tau);
        u = min(max(u,0),1);
        u_strich = 2 * u - u_old;

       
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
    Energy = 0.5 * norm(X-W*H-S,'fro').^2 + lambda2 * sum(abs(K2 * S(:))) + lambda1 * sum(abs((K1 * S(:))));
    EnergyAll = [EnergyAll,Energy];
    
    % Plot
    if var_plot == true
        set(0, 'currentfigure', fig)

        subplot(1,2,1)
        plot(EnergyAll);
        title('Whole Energy');

        subplot(1,2,2);
        S_1 = reshape(S(:,1), size(X3D,1), size(X3D,2));
        imagesc(S_1);colorbar;
        title('S');
        
        drawnow;
    end
end


%% Return Values

ret = struct('W',W,'H',H,'S',reshapeas(S,X3D),'Energy', EnergyAll);




end






%% Functions

% function D = getDerivative(n,m)
% 
% D1 = speye(n);
% D1 = D1 + sparse(1:n,max(1,0:n-1),-1,n,n);
% D1 = kron(speye(m),D1);
% 
% D2 = speye(m);
% D2 = D2 + sparse(max(1,0:m-1),1:m,-1,m,m);
% D2 = kron(D2',speye(n));
% 
% D = [D1 ; D2];
% end


