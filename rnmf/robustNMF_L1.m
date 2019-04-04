function ret = robustNMF_L1(X3D,vars)

%% Get Parameters

var_plot = vars(1);
r = vars(2);
iter = vars(3);
lambda = vars(4);




%% Run Algorithm

% Reshape Video
X = reshape(X3D,[size(X3D,1) * size(X3D,2), size(X3D,3)]);

% Precalculate W and H
[W,H] = nnmf(X,r);
S = max(X -  W * H,0);

p1 = 0;
p2 = 0;
u_strich = rand(size(S(:)));
sigma = 0.99;
tau = sigma;
E = [];

Energy = @(R,S,lambda) sum(sum(abs(R - S))) + lambda * sum(abs(S(:)));

if(var_plot)
    fig = figure('Name', 'Energy RNMF L1');
end

for i = 1:iter

    R = X -  W * H; 
    S(R < 0) = 0;
    
    if lambda < 1
        S(R > 0) = R(R > 0);
    end
    if lambda > 1
         S(R > 0) = 0;
    end
    if lambda == 1
        S(R > 0) = R(R > 0);
    end
        
        
    % Update W and H
    Xhat = (X-S);
    eps = 0.000000001;
    down = sqrt((X-W*H-S).^2 + eps^2);
    D = 1./down;
    W = W .* (Xhat .* D * H') ./ max((W * H) .* D * H',eps);
    H = H .* (W' * (Xhat .* D)) ./ max(W' * ((W * H) .* D),eps);
   


	% Energy
	ENew = Energy(X - W * H,S,lambda);
    E = [E, ENew];
    
    % Plot
    if var_plot == true
        set(0, 'currentfigure', fig)
        
        subplot(2,2,1);
        WR = reshapeas(W(:,1),X3D(:,:,1));
        imagesc(WR); colorbar; title('W1');

        subplot(2,2,2);
        SR = reshapeas(S,X3D);
        imagesc(SR(:,:,1)); colorbar; title('S');

        subplot(2,2,3);
        plot(E); title('Energy');
        
        drawnow;
    end
    
end

%% Return Values

ret = struct('W',W,'H',H,'S',reshapeas(S,X3D),'Energy', E);




function data = convConj_L21(lambda, data, split)

        p_temp = reshape(data, size(data,1) / split, split);
        p_temp_norm = sqrt(sum(p_temp.^2,1));
        beta = p_temp_norm > lambda;
        
        for j = 1:(size(data,1)/ split)
            p_temp(j, beta) = p_temp(j, beta) * lambda ./ p_temp_norm(beta);
        end
        data = p_temp(:);


