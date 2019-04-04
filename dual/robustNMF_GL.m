function ret = robustNMF_GL(X3D,vars)

%% Get Parameters

var_plot = vars(1);
r = vars(2);
iter = vars(3);
lambda = vars(4);
eps = vars(5);



%% Run Algorithm


% Reshape Video
X = reshape(X3D,[size(X3D,1) * size(X3D,2), size(X3D,3)]);


% Precalculate W and H
[W,H] = nnmf(X,r);
E = [];
S = max(X-W*H,0);

[~,Dx,Dy,~] = getDerivativeXYZ(size(X3D,1),size(X3D,2),size(X3D,3));
D = cat(1,Dx,Dy);
Energy =@(W,H,S) norm((X-W*H).*(1-S),'fro').^2+ lambda * norm(S,1) ...
                 + 1/eps * sum(sum( 0.25 * ((2*S-1).^2-1).^2)) ...
                 + eps/2 * norm(D * S(:),'fro')^2;


if(var_plot)
    fig = figure('Name', 'Energy RNMF');
end

tau = 0.0001;%0.0001;

for i = 1:iter
    
    % Calc S
    R = (X - W * H); R = R(:);
    s = S(:);
    for j = 1:5
        Deriv1 = 1/eps * 2 * (( 2 * s - 1).^2 - 1) .* (4 * (2 * s - 1));
        Deriv2 = eps * D' * (D * s);
        Deriv3 = -2 * (R - s .* R) .* R;
        
        s = s - tau * (Deriv1+Deriv2+Deriv3);
    end
    S = reshapeas(s,S);


    % Update W and H 
    epsi = 0.0001;
    P = (1 - S);
    
    % Update H 
    for col = 1:r
        wj = W(:,col);
        
        % Create Xj
        Xj = 0;
        for it = 1:r
            if it ~= col
                Xj = Xj + W(:,it) * H(it,:);
            end
        end
        Xj = X - Xj;
        top = (wj' * (P .* P .* Xj));
        bottom =  max(epsi, (wj' * (P .* P .* (wj * ones([1,size(P,2)]) ))));
        H(col,:) = top ./  bottom;
    end
    
    % Update W
    for col = 1:r
        hj = H(col,:);
        
        % Create Xj
        Xj = 0;
        for it = 1:r
            if it ~= col
                Xj = Xj + W(:,it) * H(it,:);
            end
        end
        Xj = X - Xj;
        
        W(:,col) = ((P .* P .* Xj) * hj') ./  max(epsi, ((P .* P .* ( ones([1,size(P,1)])' * hj ) ) * hj' ));
    end
    
    
    
    

	% Energy
    E = [E,Energy(W,H,S)];
    
    % Plot
    if var_plot == true
        set(0, 'currentfigure', fig)

        subplot(2,2,1)
        plot(E);

        subplot(2,2,2);
        W_1 = reshape(W(:,1),size(X3D,1), size(X3D,2));
        imagesc(W_1);colorbar;
        title('W1');

        subplot(2,2,4);
        S_5 = reshape(S(:,5), size(X3D,1), size(X3D,2));
        imagesc(S_5);colorbar;
        title('S');
        
        drawnow;
    end
    
end
%% Binarize Mask
MRes = S;
MRes(S > 0.5) = 1;
MRes(S <= 0.5) = 0;
MR = reshape(MRes, size(S,1), size(S,2));
MR = reshape(MR,size(X3D,1), size(X3D,2), size(X3D,3)); 
MR = max(min(MR,1),0);
%% Return Values


ret = struct('W',W,'H',H,'S',reshapeas(S,X3D),'M',MR,'Energy', E);


