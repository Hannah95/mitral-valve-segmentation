function ret = segmentMBO(vid, vars)

    plotFig =           vars(1);
    iterations =        vars(2);  
    alpha =             vars(3);  % Alpha of Sim-Matrix
    patchsizeXY =       vars(4); 
    patchsizeZ =        vars(5);
    nrSim =             vars(6);
    eta =               vars(7);  % Dataterm
    simType =           vars(8);%1,2,3
    splitSize =         vars(9);
    nrEV =              vars(10);

    sz = size(vid.S,3);
    numberSplit = ceil(sz/splitSize);
    
    M = [];
    E = [];
    for split = 1:numberSplit
        %% Load Video/Image
        fBeg = (split - 1) * splitSize + 1;
        fEnd = min(fBeg+splitSize - 1,sz);
        WH = reshapeas(vid.W * vid.H,vid.S);
        WH = WH(:,:,fBeg:fEnd);
        S = vid.S(:,:,fBeg:fEnd);
        
        %% Create Fidelity
        [gtAvailable,gt] = createFidelity(S,WH);
        
        %% Eigenvectors and Eigenvalues
        opt = {};
        opt.ensureRank = 1;
        opt.nrSim = nrSim;
        opt.patchsizeXY = patchsizeXY;
        opt.patchsizeZ = patchsizeZ;
        opt.method = 'dense';
        opt.alpha = alpha;
        opt.simType = simType;
        
        [simMatrix,opt] = getSimilarityMatrixNonlocalMeans(S, opt);
        XX = simMatrix.XX;
        XY = simMatrix.XY;
        [eigVal,eigVec] = nystroem(XX,XY,nrEV);
        eigVec(opt.permutation,:) = eigVec;

        
        %% MBO
        
        % Energies (used to plot)
        EnergyData = @(u) 1/2 * eta * sum(gtAvailable .* vecnorm(u - gt,2,2).^2);
        EnergyW =    @(u) 1/eps * sum(WFunc(u));
        EnergyL =    @(u)  eps * sum(dot(u, (eigVec * eigVal * (eigVec' * u))));
        Energy =     @(u) EnergyData(u) + EnergyL(u) + EnergyW(u);


        dt = 0.2;
        u = gt;
        
        % Hotone Function
        hotone = @(X)bsxfun(@eq, X(:), 1:max(X));
        
        for i = 1:iterations
            
            for j = 1:1
                a = eigVec' * u;
                d = eigVec' * eta * (gtAvailable .* (u - gt));
                a = (a - dt * d) ./ (1 + dt * diag(eigVal));
                u = eigVec * a;
            end
            
            %Binarize
            [~,I] = max(u,[],2);
            u = double(hotone(I));
            if(unique(I)) < size(gt,2)
                temp = zeros(size(gt));
                temp(:,1:size(u,2)) = u;
                u = temp;
            end
            
            % Plot
            E = [E,Energy(u)];
            if plotFig
                plot(E);
            end
            
            % Stopping Criteria
            if i > 1
                if abs(E(end) - E(end-1)) < 1e-9
                    break;
                end
            end
            
            
        end
        
        %% Reshape Mask
        [~,I] = max(u,[],2);
        I = I - 1;
        if  max(I(:)) > 0
            I = I / max(I(:));
        end
        M = cat(3,M,reshapeas(I,S));
    end
    
    % Output
    ret = struct('W',vid.W,'H',vid.H,'S',vid.S,'M', M, 'Energy', E);
end
%% Functions

function [eigValNew,eigVecNew] = plot_and_cut_EV(eigVec,eigVal, nrEV,Img,plotFig)

    nrX = 2;
    nrY = 3;
    NrEVPlot = 3;

    %Cut Number Eigenvalues
    eigValNew = eigVal(1:nrEV,1:nrEV);
    eigVecNew = eigVec(:,1:nrEV);

    % Plot Eigenvalues
    if plotFig
        figure('Name','Eigenvectors and Eigenvalues');
        subplot(nrX,nrY,1);
        % Orthogonality
        orth = eigVec' * eigVec;
        imagesc(orth);
        colorbar;
        title('Orthogonality All Eigenvectors');
        for i = 1:NrEVPlot
            subplot(nrX,nrY,i+2);
            evp = reshapeas(eigVecNew(:,i),Img);
            imagesc(evp(:,:,1));
            colorbar;
            title(['Eigenvector (frame 1) ',num2str(i)]);
        end
        subplot(nrX,nrY,NrEVPlot+2+1);
        plot(diag(eigValNew));
        title('Few Eigenvalues');
        
        % Orthogonality
        orth = eigVecNew' * eigVecNew;
        subplot(nrX,nrY,2);
        imagesc(orth);
        colorbar;
        title('Orthogonality Few Eigenvectors');
        
    end
    

end




function [gtAvailable,fb] = createFidelity(S,WH)


    S = S(:);
    WH = WH(:);
    foreground = double(S > 0.1);
    background = double(S < 0.1 );

    foreground = max(foreground - (WH > 0.1),0);

    gtAvailable = (foreground>0 | background > 0);

    foreground(~gtAvailable) = 0.5;
    background(~gtAvailable) = 0.5;

    fb = cat(2,background,foreground);


end


function mult = WFunc(u)
    e = zeros(size(u));
    mult = 1;
    for k = 1: size(u,2)
        ek = e;
        ek(:,k) = 1;
        
        r = 0.25 .* sum(abs(u-ek),2).^2;
        mult = mult .* r;
    end
end