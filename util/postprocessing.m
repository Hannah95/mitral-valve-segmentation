function segm = postprocessing(mask,display)

    % Get Centroid
    sumImg = sum(mask, 3);

    [m,n]=size(sumImg);
    [xPos,yPos]=ndgrid(1:m,1:n);
    weightMat = ones(size(sumImg));
    middelYs = [];
    middelXs = [];
    for i = 1:10
        var = 100 - i*10;
        C = regionprops(true(size(weightMat.*sumImg)), weightMat.*sumImg,  'WeightedCentroid');
        middelY = round(C.WeightedCentroid(2));
        middelX = round(C.WeightedCentroid(1));
        middelYs = [middelYs; middelY];
        middelXs = [middelXs; middelX];
        
        weightMat = gauss(xPos,yPos, middelY, middelX, var, var); 
        
    end

%     Closing per Frame
    closedM = zeros(size(mask));
    se = strel('disk',1);
    for i = 1:size(mask,3)
       closedM(:,:,i) = imclose(mask(:,:,i),se);   
    end
% 	closedM = mask;

    % BWComp per Vid
    CC = bwconncomp(closedM);
    L = labelmatrix(CC);


    % Table with Distance to Centroid
    dist = @(A,B) sqrt(sum(bsxfun(@minus, B, A).^2,2));
    tableDist = Inf(size(L,3),CC.NumObjects);
    for frameNr = 1:size(L,3)   
        includeIdx = unique(L(:,:,frameNr));
        for classNr = 2:numel(includeIdx)
            idxClass = find(L(:,:,frameNr)==includeIdx(classNr));
            posClass =  [xPos(idxClass),yPos(idxClass)];
            min_distance = min(dist(posClass,[middelY,middelX]));
            tableDist(frameNr,includeIdx(classNr)) = min_distance;
        end
    end

    % Cut L
    cutL = L;
    for frameNr = 1:size(L,3)  
        [~,idx] = min(tableDist(frameNr,:));
        frame = cutL(:,:,frameNr);
        
        if ~isempty(idx)
            bed = (frame == idx);
            frame(bed) = 1;
            frame(~bed) = 0;
        else
            frame = zeros(size(frame));
        end
        
        cutL(:,:,frameNr) = frame;
    end

    segm = cutL;
    
    
%     img = label2rgb(L(:,:,33));
% RGB = insertMarker(img, [middelX middelY],'circle');
% figure;imshow(RGB);
% 
% RGB = double(RGB);
% RGB = RGB / 255;
    % play
    if display
        figure;
        for i = 1:size(L,3) 
            subplot(2,2,1);
            imagesc(cutL(:,:,i)); title('Final Segmentation');
            axis off;
            
            subplot(2,2,2);
            img = label2rgb(L(:,:,i));
            RGB = insertMarker(img, [middelX middelY],'circle');
            imshow(RGB);
            title(['Clusters ',num2str(i)]);
            axis off;
            
            subplot(2,2,3);
            S = mask;
            imagesc(0.7 * double(cutL(:,:,i)).* S(:,:,i) + 0.3 * S(:,:,i));
            title('Final Segmentation in Origial Video');
            axis off;
            
            subplot(2,2,4);
            plotIt(sumImg,middelYs,middelXs);
            title('Calculated Centroid');
            axis off;
            
            drawnow;
            waitforbuttonpress;
        end
    end

end

%% Functions


% BWComp
function res = gauss(x,y, myX, myY, sigmaX, sigmaY)

    pre = 1 / (2 * pi * sigmaX * sigmaY);
    expoX = 1/(2 * sigmaX^2) * (x - myX).^2;
    expoY = 1/(2 * sigmaY^2) * (y - myY).^2;
    res = pre * exp( - expoX - expoY );

end

function plotIt(sumImg,middelY,middelX)
color = 'black';
    RGB = insertMarker(normalizeVal(sumImg), [middelX middelY],'circle','Color',color);
    axis off;
    imagesc(1-RGB(:,:,1));
    set(gca,'xtick',[],'ytick',[])
end

function NX = normalizeVal(X)
NX = ( X - max(X(:))) / (min(X(:)) - max(X(:)));
end

function plotCol(sumImg,middelY,middelX)
blackImg =  zeros(size(sumImg));
alpha  = insertMarker(blackImg, [middelX middelY],'circle','Color','white');
alpha = rgb2gray(alpha);
alpha = imbinarize(alpha);

% Plot
figure;

imagesc(sumImg(:,:,1));
hold on;
h =imshow(blackImg);
set(h, 'AlphaData', alpha);
end

