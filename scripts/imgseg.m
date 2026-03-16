function output = imgseg(imageName,threshold,pixel2um,saveName)
    %----------------------------------------------------------------------
    % Input:
    %     - imageName: the name of micrograph you want to analyze
    %     - threshold: the threshold used for binarize the image, e.g. 200/256
    %     - pixel2um: convert pixel to um. unit is um per pixel
    %     - saveName: the name of this image that will be used when saving the results
    %
    % Output:
    %     - output: ['No.', 'Centroid-x', 'Centroid-y', 'Area(pixel2)', 'Area(um2)', 'EQAD(um)']
    %----------------------------------------------------------------------
    

    %% Read the image and convert it to grayscale
    image = imread(imageName);
    image = rgb2gray(image);
    bw = imbinarize(image, threshold); 
    bw = imfill(bw, 'holes');
    figure
    imshow(bw)
    exportgraphics(gcf,join([saveName,'_bw.jpg']))
    
    
    %% Image segmentation. Find the pixels that belong to each blob
    [bwlabeled, numberOfBlobs] = bwlabel(bw, 8); 
    props = regionprops(bwlabeled, 'all');
    blobPixels = cell(numberOfBlobs, 1);
    [brow, bcol] = find(bwlabeled==0);
    boundaryPixels = [brow bcol];
    for k = 1:numberOfBlobs
        [row, col] = find(bwlabeled==k);
        blobPixels{k} = [row col ones(length(row),1)]; % [y, x, weight]
    end

    figure
    coloredLabels = label2rgb(bwlabeled, 'hsv', 'k', 'shuffle'); % pseudo random color labels
    figure
    imshow(coloredLabels);
    exportgraphics(gcf,join([saveName,'_blob_withGB.jpg']),'Resolution',300)
    
    
    %% Consider the boundary pixels
    while ~isempty(boundaryPixels)
        innerBoundaryPixels = [];
        oneNeighbourBoundaryPixels = []; % boundary pixels that have at leat one neighbor blobs
        for i = 1:size(boundaryPixels,1)
            p = boundaryPixels(i,:);
            allNeighbors = getNeighbors(p(1),p(2), bw);
    
            % Find the blob that the neighbors belong to
            neighborInBlob = allNeighbors(allNeighbors(:,3)>0,:);
            idx = sub2ind(size(bwlabeled),neighborInBlob(:,1), neighborInBlob(:,2));
            neighborBlobs = bwlabeled(idx);
            neighborBlobs = unique(neighborBlobs);
            if isempty(neighborBlobs)
                innerBoundaryPixels = [innerBoundaryPixels; p];            
            elseif length(neighborBlobs) == 1
                % Assign it to the unique neighbor blob
                weight = 1;
                oneNeighbourBoundaryPixels = [oneNeighbourBoundaryPixels; [p(1) p(2) neighborBlobs]];
                blobPixels{neighborBlobs} = [blobPixels{neighborBlobs}; [p(1) p(2) weight]];
            else
                % Assign it to any of the neighbor blobs randomly
                randBlob = randsample(neighborBlobs,1);
                weight = 1;
                oneNeighbourBoundaryPixels = [oneNeighbourBoundaryPixels; [p(1) p(2) randBlob]];
                blobPixels{randBlob} = [blobPixels{randBlob}; [p(1) p(2) weight]];
            end
        end
    
        % For boundary pixels that only have one neighboring blobs, modify bw and bwlabeled 
        idx = sub2ind(size(bw), oneNeighbourBoundaryPixels(:,1), oneNeighbourBoundaryPixels(:,2));
        bwlabeled(idx) = oneNeighbourBoundaryPixels(:,3);
        bw(idx) = 1;
    
        % Update boundaryPixels for next loop
        boundaryPixels = innerBoundaryPixels;
    
    end
      
    
    %% Plot the blobs
    coloredLabels = label2rgb(bwlabeled, 'hsv', 'k', 'shuffle'); % pseudo random color labels
    figure
    imshow(coloredLabels);
    exportgraphics(gcf,join([saveName,'_blob.jpg']),'Resolution',300)

    
    %% Get the area of each blob. 
    % Since the weight for every pixel is 1, we just need to count the number of pixels in each blob
    output = zeros(numberOfBlobs, 6);
    for k = 1 : numberOfBlobs           
	    blobAreaPixels = size(blobPixels{k},1);
        blobAreaMicrons = blobAreaPixels * pixel2um^2;
        blobEQAD = sqrt(blobAreaMicrons * 4/pi);
	    blobCentroid = props(k).Centroid;		
        output(k, :) = [k, blobCentroid(1), blobCentroid(2), blobAreaPixels, blobAreaMicrons, blobEQAD];
    
        % add labels for each blob
        hold on
        plot(blobCentroid(1), blobCentroid(2), 'k.')
        text(blobCentroid(1), blobCentroid(2), string(k), 'FontSize', 6)
    end
    hold off
    exportgraphics(gcf,join([saveName,'_blob_labeled.jpg']),'Resolution',300)
    savetxt(join([saveName,'_grains_all.txt']), output)

    % Get the inside and edge grains
    nrow = size(bwlabeled, 1);
    ncol = size(bwlabeled, 2);
    top = bwlabeled(1,:);
    bottom = bwlabeled(nrow,:);
    left = bwlabeled(:,1);
    right = bwlabeled(:, ncol);
    labelAtBorder = [top(:); bottom(:); left(:); right(:)];
    labelAtBorder = unique(labelAtBorder);

    edgeGrains = output(ismember(output(:,1), labelAtBorder), :);
    savetxt(join([saveName,'_grains_edge.txt']), edgeGrains)
    insideGrains = output(~ismember(output(:,1), labelAtBorder), :);
    savetxt(join([saveName,'_grains_inside.txt']), insideGrains)
       
    % Print some information about the grains
    [rows, columns, numberOfColorChannels] = size(bw);
    fprintf("The size of the image: %d %d \n", rows, columns);
    fprintf("The total area of the image in pixels: %d \n", rows*columns);
    fprintf("The total area of the grains in pixels: %d \n", sum(output(:,4)));
    fprintf('Total number of inside grains: %d \n', size(insideGrains,1))
    fprintf('Total area of inside grains: %d \n', sum(insideGrains(:,4)))
    fprintf('Total number of edge grains: %d \n', size(edgeGrains,1))
    fprintf('Total area of edge grains: %d \n', sum(edgeGrains(:,4)))

    
    
    %% functions    
    function A = getNeighbors(r, c, im)
        % ---------------------------------------------------------------------
        % Get the neighbouring pixels of a boundary pixel, only consider 4, i.e.
        % top, bottom, left and right. Special case for the corner or boundary pixels
        % Input: r, c - row and column index of the pixel
        %        im - the binary image
        % Output: A - a matrix of [row, column, pixelvalue], maximum 4 rows
        % ---------------------------------------------------------------------
        neighbors = [r c-1; r c+1; r-1 c; r+1 c];
        neighbors = neighbors((neighbors(:,1)>=1) & (neighbors(:,1)<=size(im,1)) & ...
                              (neighbors(:,2)>=1) & (neighbors(:,2)<=size(im,2)) ,:);
        ind = sub2ind(size(im), neighbors(:,1), neighbors(:,2)); % convert to linear indexing
        values = im(ind);
        A = [neighbors, values];    
    end
    
    
    function savetxt(fname, data)
        % Save the data        
        fid=fopen(fname,'w');
        header = ['No.', '\t', 'Centroid-x', '\t', 'Centroid-y', '\t',...
            'Area(pixel2)', '\t', 'Area(um2)','\t', 'EQAD(um)', '\r\n'];
        fprintf(fid, join(header));
        fprintf(fid, '%4d   %8.2f   %8.2f  %8d  %8.2f %8.2f \r\n', data'); % the transpose ' after data is important
        fclose(fid);
    end


end


