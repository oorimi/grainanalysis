%% Fig.1f
imageName = 'example.jpg';   % name of the input image
threshold = 200 / 256;       % threshold used for binarize the image
pixel2um = 200 / 122.66;     % conversion of scale: um/pixel
saveName = 'example';        % name for saving the results
fprintf('Analyzing: %s \n', imageName)
output = imgseg(imageName,threshold,pixel2um,saveName);
