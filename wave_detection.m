%% Code for wave surface detection for "High-resolution PIV measurement over wind-generated waves" paper

% the code is developed by Wagih Abu Rowin and Kevin

% date May 2023


clearvars; close all; clc

% Load data and get a list of files
load('/Volumes/Kevin/AirSea/EiF figures/Attached code/PCO4000_DI.mat')
lm_dirlist = dir(['/Volumes/Kevin/AirSea/EiF figures/Attached code/' '*Cam_' '*.b16']);

for img = 1 : length(lm_dirlist)
    disp(['Loading image #' num2str(img)])
    
    % Load image data
    clearvars -except c img lm_dirlist
    fid = fopen(lm_dirlist(img).name, 'r');
    fseek(fid, c.HeaderLength, 'bof');
    imagea = ctranspose(fread(fid, [c.ImageWidth c.ImageHeight/2], 'uint16=>double'));
    imageb = ctranspose(fread(fid, [c.ImageWidth c.ImageHeight/2], 'uint16=>double'));
    fclose(fid);
    
    % Preprocess images
    im = (imagea + imageb) / 2;
    filt = fspecial('gauss', [5 5], 2);
    imf = filter2(filt, im, 'same');
    im_hp = max(im - imf, 0);
    im_hp2 = im_hp.^2;

    % Create a binary mask to remove certain elements
    temp = zeros(size(im));
    temp(im > 1000) = 1;
    CC = bwconncomp(temp);
    for nn = 1:CC.NumObjects
        id = cell2mat(CC.PixelIdxList(nn));
        if numel(id) > 625
            im_hp2(id) = 0;
        end
    end
    
    % Apply a low-pass filter
    im_hp2_lp = filter2(fspecial('average', [10 48]), im_hp2, 'same');
    im_hp2_lp(end-7:end, :) = 0;

    % Detect intensity jump for each column of pixels
    th = mean(im_hp2_lp(2500:end-7, :), 1);
    th = repmat(th, [2672 1]);
    d1 = im_hp2_lp - 1.3 * th;
    temp = d1 > 0;

    % Apply smoothing and find surface
    d2 = max(d1, 0);
    d3 = d2 * 2;
    temp3 = zeros(size(im_hp2_lp));
    start_scan = 2200;
    end_scan = 2500;
    
    cc = 51:4008;
    for n = cc
        im_profile = d3(start_scan:end_scan, n);
        im_sm = smooth(im_profile, 10);
        [Minima, MinIdx] = findpeaks(-im_sm);
        im_diff = diff(Minima);
        
        if isempty(im_diff)
            th = abs(max(im_sm));
        else
            [F, I] = max(im_diff);
            th = abs(Minima(I, 1));
        end
        
        temp3(d3(:,n) > th, n) = 1;
    end
    
    temp3(1, :) = ones([1 size(temp3, 2)]);
    
    % Find surface
    cc = 51:4008;
    rr_3 = arrayfun(@(n) find(temp3(:, n), 1, 'last'), cc);

    temp(1, :) = ones([1 size(temp, 2)]);
    
    cc = 51:4008;
    rr_1 = arrayfun(@(n) find(temp(:, n), 1, 'last'), cc);
    
    rr = (rr_3 + rr_1) / 2;
    
    % Apply Hampel method to remove outliers
    rr = hampel(rr, 400, 1.5);
    id = find(isnan(rr));
    cc(id) = [];
    rr(id) = [];
    
    % Fit piecewise polynomials
    fitresult2 = fit(cc', rr', 'smoothingspline', 'SmoothingParam', 5e-8);
    col = 51:4008;
    row = feval(fitresult2, col);
    
    % Plot the result
    filt = fspecial('gaussian', [5 5], 0.6);
    imagea = filter2(filt, imagea, 'same');
    imageb = filter2(filt, imageb, 'same');
    
    close all
    figure()
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1])
    imagesc(imagea);
    clim([400 750]);
    colormap gray
    set(gca, 'ylim', [1500 2672], 'fontsize', 6);
    daspect([1 1 1])
    hold on
    
    plot(col, row, '--c', 'linewidth', 2)
    pause
end