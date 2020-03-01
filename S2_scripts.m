%%Embryo edge detection from images slices perlendicular to animal-vegital axis
% Matlab 2019b

tic
clear;

%% parameters
% Image location
[imFile,imDir]= uigetfile('*.btf','Locate median filtered zx stack');
imPath        = fullfile(imDir,imFile);

% Image info
frame_first    = 2;     % 1st frame
frame_last     = 2;     % last frame (Shoud not exceed frame * zx slices < 65536)

pix            = 3;     % x range to detect border, pixel
qix            = 3;     % z range to detect border, pixel

imType         = 2;     % 0: dual side view
                        % 1: right side view
                        % 2: left side view

plotShow       = 0;     % 1: show plot
                        % 0: do not show
%% job submission
% Obtain image info
t              = Tiff(imPath,'r');% Tiff object handle
zxHight        = getTag(t,'ImageLength');
zxWidth        = getTag(t,'ImageWidth');
[nImages,~,nSlices,nFrames] = getHyperstackInfo(t);
close(t);   

% Indices of slice of interest
nSliceIdx      = (1:nSlices)';

% Border detection
BWc = false(zxHight, zxWidth, nSlices, nFrames); % Empty border matrix
p = floor(pix ./ 2); q = floor(qix ./ 2);

% Border for each frame
parfor iframes = frame_first:frame_last
    % Obtain index for specific frame 
    nFrameVec = repmat(iframes,[nSlices,1]);% subscript vector
    ind_f = sub2ind([nSlices,nFrames],nSliceIdx,nFrameVec); %indecis of frame f

    % Load whole images of frame 'iframe'
    I = zeros(zxHight,zxWidth,nSlices);% Image matrix
    for iSlices = 1:nSlices  % Load all slices in frame f
        I(:,:,iSlices) = imread(imPath,ind_f(iSlices)); 
    end
    A = I;

    % Obtain border position
    y = ones(zxWidth,nSlices);% empty file of the coordinate(x,z) at the border
    BW = false(size(I));% Logocal matrix for border position, same size as image
    edgetop = cell(nSlices,2);
    
    %% Obtain coordinate at border at a specific slice
     for iSlices2 = (q + 1):(nSlices - q) % Remove edge         
        %% Border detection
            for iX = (p+1):(zxWidth-p) % x value
                % Obtain border coordinate corresponds to each x value
                Is = I(:,(iX-p):(iX+p),(iSlices2-q):(iSlices2+q));
                Is2 = reshape(Is,zxHight,[]);
                Ip = permute(Is2,[2,1]);% Transpose xy coordinate for using ifindchangepts function
                yBorder = findchangepts(Ip);% Obtain signal change points along y axis
                y(iX,iSlices2) =yBorder;           
            end
            yz = y(:,iSlices2);
       %% Obtain binary image representing border
        ind = sub2ind([zxHight,zxWidth], y(:,iSlices2), [1:zxWidth]');% Change border position to index
        BWz = BW(:,:,iSlices2);
        BWz(ind) = 1; % Make border as 1
        if plotShow == 1
        subplot(4,2,1);
        imshow(BWz); title('findchangepts');
        end
        
       %% Post-processing
        % 1.Remove edges
        BWz(:,1:(p+1)) = 0;  BWz(:,(end-p):end) = 0;
        BWz(1,:) = 0;  BWz(end,:) = 0;

        % 2.Remouve unfitted region at the edge
        
        % 2.2 Connect uncontinuous points w/ morphology closing
        se = strel('disk', 10); % Structuring Elements, circle 10 pixel
        BWzmin = imclose(BWz,se); % Morphology closing
        if plotShow == 1
        subplot(4,2,3);
        imshow(BWzmin);title('Morphology closing');
        end
        
        % 2.3 Obtain position at local minimum
         [yz22,x] = find(BWzmin > 0);               
         % Soothing border
         iMedian = 0;
         while iMedian < 0
            yz22 = medfilt1(yz22,20);  
            iMedian = iMedian + 1;
         end
         yMed = medfilt1(yz22,20);
        
         % Search for local minimum
        [ylocalmin, prominence] = islocalmin(yMed,'MinProminence',0); 
        XMin = x(ylocalmin);
        
        % Replace local minimum to 0
        % 2.1 Extract elements smaller that speficied size from smoothed image
        BWz = bwareaopen(BWz, 20);
        if plotShow == 1
        subplot(4,2,2);
        imshow(BWz); title('Extract elements smaller than 20 pix');
        end
        % 2.2 Connect uncontinuous points, morphology closing
        se = strel('disk', 10); % structual elements, circle, 3 pixel
        BWz = imclose(BWz,se); % morpholgy closing
        
        if plotShow == 1
        subplot(4,2,3);
        imshow(BWz);title('Morphology closing');
        end
        
        % Replace the value at local minimum with 0
        if ~isempty(XMin)
            BWz(:,XMin) = 0;
        end       
        
        if plotShow == 1
        subplot(4,2,4);
        imshow(BWz); title('Find local minimum and replace its value with 0');
        end
                                     
        % 2.4 Check for existence of Gap
        [yGap,xGap] = find(BWz > 0);
        dyGap = zeros(size(yGap));
        dyGap(1:(end -1)) = diff(yGap);
        
        if plotShow == 1
        subplot(4,2,5);
        plot(xGap,dyGap);
        end
        
        %Position of Gap
        Gap1a = xGap(dyGap < -20);
        if ~isempty(Gap1a) 
            if size(Gap1a,1) > 1 
                gap1aIdx = find(dyGap < -20);
               [Gap1c, gap1Idx] = min(dyGap);% min dyGap                      
                xGap1 = xGap(gap1Idx); % x with min dyGap
            else
                xGap1 = Gap1a;
            end
        else
            xGap1 = [];
        end

        % Right edge of the Gap
        Gap2a = xGap((dyGap > 20));
        if ~isempty(Gap2a) 
            if size(Gap2a,1) > 1
                gap2aIdx = find(dyGap > 20);
                [Gap2c, gap2Idx] = max(dyGap);% max dyGap
                xGap2 = xGap(gap2Idx); % x with max dyGap
            else
                xGap2 = Gap2a;
            end
        else
            xGap2 = [];
        end   
        
        if ~isempty(xGap1) && ~isempty(xGap2) && xGap1 < xGap2
           BWz(:,xGap1:xGap2) = 0; % Mask between the two points
            numComp = 2;        
        else
            numComp = 1;
        end
        
       %% 2.4 Select region close to the embryo
        % Select objects smaller than 20 pixel
        BWz = bwareaopen(BWz, 20);
               
        % Label
        CC = bwconncomp(BWz);       
        L = labelmatrix(CC);
        nObj =CC.NumObjects;
        if nObj < 1
            BWfinal = false(zxHight,zxWidth);
        else
            if plotShow == 1
            subplot(4,2,6);
            imshow(label2rgb(L));title('label components');
            end

            % criteria 2: centroid closest to x axis midpoint
            % Get x value of centroid     
            cen = regionprops(CC,'Centroid'); % getting the edge position list of the labeled region
            cen_array = [cen.Centroid];
            cen_array2 = reshape(cen_array,2,[]);
            xcen = cen_array2(1,:);% x coordinates of centroid  
            
            % Object closest to midpoint
            imagexmid = round(zxWidth/2);% x of center point
            distance_midPos = abs(imagexmid - xcen);% list of distances from x of mid point
            [minD1, idxD1] = min(distance_midPos);% Closest object
            % 2nd object closest to midpoint
            distance_midPos2 = distance_midPos;
            distance_midPos2(idxD1) = imagexmid;
            [minD2, idxD2] = min(distance_midPos2);% 2nd object closest to midpoint   
            % 3rd object closest to midpoint
            distance_midPos3 = distance_midPos2;
            distance_midPos3(idxD2) = imagexmid;
            [minD3, idxD3] = min(distance_midPos3);% 3rd object closest to midpoint 

            % Selected object, x300:400
            if imType == 0 % dual side view
                idxEdge = [];
            elseif imType == 1 % right side view
                idxEdge = [1,2];
            elseif imType == 2
                idxEdge = nObj; % left side view
            end
                
            if numComp == 1 % If no gap, select closest to the center
                if xcen(idxD1) > 300 && xcen(idxD1) < 400
                     idx = cat(1,idxD1,idxEdge);
                     %idx = idxD1;
                else
                     idx = cat(1,idxD1,idxD2,idxEdge); 
                end

            elseif numComp == 2 % if a gap exists, select 2 objects closest to x midpoint, gap region masked
                 %idx = cat(1,idx1,idx2);   
                 idx = cat(1,idxD1,idxD2,idxEdge); 
            end  
            
            BWz = ismember(labelmatrix(CC),idx);
            
            if plotShow == 1
            subplot(4,2,7);
            imshow(BWz); title('Selected components');
            end

            % Mark border at z plane
            [edgetop{iSlices2,1},edgetop{iSlices2,2}] = edgeBorderImageNofig2(I(:,:,iSlices2));
            if     imType == 0                    % dual side view
                BWz(1:5,edgetop{iSlices2,1}) = 1; % left-upper edge
                BWz(1:5,edgetop{iSlices2,2}) = 1; % right-upper edge
            elseif imType == 1                    % right side view
                BWz(1:5,edgetop{iSlices2,2}) = 1; % right-upper edge
            elseif imType == 2                    % left side view
                BWz(1:5,edgetop{iSlices2,1}) = 1; % left-upper edge
            end        
            
           %% Smoothing spline��curve fit
            [yBWz,xBWz] = find(BWz == 1);
            % Smooth spline fit
            [xData, yData] = prepareCurveData( xBWz, yBWz );

            % Fitting type
            ft = fittype( 'smoothingspline' );

            % Model
            if size(xData,1) > 10
                [fitresult, gof] = fit( xData, yData, ft, 'SmoothingParam', 0.0001);
                y2 = fitresult(1:zxWidth);

              %% Obtain curve
                % Remove y with minus values
                yr = round(y2);
                yp2 = yr(yr > 0);

                % Remove y value beyond image
                yp = zeros(size(yp2,1),1);
                for i = 1:size(yp2,1)
                    if yp2(i,1) > zxHight
                        yp(i,1) = zxHight;
                    else
                        yp(i,1) = yp2(i,1);
                    end
                end
                xp = find(yr > 0);
                new_BW = zeros(zxHight,zxWidth);
                ind = sub2ind(size(BWz),yp, xp);
                new_BW(ind) = 1;
            else
                new_BW = zeros(zxHight,zxWidth);
            end

           %% Remove signal at the edge
            new_BW(1,:) = 0;  
            new_BW(end,:) = 1;   
            new_BW(:,1) = 0;  
            new_BW(:,end) = 0; 
            BWfinal = new_BW; 
        end
      
        if plotShow == 1
        subplot(4,2,8);
        imshowpair(A(:,:,iSlices2),BWfinal);title('curvefitted');
        end

        % 2D data -> 3D matrix
        BW(:,:,iSlices2) = BWfinal; 

    end   
    %% 4D border file
    BWc(:,:,:,iframes) = BW;
end

%Remove noise at xy planes
BWxyzf2 = noiseremoval(BWc);

% save
file_s = mfilename('fullpath');
[~,name_s,~] = fileparts(file_s);
if plotShow == 0
    save([name_s '_f' num2str(frame_first) '_' num2str(frame_last) '_x' num2str(pix) '_z_' num2str(qix)],'BWc','BWxyzf2','nSlices','nFrames','-v7.3');
else
    save([name_s '_f' num2str(frame_first) '_' num2str(frame_last) '_x' num2str(pix) '_z_' num2str(qix),'.fig']);
end
toc

function [BWxyzf2] = noiseremoval(BWc)

BWxyzf = permute(BWc, [3 2 1 4]);% Change axis(y->z, x->x, z->y, f->f)

% Make subscript vector to specify a xy plane
nY = size(BWxyzf,1);
nX = size(BWxyzf,2);
nSlice = size(BWxyzf,3); % total number of slices
nFrame = size(BWxyzf,4); % total number of frames
nImage = nSlice * nFrame; % total number of xy images
 
BWxyi = reshape(BWxyzf,[nY,nX,nImage]);
BW = zeros(size(BWxyi));

% Process image
for iImage = 1:nImage
    BWxy = BWxyi(:,:,iImage);        
    % Morphology closing to connect border
    sesize1 = 1;
    se = strel('disk',sesize1);
    closeBWxy = imclose(BWxy,se);     
    closeBWxy = imclose(closeBWxy,se);  
        
    % Remove small particles
    sesize2 = 3;
    closeBWxy2 = bwareaopen(closeBWxy, sesize2);
            
    % Remove horizontal small particles
    sesize3 = 1;
    se2 = strel('line',sesize3,0);
    closeBWxy3 = imopen(closeBWxy2, se2);
    
    BW(:,:,iImage)=closeBWxy3;
end

BWxyzf2 = reshape(BW,[nY,nX,nSlice, nFrame]);

end

function [x1,x2] = edgeBorderImageNofig2(I)
%% Detect border at edge
Imean = mean(I(1:3,:),1);
y = medfilt1(Imean,5); % Upper edge, median filtered
[ylocmax, xlocmax] = findpeaks(y); % Local max
midY1 = round((min(y) + ylocmax(1))/2);% half value of 1st local max

% left upper edge
x1 =find(y>midY1, 1 );% First x above half max
y1 = y(x1);
if y1 < 20
    x1 = [];
end

% right upper edge
midYlast = round((min(y) + ylocmax(end))/2);% Last positon of last half max
k1 = find(y< midYlast);% y: index of x corresponds to x below half max
x2 = min(k1(k1 > xlocmax(end)));
y2 = y(x2);
if y2 < 20
    x2 = [];
end
end
