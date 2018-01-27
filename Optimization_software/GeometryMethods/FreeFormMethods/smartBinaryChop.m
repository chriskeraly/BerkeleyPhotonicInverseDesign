% clear; close all;

% tic;
% 
[numY, numX] = size(epsGrid0);

[xGrid, yGrid] = meshgrid(1:numX, 1:numY);

% epsGrid0 = ((xGrid - 400).^2 + (yGrid - 450).^2) < 200^2;
%epsGrid0 = epsGrid0 | ((xGrid - 700).^2 + (yGrid - 650).^2) < 200^2;

epsGrid = epsGrid0;

figure(1); imagesc(epsGrid); colormap(bluewhitered); axis equal;

% blockRad = .5*min(numX, numY);

cnt = 0;
yArr=[];
xArr=[];
blockSizeArr=[];

for blockSize = 40:-2:0
    
    filterSize = floor(blockSize/2)*2 + 1; % Must be ODD (re-enforce in case for loop is changed)
    filter = ones(filterSize);
    
    epsGridOld = zeros(size(epsGrid));
    
    while( any(any(epsGrid - epsGridOld)) )
        
        epsGridOld = epsGrid;
        
        epsGrid2 = conv2(1*epsGrid,filter,'same');% == filterSize^2;
        
        %         figure(2); imagesc(epsGrid2); colormap(bluewhitered); axis equal;
        
        epsGrid2 = epsGrid2 == filterSize^2;
        
        if(filterSize > 1)
            [yInd, xInd] = find(epsGrid2,1,'first');
            if(~isempty(yInd))
                cnt = cnt + 1;
                yArr(cnt) = yInd;
                xArr(cnt) = xInd;
                blockSizeArr(cnt) = filterSize;
                
                dEpsGrid = (abs(xGrid - xInd) <= filterSize/2) & (abs(yGrid - yInd) <= filterSize/2);
                
                %              figure(3); imagesc(~dEpsGrid); colormap(bluewhitered); axis equal;
                
                epsGrid = epsGridOld & ~dEpsGrid;
            end
        else
            [yInd, xInd] = find(epsGrid2);
            len = length(yInd);
            yArr(cnt+1:cnt+len) = yInd;
            xArr(cnt+1:cnt+len) = xInd;
            blockSizeArr(cnt+1:cnt+len) = filterSize;
            cnt = cnt + len;
        end
        %              figure(4); imagesc(epsGrid); colormap(bluewhitered); axis equal;
        
    end
end

dx = geo_x(2)-geo_x(1);
xArr = geo_x(1) + (xArr-1)*dx;
yArr = geo_y(1) + (yArr-1)*dx;
blockSizeArr = blockSizeArr*dx;
thickness = geo_z(2)-geo_z(1);
z0 = geo_z(1) + thickness/2;

% toc

% tic
% epsGridNew = zeros(size(epsGrid));
% for cnt = 1:length(yArr)
%     dEpsGrid = (abs(xGrid - xArr(cnt)) <= blockSizeArr(cnt)/2) & (abs(yGrid - yArr(cnt)) <= blockSizeArr(cnt)/2);
%     epsGridNew = epsGridNew + dEpsGrid;
% end

% toc
% figure(5); imagesc(epsGridNew); colormap(bluewhitered); axis equal;

