classdef Grating
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %   assuming z to be thickness
    % move dx to external parameter section in all shape classes
    % Get rid of thetaOut as return arg from getPadOut
    
    % New Uncommitted changes
    % 1. added dz, thickness to dFOM
    
    properties
        %External Parameters
        x0;
        y0; % bottom left corner
        z0; % center
        dx; %dx = dy
        xGrid;
        yGrid;
        zVec;
        
        %Internal Mesh
        numX;
        numY;
        numZ;
        numDepths;
        zDepths;
        
        %Boundary Conditions 0=nothing, +-1=sym, 2=periodic
        xBC;
        yBC;
        
        %Freeform shape properties
        eps; %filled-in material
        epsOut; %background material
        epsGrid; %freeform eps matrix, rows correspond to zDepths
        VBoundary;
        VBoundaryAll; % cell array of size [numFreq,numUser]
        maskGrid;
        
        %Data = cell(numMon, numFreq, numUser) [numY, numX]
        E;
        EA;
        
        numFreq;
        numUser;
        
        % Current Shapes
        currShapePad;
        maxMove; % limits growth of shapes
        eraseSize;
        velPadding; % Size of velocity region
        radiusCurv;
        
        % New shapes
        newShapeCreation;
        newShapeType; % string
        newShapePad;
        newShapeRad; % limits size of newly created shapes
        velocity;
        velocityAll; % cell array of size [numFreq,numUser]
        
        % 2D Shape Objects
        geoXY;
        
    end
    
    methods
        function obj = Grating(x0, y0, z0, xLengthReal, yLengthReal, dx, zDepths, eps, epsOut, xBC, yBC, newShapeCreation, newShapeRad, newShapePad, maxMove, eraseSize, velPadding, currShapePad, radiusCurv, numFreq, numUser)
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            
            obj.numFreq = numFreq;
            obj.numUser = numUser;
            
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx));
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx));
            obj.numDepths = length(zDepths)-1;
            
            [obj.xGrid, obj.yGrid] = meshgrid(x0+dx*(0:obj.numX-1)+dx/2, y0+dx*(0:obj.numY-1)+dx/2);
            obj.maskGrid=ones(obj.numY,obj.numX);
            
            obj.dx = dx;
            obj.numZ = floor(1e-4*round(1e4*(zDepths(end) - zDepths(1)) / dx));
            obj.zDepths = floor(1e-4*round(1e4*(zDepths-zDepths(1)) / dx)) + 1;
            obj.zVec = z0+dx*(0:obj.numZ-1)+dx/2;
            
            obj.eps = eps;
            obj.epsOut = epsOut;
            obj.epsGrid = zeros(obj.numY,obj.numX,obj.numDepths);
            obj.VBoundary = zeros(obj.numY,obj.numX);
            obj.VBoundaryAll = {zeros(obj.numY,obj.numX)};
            obj.velocity = zeros(obj.numY,obj.numX);
            obj.velocityAll = {zeros(obj.numY,obj.numX)};
            
            obj.xBC = xBC;
            obj.yBC = yBC;
            obj.newShapeCreation = newShapeCreation;
            obj.maxMove = floor(1e-4*round(1e4*maxMove / dx));
            obj.newShapeRad = floor(1e-4*round(1e4*newShapeRad / dx));
            obj.eraseSize = floor(1e-4*round(1e4*eraseSize / dx));
            obj.velPadding = floor(1e-4*round(1e4*velPadding / dx));
            obj.newShapePad  = floor(1e-4*round(1e4*newShapePad / dx));
            obj.currShapePad = floor(1e-4*round(1e4*currShapePad / dx));
            obj.radiusCurv = floor(1e-4*round(1e4*radiusCurv / dx));
            obj.newShapeType = 'circle';
            
            for i = 1:obj.numDepths
                obj.geoXY{i} = FreeForm(x0, y0, z0, xLengthReal, yLengthReal, dx, (obj.zDepths(end)-obj.zDepths(i))*obj.dx, eps, epsOut, xBC, yBC, newShapeCreation, newShapeRad, newShapePad, maxMove, eraseSize, velPadding, currShapePad, radiusCurv, numFreq, numUser);
            end
            
        end
        
        %% setEpsGrid(epsGrid, x_grid, y_grid, z_vec)
        % takes a binary matrix, 1s = obj.eps, 0s = obj.epsOut
        function obj = setEpsGrid(obj, epsGrid, x_grid, y_grid, z_vec)
            zLen = length(z_vec);
            if(zLen ~= obj.numDepths)
                error('length(z_vec) must match number of etch depths');
            else
                for i = 1:zLen
                    if(obj.numY==1) % assume epsGrid(x,zDepth)
                        obj.epsGrid(1,:,i) = interp1(x_grid, epsGrid(:,i), obj.xGrid, 'pchip');
                    else % assume epsGrid(y,x,zDepth)
                        obj.epsGrid(:,:,i) = interp2(x_grid, y_grid, epsGrid(:,:,i), obj.xGrid, obj.yGrid, 'cubic');
                    end                        
                end
                obj.epsGrid = round(obj.epsGrid)==1;
                obj.epsGrid = obj.enforceEtch(obj.epsGrid);
                for i = 1:zLen
                    obj.geoXY{i} = obj.geoXY{i}.setEpsGrid(obj.epsGrid(:,:,i),x_grid,y_grid);
                end
            end
        end
        
        %% setMaskGrid(maskGrid, x_grid, y_grid)
        % takes a binary matrix, 1s = designable, 0s = non-designable
        function obj = setMaskGrid(obj, maskGrid, x_grid, y_grid)
            if(obj.numY==1)
                maskGrid = interp1(x_grid,maskGrid,obj.xGrid);
            else
                maskGrid = interp2(x_grid, y_grid, maskGrid, obj.xGrid, obj.yGrid);
            end
            obj.maskGrid = round(maskGrid)==1;
        end
        
        %% resetData()
        function[obj] = resetData(obj)
            obj.E = cell(obj.numFreq, obj.numUser);
            obj.EA = cell(obj.numFreq, obj.numUser);
            for i = 1:obj.numDepths
               obj.geoXY{i} = obj.geoXY{i}.resetData();
            end
        end
        
        %% updateData()
        function obj = updateData(obj, data, x_vec, y_vec, z_vec, freqInd, userIndex, EAflag)
            fieldsMask = calcFieldsMask(obj);
            for i = 1:size(data,1) % Multiple frequencies may be updated at once
                Ei = obj.interpolateData(data(i,:), x_vec, y_vec, z_vec);
                if(EAflag)
                    obj.EA(freqInd(i),userIndex) = {Ei};
                else
                    obj.E(freqInd(i),userIndex) = {Ei};
                end
                for j = 1:obj.numDepths
                    Ex = Ei{1}.*fieldsMask{j};
                    Ey = Ei{2}.*fieldsMask{j};
                    Ez = Ei{3}.*fieldsMask{j};
                    Ei_array{j}(i,:)={Ex,Ey,Ez};
                end
            end
            
            for i = 1:obj.numDepths
                obj.geoXY{i} = obj.geoXY{i}.updateData(Ei_array{i}, obj.xGrid, obj.yGrid, freqInd, userIndex, EAflag);
            end
            
        end
        
        %% calcFieldsMask()
        % Allows for treating each etch mask as a separate optimization
        function fieldsMask = calcFieldsMask(obj)            
            for i = 1:obj.numDepths
                fieldMask = ones(obj.numY,obj.numX,obj.numZ);

                if(i < obj.numDepths)
                    dEps = obj.geoXY{i+1}.epsGrid-obj.geoXY{i}.epsGrid;
                    epsBnd = obj.geoXY{i}.getBoundary;
                    dEpsBnd = obj.geoXY{1}.getBoundaryOut(dEps) + obj.geoXY{1}.getBoundaryIn(dEps);
                    mask = 1* (dEpsBnd & epsBnd);
                    
                    % mask out non-overlap regions of above layers
                    mask1 = ~dEps;
                    
                    % mask out border regions of above layers
                    if(obj.numY==1)
                        radCurv = obj.eraseSize + obj.velPadding + 1;
                        filterSize = 2*radCurv + 1; % must be ODD
                        filterMid = (filterSize+1)/2;
                        filterX = 1:filterSize;
                        filter = 1* ( abs(filterX-filterMid) < radCurv );
                        mask2 = 1* ( ~conv(mask,filter,'same') );
                    else
                        radCurv = obj.eraseSize + obj.velPadding +obj.dx;
                        filterSize = 2*radCurv + 1; % must be ODD
                        filterMid = (filterSize+1)/2;
                        [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
                        filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < radCurv.^2 );
                        mask2 = 1* ( ~conv2(mask,filter,'same') );
                    end
                    mask = mask1 .* mask2;
                    fieldMask(:,:,obj.zDepths(i+1):obj.zDepths(end)-1) = repmat(mask,[1 1 obj.zDepths(end)-obj.zDepths(i+1)]);
                end
                
                % mask out everything of below layers
                if(i > 1)
                    fieldMask(:,:,obj.zDepths(1):obj.zDepths(i)-1) = 0;
                end
                
                fieldsMask{i} = fieldMask;
            end
        end
        
        %% interpolateData(field,x_grid,y_grid)
        function [E] = interpolateData(obj, field, vec1, vec2, vec3)
            newVec1 = obj.x0 + obj.dx *(0:obj.numX-1);
            newVec2 = obj.y0 + obj.dx *(0:obj.numY-1);
            newVec3 = obj.z0 + obj.dx *(0:obj.numZ-1);
            
            if(obj.numY==1)  % pivot y into z
                [grid1,grid2] = meshgrid(vec1,vec2);
                Ex = interp2(grid1,grid2, field{1}(:,:,1), newVec1,newVec3','cubic');
                Ey = interp2(grid1,grid2, field{2}(:,:,1), newVec1,newVec3','cubic');
                Ez = interp2(grid1,grid2, field{3}(:,:,1), newVec1,newVec3','cubic');
                
                Ex = permute(Ex,[3 2 1]);
                Ey = permute(Ey,[3 2 1]);
                Ez = permute(Ez,[3 2 1]);
            else
                Ex = interp3(vec1,vec2',vec3, field{1}, newVec1,newVec2',newVec3, 'cubic');
                Ey = interp3(vec1,vec2',vec3, field{2}, newVec1,newVec2',newVec3, 'cubic');
                Ez = interp3(vec1,vec2',vec3, field{3}, newVec1,newVec2',newVec3, 'cubic');
            end
            
            Ex(isnan(Ex)) = 0;
            Ey(isnan(Ey)) = 0;
            Ez(isnan(Ez)) = 0;
            
            E = {Ex,Ey,Ez};
        end
        
        %% updateShape()
        function [obj, dF] = updateShape(obj)
            dF = 0;
            maskGrid2 = ones(obj.numY,obj.numX);
            for i = 1:obj.numDepths
                obj.geoXY{i} = obj.geoXY{i}.setMaskGrid(obj.maskGrid.*maskGrid2,obj.xGrid,obj.yGrid);
                [obj.geoXY{i}, dFi] = obj.geoXY{i}.updateShape;
                epsGridNew = obj.geoXY{i}.epsGrid;
                maskGrid2 = maskGrid2 .* obj.epsGrid(:,:,i);
                obj.epsGrid(:,:,i) = epsGridNew;
                if(i<obj.numDepths) % update geometry of above layers
                    obj.epsGrid = obj.enforceEtch(obj.epsGrid);
                    obj.geoXY{i+1} = obj.geoXY{i+1}.setEpsGrid(obj.epsGrid(:,:,i+1),obj.xGrid,obj.yGrid);
                end
                obj.VBoundaryAll{i} = obj.geoXY{i}.VBoundary;
                obj.velocityAll{i} = obj.geoXY{i}.velocity;
                dF = dF + dFi;
            end
        end
        
         %% enforceEtch(epsGrid)
        function epsGrid = enforceEtch(obj,epsGrid)
            zLen = size(epsGrid,3);
            if(zLen > 1)
                for i = 1:zLen-1
                    epsGrid(:,:,i+1) = epsGrid(:,:,i+1) | epsGrid(:,:,i);
                end
            end
        end
        
        %% returnData()
        function [n,x,y,z] = returnData(obj)
            if(obj.numY == 1) %2D grating, dimensions = x, z
                x = squeeze(obj.xGrid(1,:));
                y = obj.zVec;
                z = 0;
                
                nDepth = sqrt((obj.epsGrid==1)*obj.eps + (obj.epsGrid==0)*obj.epsOut);
                zLen = length(obj.zDepths);
                for i=1:zLen-1
                    n(:,:,obj.zDepths(i):obj.zDepths(i+1)-1) = repmat(nDepth(:,:,i),[1 1 (obj.zDepths(i+1)-obj.zDepths(i))]);
                end
                n=squeeze(n);
                
            end
            
        end
        
        
        %% OLD CRAP
        
        %% updateShape2()
        function [obj, dF] = updateShape2(obj)
            if(obj.newShapeCreation)
                [obj, dF] = obj.createNewIsland;
            else
                dF=0;
            end
            
            if(~dF)
                [padOut, thetaOut] = obj.geoXY.getPadOut;
                [padIn, thetaIn] = obj.geoXY.getPadIn;
                pad = padOut + padIn;
                theta = obj.geoXY.getThetaOut + obj.geoXY.getThetaIn;
                epsPad = padOut*obj.epsOut + padIn*obj.eps;
                figure(11); imagesc(theta); caxis([-pi pi]); colormap(bluewhitered);
                figure(12); imagesc(epsPad); colormap(bluewhitered);
            else
                velCnt = 0;
                obj.VBoundary = 0;
                obj.VBoundaryAll = {};
                for freqInd = 1:obj.numFreq
                    for userInd = 1:obj.numUser
                        E = obj.E{freqInd, userInd};
                        EA = obj.EA{freqInd, userInd};
                        if(~isempty(EA))
                            [Ep, EAp, Dn, DAn] = obj.findFieldNP(E, EA, pad, theta, epsPad);
                            VBoundaryNew = obj.calcV(Ep, EAp, Dn, DAn, pad, padIn);
                            obj.VBoundaryAll(freqInd,userInd) = {VBoundaryNew};
                            obj.VBoundary = obj.VBoundary + VBoundaryNew;
                            velCnt = velCnt + 1;
                        end
                    end
                end
                obj.VBoundary = obj.VBoundary/velCnt;
                
                % Zero-Out velocity around edges
                obj.VBoundary(1:obj.currShapePad,:) = 0;
                obj.VBoundary(end-obj.currShapePad:end,:) = 0;
                obj.VBoundary(:,1:obj.currShapePad) = 0;
                obj.VBoundary(:,end-obj.currShapePad:end) = 0;
                
                obj.VBoundary = obj.fixMinDimension(obj.epsGrid, obj.VBoundary);
                
                thresh = .1;
                [epsGrid0,dF] = obj.shapeLoop(obj.VBoundary, thresh, obj.Emask);
                epsGridNew = obj.smoothRadiusCurv(epsGrid0, obj.Emask);
                
                obj.epsGrid = epsGridNew;
            end
        end
        
        %% createNewIsland()
        function [obj, dFPred] = createNewIsland(obj)
            dFPred = 0;
        end
        
        %%
        function [padOut, thetaOut] = getPadOut()
            
        end
        
        %% getBoundaryOut(grid)
        function gridOut = getBoundaryOut(obj, grid)
            Nx = obj.numX;
            Ny = obj.numY;
            
            epsilon_pad=zeros(Ny+2, Nx+2);
            
            epsilon_pad(2:Ny+1,2:Nx+1)=grid;
            
            epsilon_pad2=epsilon_pad(2:Ny+1,2:Nx+1)|epsilon_pad(1:Ny,1:Nx)...
                |epsilon_pad(1:Ny,3:Nx+2)|epsilon_pad(3:Ny+2,3:Nx+2)...
                |epsilon_pad(3:Ny+2,1:Nx)|epsilon_pad(1:Ny,2:Nx+1)...
                |epsilon_pad(2:Ny+1,1:Nx)|epsilon_pad(3:Ny+2,2:Nx+1)...
                |epsilon_pad(2:Ny+1,3:Nx+2);
            
            gridOut = epsilon_pad2 - epsilon_pad(2:Ny+1,2:Nx+1);
        end
        
    end
    
end