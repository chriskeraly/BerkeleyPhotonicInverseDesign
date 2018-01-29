classdef Spline
    
    %  Notes
    %  Assumes z to be extruded direction
    
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
        
        %Boundary Conditions 0=nothing, +-1=sym, 2=periodic
        xBC;
        yBC;
        
        %Freeform shape properties
        eps_; % single eps_ value or 'material name'
        epsOut; % single eps_ value or 'material name'
        knots; % knots(1,:) = x coordinates, knots(2,:) = y coordinates
        
        %Data = cell(numMon, numFreq, numUD) [numY, numX]
        E;
        EA;
        epsBgnd;
        epsVec; % vector of actual eps_ values
        epsOutVec; % vector of actual eps_ values
        
        % Derivative matrix
        % size = [numMon, numFreq, numUD, numY, numX]
        dFdx;
        dFdxOld;
        deltaXOld;
        velocity;
        newShape;
        
        % These should be set by the Optimizer object
        % Note that numMon is the number of monitors SIMULATIONS, and that
        % numFreq is the total number of frequencies
        numMon;
        numFreq;
        numUD;
        
        % Current Shapes
        quasiNewton; % 1 = allow for quasi-newton optimization
        maxMove; % limits growth of shapes
        maxArea; % limits growth of shapes
        eraseSize; % size of ignored field in velocity calcualtion
        velPadding; % size of region in velocity calculation
        
        % Constraints
        minDimension;
        minPadding;
        radiusCurv;
        radiusCurvHard;
        newShapePad;
        
        % New shapes
        newShapeCreation;
        newShapeType; % string
        newShapeRad; % Size of new shapes
        fastShape; % Use Sparse Perturbation for massive first iteration guess
        
        % Test flag (displaying output if 1)
        testFlag;
        
    end
    
    methods
        function obj = Spline(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eps_, epsOut, xBC, yBC, newShapeCreation, newShapeRad, newShapePad, maxMove, maxArea, eraseSize, velPadding, minPadding, radiusCurv, radiusCurvHard, minDimension, fastShape)
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.dx = dx;
            
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx));
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx));
            obj.numZ = floor(1e-4*round(1e4*thickness / dx));
            [obj.xGrid, obj.yGrid] = meshgrid(x0+dx*(0:obj.numX-1)+dx/2, y0+dx*(0:obj.numY-1)+dx/2);
            if(obj.numZ>0)
                obj.zVec = z0 + dx*(-obj.numZ/2:obj.numZ/2);
            else
                obj.zVec = z0 + dx*[-1 1];
            end
            
            obj.eps_ = eps_;
            obj.epsOut = epsOut;
            obj.epsGrid = zeros(obj.numY,obj.numX);
            obj.maskGrid=ones(obj.numY,obj.numX);
            obj.xBC = xBC;
            obj.yBC = yBC;
            
            obj.eraseSize = floor(1e-4*round(1e4*eraseSize / dx));
            obj.velPadding = floor(1e-4*round(1e4*velPadding / dx));
            obj.maxMove = floor(1e-4*round(1e4*maxMove / dx));
            obj.maxArea = floor(1e-4*round(1e4*maxArea / dx^2));
            
            obj.radiusCurv = floor(1e-4*round(1e4*radiusCurv / dx));
            obj.radiusCurvHard = floor(1e-4*round(1e4*radiusCurvHard / dx));
            obj.minDimension = floor(1e-4*round(1e4*minDimension / dx));
            obj.minPadding = floor(1e-4*round(1e4*minPadding / dx));
            obj.newShapePad  = floor(1e-4*round(1e4*newShapePad / dx));
            obj.newShapeCreation = newShapeCreation;
            obj.newShapeType = 'circle';
            obj.newShapeRad = floor(1e-4*round(1e4*newShapeRad / dx));
            
            obj.quasiNewton = 0;
            
            obj.fastShape = fastShape;
        end
        
        %% GET DATA CONSISTENT WITH OPTIMIZER
        % setOptData(testFlag,numMon,numFreq,numUD)
        function[obj] = setOptData(obj,testFlag,numMon,numFreq,numUD)
            obj.testFlag = testFlag;
            obj.numMon = numMon;
            obj.numFreq = numFreq;
            obj.numUD = numUD;
            obj = obj.initData(); % Can initialize data now
        end
        
        %% SET GEOMETRY
        % takes a binary matrix, 1s = obj.eps_, 0s = obj.epsOut
        % setEpsGrid(epsGrid, x_grid, y_grid)
        function obj = setEpsGrid(obj, epsGrid, x_grid, y_grid)
            if(obj.numY==1)
                epsGrid = interp1(x_grid,epsGrid,obj.xGrid);
            else
                epsGrid = interp2(x_grid, y_grid, epsGrid, obj.xGrid, obj.yGrid);
            end
            epsGrid = round(epsGrid)==1;
            epsGrid = obj.smoothRadiusCurv(epsGrid, obj.maskGrid);
            if(obj.xBC==1)
                epsGrid = 1*(epsGrid | fliplr(epsGrid));
            end
            if(obj.yBC==1)
                epsGrid = 1*(epsGrid | flipud(epsGrid));
            end
            obj.epsGrid = epsGrid;
        end
        
        %% SET NON-DESIGNABLE REGION
        % takes a binary matrix, 1s = designable, 0s = non-designable
        % setMaskGrid(maskGrid, x_grid, y_grid)
        function obj = setMaskGrid(obj, maskGrid, x_grid, y_grid)
            if(obj.numY==1)
                maskGrid = interp1(x_grid,maskGrid,obj.xGrid);
            else
                maskGrid = interp2(x_grid, y_grid, maskGrid, obj.xGrid, obj.yGrid);
            end
            maskGrid = round(maskGrid)==1;
            if(obj.xBC==1)
                maskGrid = 1*(maskGrid | fliplr(maskGrid));
            end
            if(obj.yBC==1)
                maskGrid = 1*(maskGrid | flipud(maskGrid));
            end
            obj.maskGrid = maskGrid;
        end
        
        %% UPDATE FIELD DATA
        % updateFieldData(data, x_grid, y_grid, monIndex, freqInd, userIndex, EAflag)
        function obj = updateFieldData(obj, data, x_vec, y_vec, z_vec, monIndex, freqInd, userIndex, EAflag)
            for i = 1:size(data,1)
                if(obj.numY==1)
                    % incoming y-data is our z-data (ie. in 'thick' direction)
                    % Swapping Ey and Ez and Swapping y-axis and z-axis for each Ex/Ey/Ez
                    [x_grid, y_grid] = meshgrid(x_vec,z_vec);
                    data_i{1} = permute(data{i,1}, [3 2 1]);
                    data_i{2} = permute(data{i,3}, [3 2 1]);
                    data_i{3} = permute(data{i,2}, [3 2 1]);
                    Ei = obj.interpolateData(data_i, x_grid, z_vec); % Also swap y_vec/z_vec
                else
                    [x_grid, y_grid] = meshgrid(x_vec,y_vec);
                    Ei = obj.interpolateData(data(i,:), x_grid, y_grid);
                end
                if(EAflag)
                    obj.EA(monIndex,freqInd(i),userIndex) = {Ei};
                else
                    obj.E(1,freqInd(i),userIndex) = {Ei};
                end
            end
        end
        
        %% UPDATE EPSILON DATA
        % no user index, assumed constant over user sim's
        % updateEpsData(obj, data, x_grid, y_grid, freqInd)
        function obj = updateEpsData(obj, data, x_vec, y_vec, z_vec, freqInd)
            for i = 1:size(data,1)
                if(obj.numY==1)
                    % Swapping epsy and epsz and Swapping y-axis and z-axis for each epsx/epsy/epsz
                    [x_grid, y_grid] = meshgrid(x_vec,z_vec);
                    data_i{1} = permute(data{i,1}, [3 2 1]);
                    data_i{2} = permute(data{i,3}, [3 2 1]);
                    data_i{3} = permute(data{i,2}, [3 2 1]);
                    obj.epsBgnd(freqInd(i)) = {obj.interpolateData(data_i, x_grid, y_grid)}; % Also swap y_vec/z_vec
                else
                    [x_grid, y_grid] = meshgrid(x_vec,y_vec);
                    obj.epsBgnd(freqInd(i)) = {obj.interpolateData(data(i,:), x_grid, y_grid)};
                end
            end
        end
        
        %% INIT DATA
        % initData()
        function[obj] = initData(obj)
            obj.E = cell(1,obj.numFreq, obj.numUD);
            obj.EA = cell(obj.numMon,obj.numFreq, obj.numUD);
            obj.epsBgnd = cell(obj.numFreq, 1);
        end
        
        %% INTERPOLATE FIELDS ONTO GEOMETRY MESH
        % DO NOT AVERAGE FIELDS IN Z, give entire array to the shape
        % Note that this is also used for eps_ = {epsx,epsy,epsz}
        % interpolateData(field, x_grid, y_grid)
        function [E] = interpolateData(obj, field, x_grid, y_grid)
            xInt = obj.xGrid;
            yInt = obj.yGrid;
            numZ = size(field{1},3);
            Ex_array = zeros(obj.numY,obj.numX,numZ);
            Ey_array = zeros(obj.numY,obj.numX,numZ);
            Ez_array = zeros(obj.numY,obj.numX,numZ);
            
            for z = 1 : numZ
                if(obj.numY==1)
                    Ex = interp1(x_grid, field{1}(:,:,z), xInt, 'pchip');
                    Ey = interp1(x_grid, field{2}(:,:,z), xInt, 'pchip');
                    Ez = interp1(x_grid, field{3}(:,:,z), xInt, 'pchip');
                else
                    Ex = interp2(x_grid, y_grid, field{1}(:,:,z), xInt, yInt, 'spline');
                    Ey = interp2(x_grid, y_grid, field{2}(:,:,z), xInt, yInt, 'spline');
                    Ez = interp2(x_grid, y_grid, field{3}(:,:,z), xInt, yInt, 'spline');
                end
                Ex_array(:,:,z) = Ex;
                Ey_array(:,:,z) = Ey;
                Ez_array(:,:,z) = Ez;
            end
            
            Ex_array(isnan(Ex_array)) = 0;
            Ey_array(isnan(Ey_array)) = 0;
            Ez_array(isnan(Ez_array)) = 0;
            
            E = {Ex_array,Ey_array,Ez_array};
        end
        
        %% UPDATE SHAPES
        % updateShapes()
        function [obj, deltaX, dF, bndUpdate] = updateShapes(obj, deltaX, F0, dFdxOpt)
            
            dF = 0;
            
            if(obj.fastShape && (obj.dFGrid>0) )
                dANew = sum(sum( abs(obj.newShapeGrid) ));
                dFNew = obj.dFGrid;
            else
                x = obj.newShape(1);
                y = obj.newShape(2);
                r = obj.newShape(3);
                if(obj.numY>1)
                    dANew = pi*r^2;
                else
                    dANew = 2*r;
                end
                dFNew = obj.newShape(4); % dANew is already multipled in calcShapeDerivs()
            end
            
            % BAD CODE!!! Define an agression parameter
            % deltaX = sign(deltaX).*abs(deltaX).^(1);
            
            dFdxOpt = permute(dFdxOpt,[2,3,1]);
            figure(100); imagesc(deltaX);
            figure(101); imagesc(dFdxOpt(:,:,1)); caxis([-5e-15 5e-15]);
            figure(102); imagesc(obj.dFdxOld(:,:,1)); caxis([-5e-15 5e-15]);
            if(obj.quasiNewton==1 && ~isempty(obj.dFdxOld))
                B = (dFdxOpt - obj.dFdxOld) ./ (-obj.deltaXOld);
                B(isnan(B)) = 0;
                B(isinf(B)) = 1;
                deltaX = dFdxOpt ./ B;
                deltaX(isnan(deltaX)) = 0;
                deltaX(isinf(deltaX)) = 0;
                figure(103); imagesc(B);
                figure(104); imagesc(deltaX);
                
                deltaX(deltaX>obj.maxMove) = obj.maxMove;
                deltaX(deltaX<-obj.maxMove) = obj.maxMove;
            end
            [epsGridNew,dFBnd, deltaXNew, dFdxNew] = obj.shapeLoop(deltaX, dFdxOpt, F0);
            if(obj.numY>1)
                dABnd = sum(sum((epsGridNew~=obj.epsGrid))) * obj.dx^2;
            else
                dABnd = sum(sum((epsGridNew~=obj.epsGrid))) * obj.dx;
            end
            if(dABnd==0), dABnd=Inf; end;
            if(isempty(dFBnd)), dFBnd=0; dABnd = Inf; end;
            Fnew = F0 + dFBnd;
            dFBnd = min(Fnew)-min(F0);
            
            if(obj.testFlag)
                fprintf(' dFNew: %g \n dFBnd: %g \n dFNew/dANew: %g \n dFBnd/dABnd: %g \n',...
                    dFNew, dFBnd, dFNew/dANew, dFBnd/dABnd);
            end
            
            % Moving boundary > adding shape?
            bndChange = ( (dFBnd/dABnd > dFNew/dANew) && (dFBnd > dFNew) );
            
            obj.dFdxOld = [];
            obj.deltaXOld = [];
            if(obj.fastShape && (obj.dFGrid>0) )
                % fastShape
                obj.epsGrid = 1*((obj.epsGrid + obj.newShapeGrid)>0);
                dF = obj.dFGrid;
                obj.fastShape = 0;
                bndUpdate = 0;
                fprintf('  Fast Shape using the Sparse Pertubation Approximation \n');
            elseif( (obj.newShapeCreation==0) || bndChange || (dFNew<.01*min(F0)) )
                % Boundary change
                epsGridNew = obj.smoothRadiusCurv(epsGridNew, obj.maskGrid);
                obj.epsGrid = epsGridNew;
                dF = dF + dFBnd;
                bndUpdate = 1;
                obj.dFdxOld = dFdxNew;
                obj.deltaXOld = deltaXNew;
                fprintf('  Moving Boundary \n');
            else
                % Add new shapes
                [xgrid, ygrid] = meshgrid(1:obj.numX, 1:obj.numY);
                if(obj.epsGrid(y,x)==0)
                    % Add circle of eps_
                    newEpsGrid = ((xgrid-x).^2 + (ygrid-y).^2) < r^2;
                    obj.epsGrid = obj.epsGrid | newEpsGrid;
                else
                    % Add circle of epsOut
                    newEpsGrid = ((xgrid-x).^2 + (ygrid-y).^2) > r^2;
                    obj.epsGrid = obj.epsGrid & newEpsGrid;
                end
                if(obj.xBC==1)
                    obj.epsGrid = 1*(obj.epsGrid | fliplr(obj.epsGrid));
                end
                if(obj.yBC==1)
                    obj.epsGrid = 1*(obj.epsGrid | flipud(obj.epsGrid));
                end
                dF = dFNew;
                bndUpdate = 0;
                fprintf('  New shape: %s, at (%g,%g) \n',obj.newShapeType,x,y);
            end
            
        end
        
        %% CALCULATE DERIVATIVES
        % Note: if numMonSims > 1, minDim = 1
        % dFdx has size [numMon, numFreq, numUD, numY, numX]
        function[obj,lb,ub,c] = calcDerivs(obj, F0, skipInd)
            obj.dFdx = zeros(obj.numMon, obj.numFreq, obj.numUD, obj.numY, obj.numX);
            
            % Find derivatives at all boundaries
            obj.dFdx = obj.calcBoundaryDerivs();
            
            % Find derivatives for new shape creation
            [obj.newShape, obj.velocity, obj.newShapeGrid, obj.dFGrid] = obj.calcShapeDerivs(F0,skipInd);
            
            ub = obj.maxMove*ones(obj.numY*obj.numX,1);
            lb = -ub;
            c = @(x)( obj.ConstraintFx(x) );
        end
        
        %% CONSTRAINT FUNCTION
        function[c,ceq] = ConstraintFx(obj,x)
            ceq = 0;
            deltaX = abs(x);
            dA = sum(sum(deltaX));
            c = dA - obj.maxArea; % enforcing that this must be less than or equal to 0
        end
        
        %% RETURN CONSTRAINTS
        function[dA] = deltaArea(obj,deltaX)
            deltaX = abs(deltaX);
            dA = sum(sum(deltaX));
        end
        
        %% CALCULATE NEW SHAPE DERIVATIVES
        % newShape = [x, y, radius, dFdA]
        function [newShape, velocity, newShapeGrid, dFGrid] = calcShapeDerivs(obj, F0, skipInd)
            erasePad = obj.getErasePad();
            
            [padIndices] = find(~erasePad);
            [pady, padx] = find(~erasePad(:,:,1));
            
            [xgrid, ygrid] = meshgrid(1:obj.numX, 1:obj.numY);
            
            % Construct mask to for adding eps_
            Nx = obj.numX;
            Ny = obj.numY;
            radCurv = obj.newShapePad;
            filterSize = 2*(radCurv-1) + 1; % must be ODD
            filterMid = (filterSize+1)/2;
            [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
            filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < (radCurv-1).^2 );
            epsGridPad = zeros(Ny+2*(radCurv-1),Nx+2*(radCurv-1));
            epsGridPad(1+(radCurv-1):end-(radCurv-1),1+(radCurv-1):end-(radCurv-1)) = obj.epsGrid;
            epsGridMask = 1 * ( conv2(1*~epsGridPad,filter,'valid') >= sum(sum(filter)) );
            
            % Construct mask to for adding epsOut
            radCurv = obj.newShapeRad + obj.minDimension;
            filterSize = 2*(radCurv-1) + 1; % must be ODD
            filterMid = (filterSize+1)/2;
            [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
            filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < (radCurv-1).^2 );
            epsGridPad = zeros(Ny+2*(radCurv-1),Nx+2*(radCurv-1));
            epsGridPad(1+(radCurv-1):end-(radCurv-1),1+(radCurv-1):end-(radCurv-1)) = ~obj.epsGrid;
            epsGridOutMask = 1 * ( conv2(1*~epsGridPad,filter,'valid') >= sum(sum(filter)) );
            
            velocity = zeros(obj.numY, obj.numX);
            for m = 1:obj.numMon
                for f = 1:obj.numFreq
                    for u = 1:obj.numUD
                        E = obj.E{1,f,u};
                        EA = obj.EA{m,f,u};
                        % Clausius Mossoti factor for a circle
                        cmx = 2*(obj.epsVec(f) - obj.epsOutVec(f))./(obj.epsVec(f) + obj.epsOutVec(f));
                        cmy = cmx;
                        cmz = (obj.epsVec(f) - obj.epsOutVec(f));
                        if(~isempty(EA))
                            dFdxShape_mfu = 2*real( cmx.* EA{1} .* E{1} ...
                                + cmy.* EA{2} .* E{2} ...
                                + cmz.* EA{3} .* E{3} );
                            dFdxShape_mfu = mean(dFdxShape_mfu,3);
                            
                            numZ = obj.numZ;
                            if(numZ==0)
                                numZ = 1;
                            end
                            dFdxShape_mfu = dFdxShape_mfu * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0)) * numZ;
                            
                            if(obj.numY==1)
                                dFdxShape_mfu = interp1(padx, dFdxShape_mfu(padIndices),xgrid);
                            else
                                velTri = TriScatteredInterp(padx, pady, dFdxShape_mfu(padIndices));
                                dFdxShape_mfu = velTri(xgrid, ygrid);
                            end
                            dFdxShape_mfu(isnan(dFdxShape_mfu)) = 0;
                            
                            if(obj.xBC==1)
                                dFdxShape_mfu = .5*(dFdxShape_mfu + fliplr(dFdxShape_mfu));
                            end
                            if(obj.yBC==1)
                                dFdxShape_mfu = .5*(dFdxShape_mfu + flipud(dFdxShape_mfu));
                            end
                            
                            if( skipInd )
                                if(skipInd==1)
                                    c = m;
                                elseif(skipInd==2)
                                    c = f;
                                elseif(skipInd==3)
                                    c = u;
                                end
                            else
                                c = 1;
                            end
                            velocity(:,:,c) = dFdxShape_mfu .* obj.maskGrid;
                            
                            if(obj.numY>1)
                                figure(30+f); imagesc(dFdxShape_mfu); colormap(bluewhitered); axis equal; title('Velocity: Sparse Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                            else
                                figure(30+f); plot(dFdxShape_mfu); title('Velocity: Sparse Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                            end
                        end
                    end
                end
            end
            
            r = obj.newShapeRad;
            if(obj.numY>1)
                dANew = pi*r^2;
            else
                dANew = 2*r;
            end
            F0 = reshape(F0,1,1,[]);
            minF0 = min(F0);
            F0_arr = repmat(F0,[obj.numY,obj.numX]);
            % Consider circle of eps_
            velocitym = velocity.*repmat(epsGridMask,[1 1 size(velocity,3)]);
            Fnew = dANew*velocitym + F0_arr;
            minFnew = min(Fnew,[],3);
            [maxMinFnew,Y] = max(minFnew,[],1);
            [maxMinFnew,x] = max(maxMinFnew);
            y = Y(x);
            dF = maxMinFnew - minF0;
            % Consider circle of epsOut
            velocitym = velocity.*repmat(epsGridOutMask,[1 1 size(velocity,3)]);
            Fnew = -dANew*velocitym + F0_arr;
            minFnew = min(Fnew,[],3);
            [maxMinFnew,Y2] = max(minFnew,[],1);
            [maxMinFnew,x2] = max(maxMinFnew);
            y2 = Y2(x2);
            dF2 = maxMinFnew - minF0;
            
            if(dF >= dF2)
                newShape = [x, y, r, dF];
            else
                newShape = [x2, y2, r, dF2];
            end
            
            % FastShape, add non-circular regions of eps_ and/or epsOut
            dF = sum(velocity,3);
            dFmax = max(max(dF));
            dFmin = max(max(dF));
            newShapeGrid = 1 * ( (epsGridMask .* (dF>.5*dFmax)) - (epsGridOutMask .* (dF<.5*dFmin)) );
            newShapeGrid = obj.smoothRadiusCurv(newShapeGrid, obj.maskGrid);
            dFGrid = sum(sum( newShapeGrid .* dF ));
        end
        
        %% CALCULATE BOUNDARY DERIVATIVES
        function dFdx = calcBoundaryDerivs(obj)
            padOut = obj.getPadOut;
            padIn = obj.getPadIn;
            pad = padOut + padIn;
            theta = obj.getThetaOut + obj.getThetaIn;
            
            dFdx = zeros(obj.numMon,obj.numFreq, obj.numUD, obj.numY, obj.numX);
            if(any(any(pad))) % only calc dFdx if there exists a boundary
                for m = 1:obj.numMon
                    for f = 1:obj.numFreq
                        for u = 1:obj.numUD
                            E = obj.E{1,f,u};
                            EA = obj.EA{m,f,u};
                            epsPad = padOut*obj.epsOutVec(f) + padIn*obj.epsVec(f);
                            if(~isempty(EA))
                                [Ep, EAp, Dn, DAn, Ez, EAz] = obj.findFieldNP(E, EA, pad, theta, epsPad);
                                
                                dFdx_mfu = obj.calcV(Ep, EAp, Dn, DAn, Ez, EAz, pad, obj.epsVec(f), obj.epsOutVec(f));
                                dFdx(m,f,u,:,:) = dFdx_mfu;
                                
                                vLim = max(max( abs(dFdx_mfu) ));
                                if(obj.numY>1)
                                    figure(40+f); imagesc(dFdx_mfu); axis equal; caxis([-vLim vLim]); colormap(bluewhitered); title('Velocity: Boundary Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                                else
                                    figure(40+f); plot(dFdx_mfu); title('Velocity: Boundary Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
        %% CALCULATE BOUNDARY DERIVATIVES - HELPER
        % Averages dFdx_mfu about the extruded dimension (z-axis)
        function dFdx_mfu = calcV(obj, Ep, EAp, Dn, DAn, Ez, EAz, pad, eps_, epsOut)
            dFdxPad_p = real( (eps_ - epsOut) .* Ep .* EAp ); % par component
            if(real(obj.epsVec(1))<0)
                dFdxPad_pz = real( Ez .* EAz .* (1./epsOut - 1./eps_) ); % perp-z component
            else
                dFdxPad_pz = real( (eps_ - epsOut) .* Ez .* EAz ); % par-z component
            end
            dFdxPad_n = real( Dn .* DAn .* (1./epsOut - 1./eps_) ); % perp component
            
            dFdxPad = dFdxPad_p + dFdxPad_pz + dFdxPad_n;
            
            figure(60); imagesc(abs(Ep(:,:,1))); colormap(bluewhitered);
            figure(61); imagesc(abs(Ez(:,:,1))); colormap(bluewhitered);
            figure(62); imagesc(abs(Dn(:,:,1))); colormap(bluewhitered);
            figure(63); imagesc(angle(Ep(:,:,1))); colormap(bluewhitered);
            figure(64); imagesc(angle(Ez(:,:,1))); colormap(bluewhitered);
            figure(65); imagesc(angle(Dn(:,:,1))); colormap(bluewhitered);
            
            figure(70); imagesc(abs(EAp(:,:,1))); colormap(bluewhitered);
            figure(71); imagesc(abs(EAz(:,:,1))); colormap(bluewhitered);
            figure(72); imagesc(abs(DAn(:,:,1))); colormap(bluewhitered);
            figure(73); imagesc(angle(EAp(:,:,1))); colormap(bluewhitered);
            figure(74); imagesc(angle(EAz(:,:,1))); colormap(bluewhitered);
            figure(75); imagesc(angle(DAn(:,:,1))); colormap(bluewhitered);

            figure(50); imagesc(dFdxPad_p(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); colormap(bluewhitered);
            figure(51); imagesc(dFdxPad_pz(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); colormap(bluewhitered);
            figure(52); imagesc(dFdxPad_n(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); colormap(bluewhitered);
            
            dimZ = size(dFdxPad,3);
            [padIndices] = find(pad);
            [pady, padx] = find(pad(:,:,1));
            [xgrid, ygrid] = meshgrid(1:obj.numX, 1:obj.numY);
            
            % Use padded derivs to interpolate deriv along actual boundary
            dFdx_mfu = 0*dFdxPad;
            for i=1:dimZ
                dFdxPad_z = dFdxPad(:,:,i);
                if(obj.numY==1)
                    dFdx_z = interp1(padx, dFdxPad_z(padIndices),xgrid);
                else
                    velTri = TriScatteredInterp(padx, pady, dFdxPad_z(padIndices));
                    dFdx_z = velTri(xgrid, ygrid);
                end
                dFdx_z(isnan(dFdx_z)) = 0;
                
                numZ = obj.numZ;
                if(numZ==0)
                    numZ = 1;
                end
                dFdx_mfu(:,:,i) = dFdx_z .* obj.getBoundary * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0)) * numZ;
            end
            dFdx_mfu = mean(dFdx_mfu,3) .* obj.maskGrid;
            
            if(obj.xBC==1)
                dFdx_mfu = .5*(dFdx_mfu +fliplr(dFdx_mfu));
            end
            if(obj.yBC==1)
                dFdx_mfu = .5*(dFdx_mfu +flipud(dFdx_mfu));
            end
            
            % Possibly Zero-Out velocity around edges (minEdgePadding?)
            %if(obj.numY > (1+2*obj.minPadding))
            %    dFdx_mfu(1:obj.minPadding,:) = 0;
            %    dFdx_mfu(end-obj.minPadding:end,:) = 0;
            %end
            %dFdx_mfu(:,1:obj.minPadding) = 0;
            %dFdx_mfu(:,end-obj.minPadding:end) = 0;
        end
        
        %% CALCULATE CONTINUOUS FIELD COMPONENTS
        % findFieldNP(E, EA, pad, theta, epsPad)
        function [Ep, EAp, Dn, DAn, Ez, EAz] = findFieldNP(obj, E, EA, pad, theta, epsPad)
            Ex = E{1};
            Ey = E{2};
            Ez = E{3};
            
            EAx = EA{1};
            EAy = EA{2};
            EAz = EA{3};
            
            dimZ = size(Ex,3);
            theta = repmat(theta,[1 1 dimZ]);
            epsPad = repmat(epsPad,[1 1 dimZ]);
            pad = repmat(pad,[1 1 dimZ]);
            
            Dn = (Ex.*cos(theta) + Ey.*sin(theta)) .* epsPad;
            Ep = (Ex.*sin(theta) - Ey.*cos(theta)) .* pad;
            if(real(obj.epsVec(1))<0)
                Ez = Ez .* epsPad;
            else
                Ez = Ez .* pad;
            end
            
            DAn = (EAx.*cos(theta) + EAy.*sin(theta)) .* epsPad;
            EAp = (EAx.*sin(theta) - EAy.*cos(theta)) .* pad;
            if(real(obj.epsVec(1))<0)
                EAz = EAz .* epsPad;
            else
                EAz = EAz .* pad;
            end
            
            
        end
        
        
        %% (NOT USED) UPDATE SHAPE FAST - MAIN HELPER
        % Currently only works for single-Adjoint optimizations
        function [epsGrid,dF] = shapeLoopFast(obj, V, dFdxOpt, F0)
            Nx = obj.numX;
            Ny = obj.numY;
            
            epsilon = 1*obj.epsGrid;
            epsilonBnd = obj.getBoundaryIn(epsilon)+obj.getBoundaryOut(epsilon);
            Emask = obj.maskGrid;
            
            dF=0;
            numMoves = obj.maxMove;
            descale=1;
            
            numMin = length(F0);
            %Vorig = permute(dFdxOpt,[2,3,1]);
            
            Vorig = obj.velocity;
            V = epsilonBnd .* Vorig;
            
            % Compute perturbation's impulse response
            % (constrained to radius of curvature)
            rIn = obj.radiusCurv;
            rOut = obj.radiusCurv;
            x2 = sqrt((rIn+rOut)^2-rOut^2);
            x1 = round( rIn*x2/(rIn+rOut) );
            x2 = round(x2);
            radCurv = rIn;
            filterSize = 2*x2+1;
            [filterX, filterY] = meshgrid(1:filterSize,1:filterSize);
            filterMid = (filterSize+1)/2;
            filRad1 = ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < x1.^2 ) .* real( (radCurv^2 - (filterX-filterMid).^2 - (filterY-filterMid).^2) .^ (.5));
            radCurv= rOut;
            filRad2 = real( radCurv - ( radCurv^2 - (x2 - ((filterX-filterMid).^2 + (filterY-filterMid).^2).^(.5)).^2).^.5);
            filRad2 = filRad2 .* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) > x1.^2 ) .* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < x2.^2 );
            fil = filRad1+filRad2;
            filBuf = round((filterSize-1)/2);
            
            % Convolve velocity with a perturbation's impulse response
            normV = 1 * (V~=0);
            VSmooth = conv2(fil,V);
            normVSmooth = conv2(fil,normV);
            VSmooth = VSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
            normVSmooth = ~normV+normV.*normVSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
            V = VSmooth ./ normVSmooth .* normV;
            
            % Zero-Out velocity to protect minDimension constraint
            V = obj.fixMinDimension(epsilon, V);
            
            % Zero-Out velocity to protect minPadding constraint
            if(obj.minPadding>0)
                V = obj.fixMinPadding(epsilon, V);
            end
            
            % Normalize V to normal angle
            % (when theta=pi/4, boundary will naturally move dx*sqrt(2) rather than dx)
            theta = obj.calcTheta(epsilon);
            thetaNorm = abs(cos(theta)) + abs(sin(theta));
            V = V ./thetaNorm .* Emask; % only used for Vmove calc, recalculated inside main loop
            
            for i = 1:20
                % Establish largest velocity (weighted to move by maxMove)
                Vmax = max(max( V.*obj.getBoundaryOut(epsilon) ) );
                Vmin = min(min( V.*obj.getBoundaryIn(epsilon) ) );
                Vmove = max(Vmax,-Vmin);
                
                %figure(100); imagesc(V); caxis([-Vmove Vmove]); colormap(bluewhitered);
                
                Vup = V .* (1*(V>0));
                normV = 1 * (V~=0);
                VSmooth = conv2(fil,Vup);
                normVSmooth = conv2(fil,normV);
                VSmooth = VSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                normVSmooth = ~normV+normV.*normVSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                Vup = VSmooth ./ normVSmooth .* normV;
                
                Vdown = V .* (1*(V<0));
                normV = 1 * (V~=0);
                VSmooth = conv2(fil,Vdown);
                normVSmooth = conv2(fil,normV);
                VSmooth = VSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                normVSmooth = ~normV+normV.*normVSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                Vdown = VSmooth ./ normVSmooth .* normV;
                
                binIn = ceil(i/2)/(numMoves+1);
                binOut = ceil(i/2)/(numMoves+1);
                epsilonRem=obj.getBoundaryIn(epsilon).*(Vdown<(binIn*-Vmove));
                epsilonAdd=obj.getBoundaryOut(epsilon).*(Vup>(binOut*Vmove));
                if(mod(i+1,2))
                    epsilon = epsilon + epsilonAdd;
                    dF=dF+sum(sum(epsilonAdd.*Vorig)); % mutiply by actual velocity, not the contstrained one
                    flag=1;
                else
                    epsilon = epsilon - epsilonRem;
                    dF=dF-sum(sum(epsilonRem.*Vorig));
                    flag=-1;
                end
                epsilon = obj.smoothGrid2(epsilon,flag);
                epsilon = obj.smoothGrid(epsilon,flag);
                
                epsilonBnd = obj.getBoundaryIn(epsilon)+obj.getBoundaryOut(epsilon);
                
                V = obj.shapeLoopVelocity(V).*epsilonBnd.*Emask/descale;
                Vorig = obj.shapeLoopVelocity(Vorig).*repmat(epsilonBnd,[1 1 numMin]).*repmat(Emask/descale,[1 1 numMin]);
                
                figure(101); imagesc(epsilon); colormap(bluewhitered); axis equal;
                
            end
            epsilon = obj.smoothRadiusCurv(epsilon, Emask);
            figure(102); imagesc(epsilon); colormap(bluewhitered); axis equal;
            
            V = epsilonBnd.*Vorig; % BAD-CODE, no check on size of velocity
            epsGrid = Emask.*epsilon + (~Emask).*obj.epsGrid;
        end
        
        
        %% UPDATE SHAPE - MAIN HELPER
        function [epsGrid,dF,V,Vorig] = shapeLoop(obj, V, dFdxOpt, F0)
            Nx = obj.numX;
            Ny = obj.numY;
            
            epsilon = 1*obj.epsGrid;
            Emask = obj.maskGrid;
            
            dF=0;
            descale = 1;
            numMoves = obj.maxMove;
            
            numMin = length(F0);
            Vorig = dFdxOpt;
            
            if(obj.xBC==1)
                V = .5*(V + fliplr(V));
            end
            if(obj.yBC==1)
                V = .5*(V + flipud(V));
            end
            
            % Compute perturbation's impulse response
            % (constrained to radius of curvature)
            if(obj.numY>1)
                rIn = obj.radiusCurv;
                rOut = obj.radiusCurv;
                x2 = sqrt((rIn+rOut)^2-rOut^2);
                x1 = round( rIn*x2/(rIn+rOut) );
                x2 = round(x2);
                radCurv = rIn;
                filterSize = 2*x2+1;
                [filterX, filterY] = meshgrid(1:filterSize,1:filterSize);
                filterMid = (filterSize+1)/2;
                filRad1 = ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < x1.^2 ) .* real( (radCurv^2 - (filterX-filterMid).^2 - (filterY-filterMid).^2) .^ (.5));
                radCurv= rOut;
                filRad2 = real( radCurv - ( radCurv^2 - (x2 - ((filterX-filterMid).^2 + (filterY-filterMid).^2).^(.5)).^2).^.5);
                filRad2 = filRad2 .* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) > x1.^2 ) .* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < x2.^2 );
                fil = filRad1+filRad2;
                filBuf = round((filterSize-1)/2);
                
                % Convolve velocity with a perturbation's impulse response
                normV = 1 * (V~=0);
                VSmooth = conv2(fil,V);
                normVSmooth = conv2(fil,normV);
                VSmooth = VSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                normVSmooth = ~normV+normV.*normVSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                V = VSmooth ./ normVSmooth .* normV;
            end
            
            % Zero-Out velocity to protect minDimension constraint
            V = obj.fixMinDimension(epsilon, V);
            
            % Zero-Out velocity to protect minPadding constraint
            if(obj.minPadding>0)
                V = obj.fixMinPadding(epsilon, V);
            end
            vLim = max(max( abs(V) ));
            
            % Normalize V to normal angle
            % (when theta=pi/4, boundary will naturally move dx*sqrt(2) rather than dx)
            theta = obj.calcTheta(epsilon);
            thetaNorm = abs(cos(theta)) + abs(sin(theta));
            V = V ./thetaNorm .* Emask; % only used for Vmove calc, recalculated inside main loop
            
            % Establish largest velocity (weighted to move by maxMove)
            Vmax = max(max( V.*obj.getBoundaryOut(epsilon) ) );
            Vmin = min(min( V.*obj.getBoundaryIn(epsilon) ) );
            Vmove = max(Vmax,-Vmin);
            if(obj.quasiNewton==1 && ~isempty(obj.dFdxOld))
                Vmove = obj.maxMove;
            end
            
            for i=1:2*numMoves
                Vup = V .* (1*(V>0));
                if(obj.numY>1)
                    normV = 1 * (V~=0);
                    VSmooth = conv2(fil,Vup);
                    normVSmooth = conv2(fil,normV);
                    VSmooth = VSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                    normVSmooth = ~normV+normV.*normVSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                    Vup = VSmooth ./ normVSmooth .* normV;
                end
                
                Vdown = V .* (1*(V<0));
                if(obj.numY>1)
                    normV = 1 * (V~=0);
                    VSmooth = conv2(fil,Vdown);
                    normVSmooth = conv2(fil,normV);
                    VSmooth = VSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                    normVSmooth = ~normV+normV.*normVSmooth(1+filBuf:Ny+filBuf, 1+filBuf:Nx+filBuf);
                    Vdown = VSmooth ./ normVSmooth .* normV;
                end
                % Normalize V to normal angle
                % (when theta=pi/4, boundary will naturally move dx*sqrt(2) rather than dx)
                %                 theta = obj.calcTheta(epsilon);
                %                 thetaNorm = abs(cos(theta)) + abs(sin(theta));
                %                 Vup = Vup ./thetaNorm .* Emask;
                %                 Vdown = Vdown ./thetaNorm .* Emask;
                
                binIn = ceil(i/2)/(numMoves+1);
                binOut = ceil(i/2)/(numMoves+1);
                epsilonRem=obj.getBoundaryIn(epsilon).*(Vdown<(binIn*-Vmove));
                epsilonAdd=obj.getBoundaryOut(epsilon).*(Vup>(binOut*Vmove));
                if(mod(i+1,2))
                    epsilon = epsilon + epsilonAdd;
                    dF=dF+squeeze(sum(sum(repmat(epsilonAdd,[1 1 numMin]).*Vorig,1),2)); % mutiply by actual velocity, not the contstrained one
                    flag=1;
                else
                    epsilon = epsilon - epsilonRem;
                    dF=dF-squeeze(sum(sum(repmat(epsilonRem,[1 1 numMin]).*Vorig,1),2));
                    flag=-1;
                end
                if(obj.numY>1)
                    epsilon = obj.smoothGrid2(epsilon,flag);
                    epsilon = obj.smoothGrid(epsilon,flag);
                end
                
                epsilonMask = obj.getBoundaryIn(epsilon)+obj.getBoundaryOut(epsilon);
                
                V = obj.shapeLoopVelocity(V).*epsilonMask.*Emask/descale;
                Vorig = obj.shapeLoopVelocity(Vorig).*repmat(epsilonMask,[1 1 numMin]).*repmat(Emask/descale,[1 1 numMin]);
            end
            
            epsGrid = Emask.*epsilon + (~Emask).*obj.epsGrid;
            
            if(obj.xBC==1)
                epsGrid = 1*(epsGrid | fliplr(epsGrid));
            end
            if(obj.yBC==1)
                epsGrid = 1*(epsGrid | flipud(epsGrid));
            end
            
            minF0 = min(F0);
            Fnew = F0 + dF;
            minFnew = min(Fnew);
            dF = minFnew - minF0;
        end
        
        %% UPDATE SHAPE - SUB HELPER
        % approximates velocity 1 cell away from the boundary
        function [Vshape] = shapeLoopVelocity(obj, V)
            Nx = obj.numX;
            Ny = obj.numY;
            
            Vshape=zeros(size(V));
            for i=1:size(V,3)
                Vi = V(:,:,i);
                Vpad = zeros(size(Vi)+2);
                Vpad(2:Ny+1,2:Nx+1) = Vi;
                norm_pad = 1*(Vpad~=0);
                
                Vpad2=Vpad(2:Ny+1,2:Nx+1)+Vpad(1:Ny,1:Nx)...
                    +Vpad(1:Ny,3:Nx+2)+Vpad(3:Ny+2,3:Nx+2)...
                    +Vpad(3:Ny+2,1:Nx)+Vpad(1:Ny,2:Nx+1)...
                    +Vpad(2:Ny+1,1:Nx)+Vpad(3:Ny+2,2:Nx+1)...
                    +Vpad(2:Ny+1,3:Nx+2);
                
                norm_pad2=norm_pad(2:Ny+1,2:Nx+1)+norm_pad(1:Ny,1:Nx)...
                    +norm_pad(1:Ny,3:Nx+2)+norm_pad(3:Ny+2,3:Nx+2)...
                    +norm_pad(3:Ny+2,1:Nx)+norm_pad(1:Ny,2:Nx+1)...
                    +norm_pad(2:Ny+1,1:Nx)+norm_pad(3:Ny+2,2:Nx+1)...
                    +norm_pad(2:Ny+1,3:Nx+2);
                
                norm_pad2(norm_pad2==0) = 1;
                Vshape(:,:,i) = Vpad2 ./ norm_pad2;
            end
        end
        
        %% getThetaOut()
        function [thetaOut] = getThetaOut(obj)
            Nx = obj.numX;
            Ny = obj.numY;
            
            thetaOut = zeros(obj.numY, obj.numX);
            gridNew = obj.epsGrid;
            thetaNew = obj.calcTheta(gridNew);
            
            num = obj.eraseSize + obj.velPadding;
            
            for i=1:num
                thetaPad = zeros(Ny+2, Nx+2);
                boundary = obj.getBoundaryOut(gridNew) + obj.getBoundaryIn(gridNew);
                thetaPad(2:Ny+1,2:Nx+1) = boundary.*exp(1i*thetaNew); % exponential, so angles can be averaged properly
                
                %                 thetaNorm = abs(cos(angle(thetaPad))) + abs(sin(angle(thetaPad)));
                %                 thetaPad = thetaPad ./thetaNorm;
                
                thetaPad2=thetaPad(2:Ny+1,2:Nx+1)+thetaPad(1:Ny,1:Nx)...
                    +thetaPad(1:Ny,3:Nx+2)+thetaPad(3:Ny+2,3:Nx+2)...
                    +thetaPad(3:Ny+2,1:Nx)+thetaPad(1:Ny,2:Nx+1)...
                    +thetaPad(2:Ny+1,1:Nx)+thetaPad(3:Ny+2,2:Nx+1)...
                    +thetaPad(2:Ny+1,3:Nx+2);
                thetaPad2 = 1e-4*round(1e4*thetaPad2); % to remove miniscule numbers which should be zero
                thetaNew = angle(thetaPad2);
                
                gridOut = obj.getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                
                if(i > obj.eraseSize)
                    thetaOut = thetaOut + thetaNew.*(1*gridOut);
                end
            end
        end
        
        %% getThetaIn()
        function [thetaIn] = getThetaIn(obj)
            Nx = obj.numX;
            Ny = obj.numY;
            
            thetaIn = zeros(obj.numY, obj.numX);
            gridNew = obj.epsGrid;
            thetaNew = obj.calcTheta(gridNew);
            
            num = obj.eraseSize + obj.velPadding;
            
            for i=1:num
                thetaPad = zeros(Ny+2, Nx+2);
                boundary = obj.getBoundaryOut(gridNew) + obj.getBoundaryIn(gridNew);
                thetaPad(2:Ny+1,2:Nx+1) = boundary.*exp(1i*thetaNew); % exponetial, so angles can be averaged properly
                
                thetaPad2=thetaPad(2:Ny+1,2:Nx+1)+thetaPad(1:Ny,1:Nx)...
                    +thetaPad(1:Ny,3:Nx+2)+thetaPad(3:Ny+2,3:Nx+2)...
                    +thetaPad(3:Ny+2,1:Nx)+thetaPad(1:Ny,2:Nx+1)...
                    +thetaPad(2:Ny+1,1:Nx)+thetaPad(3:Ny+2,2:Nx+1)...
                    +thetaPad(2:Ny+1,3:Nx+2);
                thetaPad2 = 1e-4*round(1e4*thetaPad2); % to remove miniscule numbers which should be zero
                thetaNew = angle(thetaPad2);
                
                gridIn = obj.getBoundaryIn(gridNew);
                gridNew = gridNew & ~gridIn;
                
                if(i > obj.eraseSize)
                    thetaIn = thetaIn + thetaNew.*(1*gridIn);
                end
            end
        end
        
        %% getBoundaryApproxPad()
        function bndPad = getBoundaryApproxPad(obj)
            bndPad = zeros(obj.numY, obj.numX);
            gridNew = obj.epsGrid;
            for i=1:obj.maxMove+obj.newShapeRad
                gridOut = obj.getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                bndPad = bndPad | gridOut;
            end
            gridNew = obj.epsGrid;
            for i=1:obj.maxMove+obj.newShapeRad
                gridIn = obj.getBoundaryIn(gridNew);
                gridNew = ~gridIn & gridNew;
                bndPad = bndPad | gridIn;
            end
        end
        
        %% getErasePad()
        function erasePad = getErasePad(obj)
            erasePad = zeros(obj.numY, obj.numX);
            gridNew = obj.epsGrid;
            for i=1:obj.eraseSize
                gridOut = obj.getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                erasePad = erasePad | gridOut;
            end
            gridNew = obj.epsGrid;
            for i=1:obj.eraseSize
                gridIn = obj.getBoundaryIn(gridNew);
                gridNew = ~gridIn & gridNew;
                erasePad = erasePad | gridIn;
            end
        end
        
        %% getPadOut()
        function padOut = getPadOut(obj)
            num = obj.eraseSize + obj.velPadding;
            padOut = zeros(obj.numY, obj.numX);
            gridNew = obj.epsGrid;
            for i=1:num
                gridOut = obj.getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                if(i > obj.eraseSize)
                    padOut = padOut | gridOut;
                end
            end
        end
        
        %% getPadIn()
        function padIn = getPadIn(obj)
            num = obj.eraseSize + obj.velPadding;
            padIn = zeros(obj.numY, obj.numX);
            gridNew = obj.epsGrid;
            for i=1:num
                gridIn = obj.getBoundaryIn(gridNew);
                gridNew = ~gridIn & gridNew;
                if(i > obj.eraseSize)
                    padIn = padIn | gridIn;
                end
            end
        end
        
        %% calcTheta()
        function theta = calcTheta(obj, grid)
            Nx = obj.numX;
            Ny = obj.numY;
            
            %flipud because y increases downward
            filterSize = 11; % must be ODD
            filterMid = (filterSize+1)/2;
            
            [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
            filterX = filterX - filterMid;
            filterY = filterY - filterMid;
            filter = exp(1i*angle(filterX+1i*filterY));
            filter(filterMid,filterMid)=0;
            %filter = filter .* (abs(filterX+1i*filterY)<filterMid);
            
            theta = mod(2*pi+angle(1e-4*round(1e4*conv2(filter,1*grid))),2*pi);
            theta = theta(filterMid:filterMid+Ny-1, filterMid:filterMid+Nx-1).*(obj.getBoundaryIn(grid)+obj.getBoundaryOut(grid));
        end
        
        %% getBoundary()
        function epsGridBound = getBoundary(obj)
            epsGridBound = obj.getBoundaryIn(obj.epsGrid) +obj.getBoundaryOut(obj.epsGrid);
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
        
        %% getBoundaryIn(grid)
        function gridIn = getBoundaryIn(obj, grid)
            Nx = obj.numX;
            Ny = obj.numY;
            
            epsilon_pad=zeros(Ny+2, Nx+2);
            
            epsilon_pad(2:Ny+1,2:Nx+1)=(grid==1);
            
            if(Ny==1)
                epsilon_pad2=epsilon_pad(2:Ny+1,2:Nx+1)...
                    &epsilon_pad(2:Ny+1,1:Nx)...
                    &epsilon_pad(2:Ny+1,3:Nx+2);
            else
                epsilon_pad2=epsilon_pad(2:Ny+1,2:Nx+1)&epsilon_pad(1:Ny,1:Nx)...
                    &epsilon_pad(1:Ny,3:Nx+2)&epsilon_pad(3:Ny+2,3:Nx+2)...
                    &epsilon_pad(3:Ny+2,1:Nx)&epsilon_pad(1:Ny,2:Nx+1)...
                    &epsilon_pad(2:Ny+1,1:Nx)&epsilon_pad(3:Ny+2,2:Nx+1)...
                    &epsilon_pad(2:Ny+1,3:Nx+2);
            end
            
            gridIn = epsilon_pad(2:Ny+1,2:Nx+1) - epsilon_pad2;
        end
        
        %% 1*dx SMOOTHING FUNCTION
        % gets rid of tiny 1 mesh cell features
        function gridNew = smoothGrid(obj, grid, flag)
            Nx = obj.numX;
            Ny = obj.numY;
            
            grid_pad=ones(size(grid)+2);
            
            for i=1:10
                grid_pad(2:Ny+1,2:Nx+1)=grid;
                
                grid_pad2=grid_pad(1:Ny,2:Nx+1)...
                    +grid_pad(2:Ny+1,1:Nx)+grid_pad(3:Ny+2,2:Nx+1)...
                    +grid_pad(2:Ny+1,3:Nx+2);
                
                if(flag>0)
                    gridNew = 1*(grid | (grid_pad2>2));
                else
                    gridNew = 1*(grid & (grid_pad2>1));
                end
                
                grid_pad=ones(size(grid)+2);
                grid_pad(2:Ny+1,2:Nx+1)=gridNew;
                
                grid_pad2=grid_pad(1:Ny,2:Nx+1)...
                    +grid_pad(2:Ny+1,1:Nx)+grid_pad(3:Ny+2,2:Nx+1)...
                    +grid_pad(2:Ny+1,3:Nx+2);
                
                if(flag<0)
                    gridNew = 1*(gridNew | (grid_pad2>2));
                else
                    gridNew = 1*(gridNew & (grid_pad2>1));
                end
                
                grid=gridNew;
            end
            
        end
        
        %% 2*dx SMOOTHING FUNCTION
        % gets rid of tiny 2 mesh cell features
        function gridNew=smoothGrid2(obj, grid, flag)
            Nx = obj.numX;
            Ny = obj.numY;
            
            grid_pad=ones(size(grid)+2);
            grid_pad3=ones(size(grid)+2);
            
            for i=1:10
                grid_pad(2:Ny+1,2:Nx+1)=grid;
                
                grid_pad2=grid_pad(1:Ny,2:Nx+1)...
                    +grid_pad(2:Ny+1,1:Nx)+grid_pad(3:Ny+2,2:Nx+1)...
                    +grid_pad(2:Ny+1,3:Nx+2);
                
                grid_pad3(2:Ny+1,2:Nx+1)=grid_pad2;
                grid_pad4=grid_pad3(1:Ny,2:Nx+1)...
                    +grid_pad3(2:Ny+1,1:Nx)+grid_pad3(3:Ny+2,2:Nx+1)...
                    +grid_pad3(2:Ny+1,3:Nx+2);
                
                if(flag>0)
                    gridNew = 1*(grid | (grid_pad4>=7 & grid_pad4<=9) );
                else
                    gridNew = 1*(grid & (grid_pad4<7 | grid_pad4>9) );
                end
                
                grid=gridNew;
                
                grid_pad(2:Ny+1,2:Nx+1)=grid;
                
                grid_pad2=grid_pad(1:Ny,2:Nx+1)...
                    +grid_pad(2:Ny+1,1:Nx)+grid_pad(3:Ny+2,2:Nx+1)...
                    +grid_pad(2:Ny+1,3:Nx+2);
                
                grid_pad3(2:Ny+1,2:Nx+1)=grid_pad2;
                grid_pad4=grid_pad3(1:Ny,2:Nx+1)...
                    +grid_pad3(2:Ny+1,1:Nx)+grid_pad3(3:Ny+2,2:Nx+1)...
                    +grid_pad3(2:Ny+1,3:Nx+2);
                
                if(flag<0)
                    gridNew = 1*(grid | (grid_pad4>=7 & grid_pad4<=9) );
                else
                    gridNew = 1*(grid & (grid_pad4<7 | grid_pad4>9) );
                end
                
                grid=gridNew;
            end
        end
        
        
        %% HARD RADIUS OF CURVATURE CONSTRAINT
        % Adds/Subtracts material (blind to dFdx) to enforce radiusCurvHard
        function epsGridNew = smoothRadiusCurv(obj, epsGrid0, maskGrid)
            Nx = obj.numX;
            Ny = obj.numY;
            radCurv = obj.radiusCurvHard;
            
            if(obj.numY > 1) % Only applies to >1D structures
                filterSize = 2*radCurv + 1; % must be ODD
                filterMid = (filterSize+1)/2;
                [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
                filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < radCurv.^2 );
                
                epsGrid1 = (~maskGrid).*epsGrid0;
                
                epsGridPad = zeros(Ny+4*radCurv,Nx+4*radCurv);
                epsGridPad(1+2*radCurv:end-2*radCurv,1+2*radCurv:end-2*radCurv) = epsGrid1;
                epsGridOut1 = 1 * ( conv2(1*~epsGridPad,filter,'valid') >= sum(sum(filter)) );
                
                epsGridPad = zeros(Ny+4*radCurv,Nx+4*radCurv);
                epsGridPad(1+2*radCurv:end-2*radCurv,1+2*radCurv:end-2*radCurv) = epsGrid0;
                epsGridOut = 1 * ( conv2(1*~epsGridPad,filter,'valid') >= sum(sum(filter)) );
                
                maskGridPad = ones(Ny+4*radCurv,Nx+4*radCurv);
                maskGridPad(1+2*radCurv:end-2*radCurv,1+2*radCurv:end-2*radCurv) = maskGrid;
                maskGridPad = 1 * ( conv2(1*maskGridPad,filter,'valid') >= sum(sum(filter)) );
                
                epsGridOut = 1*(epsGridOut.*epsGridOut1);
                
                epsGridOut = 1 * ( conv2(epsGridOut,filter,'same') == 0 );
                epsGridOut = epsGridOut(1+radCurv:end-radCurv,1+radCurv:end-radCurv);
                
                epsGridPad = ones(Ny+4*radCurv,Nx+4*radCurv);
                epsGridPad(1+2*radCurv:end-2*radCurv,1+2*radCurv:end-2*radCurv) = epsGridOut;
                epsGridIn = 1 * ( conv2(epsGridPad,filter,'valid') >= sum(sum(filter)) );
                epsGridIn = 1 * ( conv2(epsGridIn,filter,'valid') ~= 0 );
                
                epsGridNew = epsGridIn;
                if(obj.testFlag)
                    %figure(21); imagesc(epsGrid0); colormap(bluewhitered); axis equal; title('Before Hard RoC Constraint','FontWeight','bold','FontSize',12,'FontName','Arial');
                    figure(22); imagesc(epsGridNew-epsGrid0); colormap(bluewhitered); axis equal; title('Delta from Hard RoC Constraint','FontWeight','bold','FontSize',12,'FontName','Arial');
                    %figure(23); imagesc(epsGridNew); colormap(bluewhitered); axis equal; title('After Hard RoC Constraint','FontWeight','bold','FontSize',12,'FontName','Arial');
                end
            else
                epsGridNew=epsGrid0;
            end
            
        end
        
        %% SOFT MINIMUM DIMENSION CONSTRAINT
        % Zeroes any negative boundary velocity impinging on a feature of
        % minimum dimension
        function velocityMinDim = fixMinDimension(obj, epsGrid0, V)
            Nx = obj.numX;
            Ny = obj.numY;
            minDim = round(obj.minDimension/2); % ie. radius of min dimension
            
            filterSize = 2*minDim + 1; % must be ODD
            filterMid = (filterSize+1)/2;
            
            if(obj.numY==1)
                filterX = 1:filterSize;
                filter = 1* ( (filterX-filterMid).^2 < minDim.^2 );
                epsGridPad = ones(1,Nx+2*minDim);
                epsGridPad(1,1+minDim:end-minDim) = epsGrid0;
                epsGridIn1 = 1 * ( conv(epsGridPad,filter,'valid') >= sum(filter) );
                epsGridIn = 1 * ( conv(epsGridIn1,filter,'same') ~= 0 );
            else
                [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
                filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < minDim.^2 );
                epsGridPad = ones(Ny+2*minDim,Nx+2*minDim);
                epsGridPad(1+minDim:end-minDim,1+minDim:end-minDim) = epsGrid0;
                epsGridIn = 1 * ( conv2(epsGridPad,filter,'valid') >= sum(sum(filter)) );
                epsGridIn = 1 * ( conv2(epsGridIn,filter,'same') ~= 0 );
            end
            
            minDimMask = obj.getBoundaryIn(epsGridIn)+obj.getBoundaryOut(epsGridIn);
            
            velocityMinDim = V.*(minDimMask.*(V<0)+(V>0));
        end
        
        %% SOFT MINIMUM PADDING CONSTRAINT
        % Zeroes any positive boundary velocity impinging on the padding
        % between 2 shapes
        function velocityMinPad = fixMinPadding(obj, epsGrid0, V)
            Nx = obj.numX;
            Ny = obj.numY;
            minDim = round(obj.minPadding/2); % ie. radius of min padding
            
            filterSize = 2*minDim + 1; % must be ODD
            filterMid = (filterSize+1)/2;
            
            if(obj.numY==1)
                filterX = 1:filterSize;
                filter = 1* ( (filterX-filterMid).^2 < minDim.^2 );
                epsGridPad = ones(1,Nx+2*minDim);
                epsGridPad(1,1+minDim:end-minDim) = 1*(~epsGrid0);
                epsGridIn1 = 1 * ( conv(epsGridPad,filter,'valid') >= sum(filter) );
                epsGridIn = 1 * ( conv(epsGridIn1,filter,'same') ~= 0 );
            else
                [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
                filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < minDim.^2 );
                epsGridPad = ones(Ny+2*minDim,Nx+2*minDim);
                epsGridPad(1+minDim:end-minDim,1+minDim:end-minDim) = 1*(~epsGrid0);
                epsGridIn = 1 * ( conv2(epsGridPad,filter,'valid') >= sum(sum(filter)) );
                epsGridIn = 1 * ( conv2(epsGridIn,filter,'same') ~= 0 );
            end
            
            minDimMask = obj.getBoundaryIn(epsGridIn)+obj.getBoundaryOut(epsGridIn);
            
            velocityMinPad = V.*(minDimMask.*(V>0)+(V<0));
        end
        
        %% returnData()
        % return nk data for import to Lumerical
        function data = returnData(obj)  
            % FreeForm Binary Import
            if( (obj.numX>1) && isreal(obj.eps_) && isreal(obj.epsOut) )
                numY = obj.numY;
                numX = obj.numX;
                [xGrid, yGrid] = meshgrid(1:numX, 1:numY);
                epsGrid = obj.epsGrid;
                cnt = 0;
                yArr=[];
                xArr=[];
                blockSizeArr=[];
                
                if(obj.numY==1)
                    bndIn = obj.getBoundaryIn(obj.epsGrid);
                    ind = find(bndIn);
                    x1 = obj.xGrid(ind(1:2:end-1));
                    x2 = obj.xGrid(ind(2:2:end));
                    % eps_ center coordinates and square dimensions
                    xArr = .5*(x2+x1);
                    yArr = obj.yGrid(1,1);
                    blockSizeArr = x2-x1;
                else
                    for blockSize = 150:-2:0 %min(numX,numY):-2:0
                        filterSize = floor(blockSize/2)*2 + 1; % Must be ODD (re-enforce in case for loop is changed)
                        filter = ones(filterSize);
                        epsGridOld = zeros(size(epsGrid));
                        while( any(any(epsGrid - epsGridOld)) )
                            epsGridOld = epsGrid;
                            epsGrid2 = conv2(1*epsGrid,filter,'same') == filterSize^2;
                            if(filterSize > 1)
                                [yInd, xInd] = find(epsGrid2,1,'first');
                                if(~isempty(yInd))
                                    cnt = cnt + 1;
                                    yArr(cnt) = yInd;
                                    xArr(cnt) = xInd;
                                    blockSizeArr(cnt) = filterSize;
                                    dEpsGrid = (abs(xGrid - xInd) <= filterSize/2) & (abs(yGrid - yInd) <= filterSize/2);
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
                        end
                    end
                    % eps_ center coordinates and square dimensions
                    xArr = obj.xGrid(1,xArr);
                    yArr = obj.yGrid(yArr,1);
                    blockSizeArr = blockSizeArr*obj.dx;
                end
                

                % epsOut center coordinates and dimensions
                geo_x = mean([obj.xGrid(1,1) obj.xGrid(1,end)]);
                geo_xspan = (obj.numX-1)*obj.dx;
                geo_y = mean([obj.yGrid(1,1) obj.yGrid(end,1)]);
                geo_yspan = (obj.numY-1)*obj.dx;
                % center z coordinate and extrusion depth
                geo_z = obj.z0;
                geo_zspan = (obj.numZ-1) * obj.dx;
                
                if( ischar(obj.eps_) )
                    geo_n = 0;
                    geo_mat = obj.eps_;
                else
                    geo_n = sqrt(obj.eps_);
                    geo_mat = '';
                end
                if( ischar(obj.epsOut) )
                    geo_nClad = 0;
                    geo_matClad = obj.epsOut;
                else
                    geo_nClad = sqrt(obj.epsOut);
                    geo_matClad = '';
                end
                
                data = struct('xArr',xArr,'yArr',yArr,'blockSizeArr',blockSizeArr,'geo_n',geo_n,'geo_mat',geo_mat,'geo_nClad',geo_nClad,'geo_matClad',geo_matClad,'geo_x',geo_x,'geo_xspan',geo_xspan,'geo_y',geo_y,'geo_yspan',geo_yspan,'geo_z',geo_z,'geo_zspan',geo_zspan);
                
                
            % FreeForm importnk Matrix
            elseif( ~isreal(obj.eps_) && ~isreal(obj.epsOut) && ~ischar(obj.eps_) && ~ischar(obj.epsOut) )
                n = sqrt((obj.epsGrid==1)*obj.eps_ + (obj.epsGrid==0)*obj.epsOut);
                n = n.'; % transpose for Lumerical
                zVec = [obj.zVec(1),obj.zVec(end)];
                len = length(zVec);
                n = repmat(n,[1 1 len]);
                x = obj.xGrid(1,:);
                y = obj.yGrid(:,1);
                z = zVec;
                if(obj.numY==1)
                    n = permute(n,[1 3 2]);
                    data = struct('geo_n',n,'geo_x',x,'geo_y',z,'geo_z',y);
                else
                    data = struct('geo_n',n,'geo_x',x,'geo_y',y,'geo_z',z);
                end
                
            % Error
            else
                error('No export option available for material types.');
            end
            
            %             FreeForm_image Import
            %             if( ischar(obj.eps_) )
            %                 geo_n = 0;
            %                 geo_mat = obj.eps_;
            %             else
            %                 geo_n = sqrt(obj.eps_);
            %                 geo_mat = '';
            %             end
            %
            %             geo_x = obj.x0;
            %             geo_xspan = (obj.numX-1) * obj.dx;
            %             geo_y = obj.y0;
            %             geo_yspan = (obj.numY-1) * obj.dx;
            %             geo_z = obj.z0;
            %             geo_zspan = (obj.numZ-1) * obj.dx;
            %
            %             geo_image = 'geo.jpg';
            %             imwrite(obj.epsGrid,geo_image,'jpg')
            %
            %             data = struct('geo_mat',geo_mat,'geo_n',geo_n,'geo_x',geo_x,'geo_xspan',geo_xspan,'geo_y',geo_y,'geo_yspan',geo_yspan,'geo_z',geo_z,'geo_zspan',geo_zspan,'geo_image',geo_image);
        end
        
        %% getDrawing()
        % return permitivity data for Matlab display to user
        function drawing = getDrawing(obj, x, y, pad)
            drawing = (obj.epsGrid==1)*obj.epsVec(1) + (obj.epsGrid==0)*obj.epsOutVec(1);
        end
        
    end
    
end
