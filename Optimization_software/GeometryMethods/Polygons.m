classdef Polygons
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %   assuming z to be thickness
    % move dx to external parameter section in all shape classes
    % Combine update and acceptUpdates
    
    %% PROPERTIES
    properties
        %External Parameters
        x0;
        y0; % bottom left corner
        z0; % center
        dx; %dx = dy
        thickness;
        epsClad; % single eps value or 'material name'

        %Internal Mesh 
        numX;
        numY;
        zVec;
        
        %Boundary Conditions 0=nothing, +-1=sym, 2=periodic
        xBC;
        yBC;
        
        %Array of shape objects
        shapes;
        numShapes;
        numDynamicShapes;
        
        % E, epsBgnd are cell arrays of size [1 x numFreq x numUser]
        % EA is a cell array of size [numMon x numFreq x numUser]
        E;
        EA; 
        epsBgnd;
        epsCladVec; %vector of actual epsClad values
        newShapeEpsVec; %vector of actual newShapeEps values
        
        % Derivative matrix
        dFdx;
        dFdxOld;
        deltaXOld;
        FOld;
        B;
        
        % These should be set by the Optimizer object
        % Note that numMon is the number of monitors SIMULATIONS, and that
        % numFreq is the total number of frequencies
        numMon;
        numFreq;
        numUD;
        
        %Line data (i.e. boundary data)
        xLine;
        yLine;
        nxLine;
        nyLine;
        VLine;

        % Current Shapes
        currShapePad;
        maxMove; % limits max abs change in a single parameter
        maxArea; % limits total area change
        eraseSize;
        velPadding; % Size of velocity region
        lineStep;
        
        % New shapes
        newShapeCreation;
        newShapeType; % string
        newShapeEps;
        newShapePad;
        newShapeRules;
        newShapeRad; % limits size of newly created shapes
        newShape; % Calculated derivative and location of new shape
        velocity; % unmasked velocity
        velocitym; % masked velocity
        illegalMask;
		
		% Test flag (displaying output if 1)
		testFlag;
        
        quasiNewton;

    end
    
    methods
        %% CONSTRUCTOR
        function obj = Polygons(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, epsClad, xBC, yBC, newShapeCreation, maxMove, maxArea, newShapeRad, eraseSize, velPadding, lineStep, newShapePad, currShapePad, newShapeType, newShapeEps, newShapeRules)
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx))+1;
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx))+1;
            obj.dx = dx;
            obj.thickness = thickness;
            obj.epsClad = epsClad;
            obj.xBC = xBC;
            obj.yBC = yBC;
            obj.newShapeCreation = newShapeCreation;
            obj.maxMove = maxMove / dx;
            obj.maxArea = maxArea / dx^2;
            obj.newShapeRad = newShapeRad / dx;
            obj.eraseSize = floor(1e-4*round(1e4*eraseSize / dx));
            obj.velPadding = floor(1e-4*round(1e4*velPadding / dx));
            obj.numShapes = 0;
            obj.numDynamicShapes = 0;
            obj.shapes = {};
            obj.lineStep = floor(1e-4*round(1e4*lineStep / dx));
            obj.newShapePad  = newShapePad / dx;
            obj.currShapePad = currShapePad / dx;
            obj.newShapeType = newShapeType;
            obj.newShapeEps = newShapeEps;
            obj.newShapeRules = newShapeRules;
            obj.illegalMask = zeros(obj.numY,obj.numX);
            
            obj.quasiNewton = 1;

            if( (obj.numX<=1) && (obj.numY<=1) )
                error('Polygons area is not large enough. Both xLength and yLength must be greater than dx');
            end
        end
        
        %% GET DATA CONSISTENT WITH OPTIMIZER
        function[obj] = setOptData(obj,testFlag,numMon,numFreq,numUD)
            obj.testFlag = testFlag;
            obj.numMon = numMon;
            obj.numFreq = numFreq;
            obj.numUD = numUD;
            obj = obj.initData; % Can initialize data now
        end
        
        %% ADD SHAPE TO GEOMETRY
		% functional form:
        % obj = Rectangle(eps, dynamic, meshOrder, rules, params, testFlag)
        function obj = addShape(obj,newShapeObj)
            if( strcmp('rectangle', newShapeObj{1}) )
				obj.numShapes = obj.numShapes + 1;
                obj.shapes{obj.numShapes} = Rectangle( newShapeObj{2}, newShapeObj{3}, newShapeObj{4}, newShapeObj{5}, newShapeObj{6}, obj.testFlag);
                if(newShapeObj{3}), obj.numDynamicShapes = obj.numDynamicShapes + 1; end;
            elseif( strcmp('ellipse', newShapeObj{1}) || strcmp('circle', newShapeObj{1}) )
                obj.numShapes = obj.numShapes + 1;
                obj.shapes{obj.numShapes} = Ellipse( newShapeObj{2}, newShapeObj{3}, newShapeObj{4}, newShapeObj{5}, newShapeObj{6}, obj.testFlag);
                if(newShapeObj{3}), obj.numDynamicShapes = obj.numDynamicShapes + 1; end;
            else
                error('New shape not understood.');
            end
        end
        
        %% UPDATE FIELD DATA
        function obj = updateFieldData(obj, data, x_vec, y_vec, z_vec, monIndex, freqInd, userIndex, EAflag)
            [x_grid, y_grid] = meshgrid(x_vec,y_vec);
            for i = 1:size(data,1)
                Ei = obj.interpolateData(data(i,:), x_grid, y_grid);
                if(EAflag)
                    obj.EA(monIndex,freqInd(i),userIndex) = {Ei};
                else
                    obj.E(1,freqInd(i),userIndex) = {Ei};
                end
            end
        end
        
        %% UPDATE EPSILON DATA 
        % no user index, assumed constant over user sim's
        function obj = updateEpsData(obj, data, x_vec, y_vec, z_vec, freqInd)
            [x_grid, y_grid] = meshgrid(x_vec,y_vec);
            for i = 1:size(data,1)
                obj.epsBgnd(freqInd(i)) = {obj.interpolateData(data(i,:), x_grid, y_grid)};
            end
        end
        
        %% INIT FIELD DATA
        function[obj] = initData(obj)
            obj.E = cell(1,obj.numFreq, obj.numUD);
            obj.EA = cell(obj.numMon,obj.numFreq, obj.numUD);
            obj.epsBgnd = cell(obj.numFreq, 1);
        end
        
        %% INTERPOLATE FIELDS ONTO GEOMETRY MESH
        % DO NOT AVERAGE FIELDS IN Z, give entire array to the shape
        % Note that this is also used for eps = {epsx,epsy,epsz}
        function field = interpolateData(obj, field, x_grid, y_grid)
            xi = IntToExt(obj,1:obj.numX,ones(1,obj.numX),ones(1,obj.numX));
            yi = IntToExt(obj,1:obj.numY,2*ones(1,obj.numY),ones(1,obj.numY));
            [xInt,yInt] = meshgrid(xi,yi);          
            for z = 1 : size(field{1},3)
                Ex = interp2(x_grid, y_grid, field{1}(:,:,z), xInt, yInt, 'cubic');
                Ey = interp2(x_grid, y_grid, field{2}(:,:,z), xInt, yInt, 'cubic');
                Ez = interp2(x_grid, y_grid, field{3}(:,:,z), xInt, yInt, 'cubic');
                
                Ex(isnan(Ex)) = 0;
                Ey(isnan(Ey)) = 0;
                Ez(isnan(Ez)) = 0;
                Ex_array(:,:,z) = Ex;
                Ey_array(:,:,z) = Ey;
                Ez_array(:,:,z) = Ez;
            end

            field = {Ex_array,Ey_array,Ez_array};
        end
        
        %% UPDATE ALL SHAPES
        function [obj, deltaX, dF, bndUpdate] = updateShapes(obj, deltaX, F0, dFdxOpt)
            
            xNew = obj.newShape(1);
            yNew = obj.newShape(2);
            rNew = obj.newShape(3);
            dANew = pi*rNew^2*obj.dx^2;
            dFNew = obj.newShape(4);
            
            dFdxOpt = permute(dFdxOpt,[2,1]);
            if(obj.quasiNewton && ~isempty(obj.dFdxOld) && length(dFdxOpt)==length(obj.dFdxOld))
                y = (dFdxOpt - obj.dFdxOld);
                s = obj.deltaXOld;
                B = obj.B + (y*y')/(y'*s) - (obj.B*(s*s')*obj.B)/(s'*obj.B*s);
                %B = (dFdxOpt - obj.dFdxOld) ./ (-obj.deltaXOld);
                p = -inv(B)*dFdxOpt;
                p(isnan(p)) = 0;
                p(isinf(p)) = 1;
                %deltaX = (dFdxOpt ./ B .* (B>0)) + (deltaX .* (B<=0));
                if(F0<obj.FOld)
                    obj.maxMove = 0.5*obj.maxMove;
                end
                deltaX = obj.maxMove * p /max(abs(p));
                deltaX(isnan(deltaX)) = 0;
                deltaX(isinf(deltaX)) = 0;
                figure(100); plot(obj.dFdxOld);
                figure(101); plot(obj.deltaXOld);
                figure(102); plot(dFdxOpt);
                figure(103); plot(deltaX);
                %deltaX(deltaX>2*obj.maxMove) = 2*obj.maxMove;
                %deltaX(deltaX<-2*obj.maxMove) = -2*obj.maxMove;
                %norm = max(abs(deltaX));
                %deltaX = deltaX / norm * obj.maxMove;
            else
                B = eye(length(dFdxOpt));
            end
            
            [deltaX,removeFlags] = obj.simulateUpdate(deltaX);
            dFBnd = sum(dFdxOpt.*deltaX);
            dABnd = sum(abs(deltaX))*obj.dx^2;
            if(dABnd==0), dABnd=Inf; end;
            if(isempty(dFBnd)), dFBnd=0; dABnd = Inf; end;
            Fnew = F0 + dFBnd;
            dFBnd = min(Fnew)-min(F0);
            
            % For Testing
            if(obj.testFlag)
                fprintf(' dFNew: %g \n dFBnd: %g \n dFNew/dANew: %g \n dFBnd/dABnd: %g \n',...
                    dFNew, dFBnd, dFNew/dANew, dFBnd/dABnd);
            end
            
            obj.dFdxOld = [];
            obj.deltaXOld = [];
            obj.FOld = F0;
            obj.B = B;
            % Note this assumes rectangle accepts radius, not total length
            if( (obj.newShapeCreation==0) || ( (dFBnd/dABnd > dFNew/dANew) && (dFBnd>dFNew) ) || (dFNew==0) )
                obj = obj.removeShapes(removeFlags~=0);
                obj = obj.acceptUpdates(deltaX);
                dF = dFBnd;
                bndUpdate = 1;
                obj.dFdxOld = dFdxOpt;
                obj.deltaXOld = deltaX;
                fprintf('  Moving Boundary');
            else
                newShapeObj = {obj.newShapeType, obj.newShapeEps, 1, 1, obj.newShapeRules, [xNew, yNew, obj.newShapeRad, obj.newShapeRad]};
                obj = obj.addShape(newShapeObj);
                obj.shapes{end}.epsVec = obj.newShapeEpsVec;
                dF = dFNew;
                bndUpdate = 0;
                fprintf('  New shape: %s, at (%g,%g) \n',obj.newShapeType,xNew,yNew);
            end

        end
        
        %% CALCULATE BOUNDARY DERIVATIVES
        % Note: if numMonSims > 1, minDim = 1
        % dFdx has size [numMonSimsxnumFreqSimsxnumUDx4*numDynamicShapes]
        function[obj,lb,ub,c] = calcDerivs(obj,F0,skipInd,freq)
            if(obj.numDynamicShapes > 0)
                obj.dFdx = zeros(obj.numMon, obj.numFreq, obj.numUD, 4*obj.numDynamicShapes);
                
                thickness = obj.thickness;
                if(thickness==0)
                    thickness = 1;
                end
                            
                % Find derivatives with respect to all shape parameters
                dynamicCnt = 0;
                for j = 1:obj.numShapes
                    if(obj.shapes{j}.dynamic)
                        dynamicCnt = dynamicCnt + 1;
                        [obj.shapes{j},dFdxj] = obj.shapes{j}.calcDerivs(obj.E, obj.EA, obj.epsBgnd, obj.eraseSize, obj.velPadding, obj.lineStep, obj.dx);
                        dFdxj = dFdxj * thickness;
                        obj.dFdx(:,:,:,4*dynamicCnt-3:4*dynamicCnt) = dFdxj;
                    end
                end
            else
                obj.dFdx = zeros(obj.numMon, obj.numFreq, obj.numUD, 1);
            end
            
            % Find velocity for new shape creation
            obj = obj.shapeDeriv(F0,skipInd);
            
            ub = obj.maxMove*ones(4*obj.numDynamicShapes,1);
            lb = -ub;
            c = @(x)( obj.ConstraintFx(x) );
        end
        
        %% CONSTRAINT FUNCTION
        function[c,ceq] = ConstraintFx(obj,x)
            ceq = 0;
            dA = 0;
            for i=1:obj.numShapes
                if(obj.shapes{i}.dynamic)
                    dA = dA + obj.shapes{i}.deltaArea(x);
                end
            end
            c = dA - obj.maxArea; % enforcing that this must be less than or equal to 0
        end
        
        %% CALCULATE SHAPE DERIVATIVE
        % x,y in internal mesh units
        function obj = shapeDeriv(obj,F0,skipInd)
            
            mask = obj.getDrawing(obj.newShapePad,1); % only mesh order 1 objects
            mask(1:obj.newShapePad,:) = 1e6; % Mask out borders
            mask(end-obj.newShapePad+1:end,:) = 1e6;
            mask(:,1:obj.newShapePad) = 1e6;
            mask(:,end-obj.newShapePad+1:end) = 1e6;
            
            szE = size(obj.E{1,1,1}{1});
            szEA = size(obj.EA);
            numMonSims = szEA(1);
            if( skipInd )
                v = zeros(szE(1),szE(2),max(szEA(skipInd),1));
            else
                v = zeros(szE(1),szE(2),1);
            end
            vm = v;
            for m = 1:obj.numMon
                for f = 1:obj.numFreq
                    
                    % Clausius Mossoti factor for a circle
                    cmx = 3*(obj.newShapeEpsVec(f) - obj.epsCladVec(f))./(obj.newShapeEpsVec(f) + obj.epsCladVec(f));
                    cmy = cmx;
                    cmz = (obj.newShapeEpsVec(f) - obj.epsCladVec(f));
                    
                    if(obj.testFlag)
                        fprintf('Clausius-Mossotti factor = [%g,%g,%g] \n',cmx,cmy,cmz);
                    end
                    
                    for u = 1:obj.numUD
                        E = obj.E{1,f,u};
                        EA = obj.EA{m,f,u};
                        if( ~isempty(EA) )
                            
                            vi = 2*real( cmx.* EA{1} .* E{1} ...
                                + cmy.* EA{2} .* E{2} ...
                                + cmz.* EA{3} .* E{3} );
                            vi = mean(vi,3); % Average velocity along extruded z-direction
                            
                            thickness = obj.thickness;
                            if(thickness==0)
                                thickness = 1;
                            end
                            vi = vi * obj.dx^2 * thickness;
                            
                            vim = vi.*(mask==0);
                            
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
                            v(:,:,c) = v(:,:,c) + vi;
                            vm(:,:,c) = vm(:,:,c) + vim;
                            
                            figure(32); imagesc(vm); colormap(bluewhitered); axis equal; title('Velocity: Sparse Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                        end
                    end
                end
            end
            obj.velocity = v;
            obj.velocitym = vm;
            
            dANew = pi*obj.newShapeRad^2;
            F0 = reshape(F0,1,1,[]);
            minF0 = min(F0);
            Fnew = dANew*vm + repmat(F0,[obj.numY,obj.numX]);
            minFnew = min(Fnew,[],3);
            [maxMinFnew,Y] = max(minFnew,[],1);
            [maxMinFnew,x] = max(maxMinFnew);
            y = Y(x);
            dF = maxMinFnew - minF0;
            
            obj.newShape = [x, y, obj.newShapeRad, dF];
        end
        
        %% SIMULATE UPDATE
        function[deltaX, removeFlags] = simulateUpdate(obj,deltaX)
            % Mask out boundary area (according to obj.currShapePad)
            maskZero = zeros(obj.numY, obj.numX);
            bPad = obj.eraseSize+obj.velPadding-1;
            maskZero(1:bPad,:) = 1e6;
            maskZero(end-bPad+1:end,:) = 1e6;
            maskZero(:,1:bPad) = 1e6;
            maskZero(:,end-bPad+1:end) = 1e6;
            
            removeFlags = zeros(1,obj.numShapes);
            i = 1; reDo = 0;
            dynamicCnt = 0;
            while(i<=obj.numShapes)
                if(obj.shapes{i}.dynamic)
                    if(~reDo)
                        dynamicCnt = dynamicCnt + 1;
                    end
                    deltaParams = deltaX(4*dynamicCnt-3:4*dynamicCnt);
                    mask = maskZero;
                    for j = 1:obj.numShapes
                        if( (obj.shapes{j}.meshOrder == 1) && (j~=i) )
                            mask = mask + j * obj.shapes{j}.getDrawing(obj.numX, obj.numY, obj.currShapePad);
                        end
                    end
                    mask = mask .* obj.shapes{i}.getNewDrawing(obj.numX, obj.numY, obj.currShapePad, deltaParams);
                    
                    if( sum(sum(mask)) ~= 0 )
                        if(reDo) %Trying again didn't work, set to 0
                            ind = find(mask);
                            val = mask(ind(1));
                            if(val~=1e6)
                                fprintf('  Collision avoided between shapes: %g, %g \n',i,val);
                            else
                                fprintf('  Collision avoided between shape %g and boundary \n',i);
                            end
                            deltaX(4*dynamicCnt-3:4*dynamicCnt) = 0;
                            reDo = 0;
                        else % Try again once, with half the step
                            reDo = 1;
                            deltaX(4*dynamicCnt-3:4*dynamicCnt) = deltaX(4*dynamicCnt-3:4*dynamicCnt)/2;
                            i = i-1; % re-try
                        end
                    else
                        snew = obj.shapes{i}.accept(deltaX(4*dynamicCnt-3:4*dynamicCnt));
                        removeFlags(i) =  snew.tooSmall(obj.eraseSize + obj.velPadding);
                        reDo = 0;
                        if(obj.testFlag)
                            if(removeFlags(i))
                                fprintf('  Preparing to remove object: %g \n',i);
                            end
                        end
                    end
                end
                
                i = i+1;
            end
        end
        
        %% ACCEPT UPDATES IF PERMISSIBLE
        function [obj] = acceptUpdates(obj,deltaX)            
            dynamicCnt = 0;
            for i = 1:obj.numShapes
                if(obj.shapes{i}.dynamic)
                    if(obj.testFlag)
                        fprintf('  Accepting updates for shape: %g \n',i);
                    end
                    dynamicCnt = dynamicCnt + 1;
                    obj.shapes{i} = obj.shapes{i}.accept(deltaX(4*dynamicCnt-3:4*dynamicCnt));
                end
            end
        end
        
        %% REMOVE SHAPE FROM GEOMETRY
        function obj = removeShapes(obj, indices)
            numRemovedShapes = length(obj.shapes);
            obj.shapes(indices) = [];
            obj.numShapes = length(obj.shapes);
            numRemovedShapes = numRemovedShapes - obj.numShapes;
            obj.numDynamicShapes = obj.numDynamicShapes - numRemovedShapes;
        end
        
        %% DRAW ALL SHAPES
        function drawing = getDrawing(obj, padding, desiredMesh)

            meshOrder = [];
            for i = 1:obj.numShapes
                meshOrder(i) = obj.shapes{i}.meshOrder;
            end
            
            if(nargin>2)
                mCnt = desiredMesh; 
                drawing = zeros(obj.numY, obj.numX);
            else
                mCnt = max(meshOrder):-1:1;
                if(isempty(obj.epsCladVec))
                    drawing = 0*ones(obj.numY, obj.numX);
                else
                    drawing = obj.epsCladVec(1)*ones(obj.numY, obj.numX);
                end
            end
            for mCount = mCnt
                for i = 1:obj.numShapes
                    if(meshOrder(i) == mCount)
                        shapeDrawing = obj.shapes{i}.getDrawing(obj.numX, obj.numY, padding);
                        drawing = drawing .* 1.*(~shapeDrawing);
                        if( isempty(obj.shapes{i}.epsVec) )
                            drawing = drawing + shapeDrawing;
                        else
                            drawing = drawing + obj.shapes{i}.epsVec(1)*shapeDrawing;
                        end
                    end
                end
            end                 
        end
        
        %% CONCATENATE ALL (X,Y,NX,NY) DATA
        function[obj] = updateLineData(obj)
            obj.xLine = [];
            obj.yLine = [];
            obj.nxLine = [];
            obj.nyLine = [];
            
            for i = 1:obj.numShapes
                if(obj.shapes{i}.dynamic)
                    obj.shapes{i} = obj.shapes{i}.padding(obj.eraseSize, obj.velPadding, obj.lineStep);
                    obj.xLine = cat(2, obj.xLine, obj.shapes{i}.xLine);
                    obj.yLine = cat(2, obj.yLine, obj.shapes{i}.yLine);
                    obj.nxLine = cat(2, obj.nxLine, obj.shapes{i}.nxLine);
                    obj.nyLine = cat(2, obj.nyLine, obj.shapes{i}.nyLine);
                end
            end
        end
        
        %% RETURN SHAPE DATA IN DATAMATRIX
        function data = returnData(obj)
            dataMatrix = cell(obj.numShapes,5);
            cnt = 0;
            for i = 1:obj.numShapes
                if(1 || obj.shapes{i}.dynamic)
                    objShape = obj.shapes{i};
                    cnt = cnt + 1;
                    if( strcmp( class(objShape), 'Rectangle' ) )
                        coord = IntToExt(obj, [objShape.x0, objShape.y0, objShape.xLength, objShape.yLength], [1,2,1,2], [1,1,0,0]);
                        dataMatrix(cnt,:) = [{1}, {objShape.eps}, {coord}, {obj.z0}, {obj.thickness}];
                    elseif( strcmp( class(objShape), 'Ellipse' ) )
                        coord = IntToExt(obj, [objShape.x0, objShape.y0, objShape.xr, objShape.yr], [1,2,1,2], [1,1,0,0]);
                        dataMatrix(cnt,:) = [{2}, {objShape.eps}, {coord}, {obj.z0}, {obj.thickness}];
                    end
                end
            end
            data = dataMatrix(1:cnt,:);
        end
    end
    
end
