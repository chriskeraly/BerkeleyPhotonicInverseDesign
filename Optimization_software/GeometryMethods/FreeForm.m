classdef FreeForm
    
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
        
        %Boundary Conditions
        xBC; % 0=nothing, 1=symmetric along x
        yBC; % 0=nothing, 1=symmetric along y
        diagBC; % 0=nothing, 1=symmetric along xy diagonal
        
        %Freeform shape properties
        eps_; % single eps_ value or 'material name'
        epsOut; % single eps_ value or 'material name'
        epsGrid; % freeform eps_ matrix (1 = eps_, 0 = epsOut)
        maskGrid; % 0s indicate non-designable regions
        
        % Current Shapes
        maxMove; % limits growth of shapes
        maxArea; % limits growth of shapes
        quasiNewton = 1; % 1 = allow for quasi-newton optimization
        dFdxOld;
        deltaXOld;
        bOld;
        
        % Constraints
        minDimension;
        minPadding;
        radiusCurv;
        radiusCurvHard;
        newShapePad;
        
        % New shapes
        newShapeCreation;
        newShapeRad; % Size of new shapes
        initialShape; % Use Sparse Perturbation for massive + random-ish first iteration guess

        % Test Flag
        testFlag = 0;
    end
    
    methods

        function obj = FreeForm(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eps_, epsOut, xBC, yBC, diagBC, newShapeCreation, newShapeRad, newShapePad, maxMove, maxArea, minPadding, radiusCurv, radiusCurvHard, minDimension, initialShape)
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.dx = dx;
            
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx));
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx));
            obj.numZ = floor(1e-4*round(1e4*thickness / dx));
            if(obj.numX<1)
                error('x_span must be > dx');
            end
            if(obj.numY<1)
                error('y_span must be > dx');
            end
            if((obj.numZ~=0) && (obj.numZ<1))
                error('For 2D geometries: thickness must be 0. For 3D geometries: thickness must be > dx');
            end
            
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
            obj.diagBC=diagBC;
            
            obj.maxMove = floor(1e-4*round(1e4*maxMove / dx));
            obj.maxArea = floor(1e-4*round(1e4*maxArea / dx^2));
            if(obj.maxMove<1)
                error('maxMove must be > dx');
            end
            
            obj.radiusCurv = floor(1e-4*round(1e4*radiusCurv / dx));
            obj.radiusCurvHard = floor(1e-4*round(1e4*radiusCurvHard / dx));
            obj.minDimension = floor(1e-4*round(1e4*minDimension / dx));
            if(obj.minDimension < 2*obj.radiusCurvHard)
                error('minDimension must be >= 2*radiusCurvHard');
            end
            
            obj.minPadding = floor(1e-4*round(1e4*minPadding / dx));
            if( (obj.minPadding < 2) && (obj.minPadding ~=0) )
                error('minPadding must either be 0 or > 2*dx');
            end
    
            obj.newShapePad  = floor(1e-4*round(1e4*newShapePad / dx));
            obj.newShapeCreation = newShapeCreation;
            obj.newShapeRad = floor(1e-4*round(1e4*newShapeRad / dx));
            if(newShapeCreation && (obj.newShapeRad < obj.radiusCurvHard) )
                error('newShapeRad must be >= radiusCurvHard');
            end
            
            obj.initialShape = initialShape;
        end
        
         %% GET GEOMETRY DATA
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % (x0,y0) is the bottom left corner coordinate (meters)
        % z0 is the center z coordinate (meters)
        % thickness of extruded planar geometry (meters)
        % dx is the uniform spacing of epsGrid (meters)
        % epsGrid is a binary bitmap of permittivity (1 = eps_, 0 = epsOut)
        % eps_ = relative permittivity or 'Lumerical material name'
        % epsOut = relative permittivity or 'Lumerical material name'
        function [epsGrid, eps_, epsOut, x0, y0, z0, dx, thickness] = getGeometry(obj)
            epsGrid = obj.epsGrid;
            eps_ = obj.eps_;
            epsOut = obj.epsOut;
            x0 = obj.x0;
            y0 = obj.y0;
            z0 = obj.z0;
            dx = obj.dx;
            thickness = obj.numZ * dx;
        end
        
        %% UPDATE SHAPES
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % F0 = current Merit Function
        % dFdxBnd = Gradient with respect to a boundary perturbation
        % dFdxSpace = Gradient w.r.t a pertubation anywhere in space
        function [obj, deltaX, dF, bndUpdate] = updateShapes(obj, deltaX, F0, dFdxBnd, dFdxSpace)
            
            [newShape, dF, newShapeGrid, dFGrid] = obj.calcNewShapePerturbation(dFdxSpace, F0);
            [epsGridNew, dFBnd, deltaXNew, dFdxNew, bNew] = obj.calcBoundaryPerturbation(dFdxBnd, deltaX, F0);

            % Compare new shape vs. boundary perturbation
            if(obj.initialShape && (dFGrid>0) )
                dANew = sum(sum( abs(newShapeGrid) ));
                dFNew = dFGrid;
            else
                r = newShape(3);
                if(obj.numY>1)
                    dANew = pi*r^2;
                else
                    dANew = 2*r;
                end
                dFNew = newShape(4);
            end
            dABnd = sum(sum((epsGridNew~=obj.epsGrid)));
            if(dABnd==0), dABnd=Inf; end;
            if(isempty(dFBnd)), dFBnd=0; dABnd = Inf; end;
            Fnew = F0 + dFBnd;
            dFBnd = min(Fnew)-min(F0);
            
            if(obj.testFlag)
                fprintf(' dFNew: %g \n dFBnd: %g \n dFNew/dANew: %g \n dFBnd/dABnd: %g \n',...
                    dFNew, dFBnd, dFNew/dANew, dFBnd/dABnd);
            end
            
            % Is dF(newShape) > dF(movedBoundary)?
            bndChange = ( (dFBnd/dABnd > dFNew/dANew) && (dFBnd > dFNew) );
            
            obj.dFdxOld = [];
            obj.deltaXOld = [];
            obj.bOld = [];
            if(obj.initialShape && (dFGrid>0) )
                % initialShape
                obj.epsGrid = 1*((obj.epsGrid + newShapeGrid)>0);
                dF = dFGrid;
                obj.initialShape = 0;
                bndUpdate = 0;
                fprintf('  Fast Shape using the Sparse Pertubation Approximation \n');
            elseif( (obj.newShapeCreation==0) || bndChange || (dFNew<.01*min(F0)) )
                % Boundary change
                epsGridNew = obj.smoothRadiusCurv(epsGridNew, obj.maskGrid);
                obj.epsGrid = epsGridNew;
                dF = dF + dFBnd;
                bndUpdate = 1;
                obj.bOld = bNew;
                obj.dFdxOld = dFdxNew;
                obj.deltaXOld = deltaXNew;
                fprintf('  Moving Boundary \n');
            else
                % Add new shapes
                x = newShape(1);
                y = newShape(2);
                r = newShape(3);
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
                if(obj.diagBC==1)
                    obj.epsGrid = 1*(obj.epsGrid | obj.epsGrid');
                end
                dF = dFNew;
                bndUpdate = 0;
                fprintf('  New shape: %s, at (%g,%g) \n','circle',x,y);
            end
            
        end
        
        %% CHANGE STEP-SIZE
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % Optimization wil automatically change the step-size
        % factor < 1 for step-size reduction
        % factor >1 for step-size increase
        function obj = changeStepSize(obj, factor)
            obj.maxMove = obj.maxMove * factor;
        end
        
        %% FREEFORM.M PRIVATE FUNCTIONS BELOW....
        
        
        %% RETURN CONSTRAINT PARAMS
        function [ub,lb,c] = getConstraintParams(obj)
            % Upper Bound for every parameter change
            ub = obj.maxMove*ones(numel(obj.epsGrid),1);
            % Lower Bound for every parameter change
            lb = -ub;
            % Constraint functions
            c = @(x)( obj.ConstraintFx(x) );
        end
        
        %% MAX AREA CHANGE CONSTRAINT FUNCTION
        function[c,ceq] = ConstraintFx(obj,x)
            ceq = 0;
            deltaX = abs(x);
            dA = sum(sum(deltaX));
            c = dA - obj.maxArea; % enforcing that this must be less than or equal to 0
        end
        
        %% GET DRAWING
        % return permitivity data for Matlab display to user
        %function drawing = getDrawing(obj, x, y, pad)
        %    drawing = obj.epsGrid;
            %drawing = (obj.epsGrid==1)*obj.epsVec(1) + (obj.epsGrid==0)*obj.epsOutVec(1);
        %end
        
        %% SET GEOMETRY
        % takes a binary matrix, 1s = obj.eps_, 0s = obj.epsOut
        % setGeometry(epsGrid, x_grid, y_grid)
        function obj = setGeometry(obj, epsGrid, x_grid, y_grid)
            if(obj.numY==1)
                epsGrid = interp1(x_grid,1*epsGrid,obj.xGrid);
            else
                epsGrid = interp2(x_grid, y_grid, 1*epsGrid, obj.xGrid, obj.yGrid);
            end
            epsGrid = round(epsGrid)==1;
            epsGrid = obj.smoothRadiusCurv(epsGrid, obj.maskGrid);
            if(obj.xBC==1)
                epsGrid = 1*(epsGrid | fliplr(epsGrid));
            end
            if(obj.yBC==1)
                epsGrid = 1*(epsGrid | flipud(epsGrid));
            end
            if(obj.diagBC==1)
                epsGrid = 1*(epsGrid | epsGrid');
            end
            obj.epsGrid = epsGrid;
        end
        
        %% SET DESIGNABLE REGION
        % takes a binary matrix, 1s = designable, 0s = non-designable
        % setDesignableRegion(maskGrid, x_grid, y_grid)
        function obj = setDesignableRegion(obj, maskGrid, x_grid, y_grid)
            if(obj.numY==1)
                maskGrid = interp1(x_grid,1*maskGrid,obj.xGrid);
            else
                maskGrid = interp2(x_grid, y_grid, 1*maskGrid, obj.xGrid, obj.yGrid);
            end
            maskGrid = round(maskGrid)==1;
            if(obj.xBC==1)
                maskGrid = 1*(maskGrid | fliplr(maskGrid));
            end
            if(obj.yBC==1)
                maskGrid = 1*(maskGrid | flipud(maskGrid));
            end
            if(obj.diagBC==1)
                maskGrid = 1*(maskGrid | maskGrid');
            end
            obj.maskGrid = maskGrid;
        end
        
        %% CALCULATE NEW SHAPE PERTUBATION
        function [newShape, dF, newShapeGrid, dFGrid] = calcNewShapePerturbation(obj, dFdxSpace, F0)

            Nx = obj.numX;
            Ny = obj.numY;
            
            dFdxSpace = dFdxSpace .* obj.maskGrid;
            if(obj.xBC==1)
                dFdxSpace = .5*(dFdxSpace + fliplr(dFdxSpace));
            end
            if(obj.yBC==1)
                dFdxSpace = .5*(dFdxSpace + flipud(dFdxSpace));
            end
            if(obj.diagBC==1)
                dFdxSpace = .5*(dFdxSpace + dFdxSpace');
            end
            
            % Construct mask for adding eps_
            radCurv = obj.newShapePad;
            filterSize = 2*(radCurv-1) + 1; % must be ODD
            filterMid = (filterSize+1)/2;
            [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
            filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < (radCurv-1).^2 );
            epsGridPad = zeros(Ny+2*(radCurv-1),Nx+2*(radCurv-1));
            epsGridPad(1+(radCurv-1):end-(radCurv-1),1+(radCurv-1):end-(radCurv-1)) = obj.epsGrid;
            epsGridMask = 1 * ( conv2(1*~epsGridPad,filter,'valid') >= sum(sum(filter)) );
            
            % Construct mask for adding epsOut
            radCurv = obj.newShapeRad + obj.minDimension;
            filterSize = 2*(radCurv-1) + 1; % must be ODD
            filterMid = (filterSize+1)/2;
            [filterX, filterY] = meshgrid(1:filterSize, 1:filterSize);
            filter = 1* ( ((filterX-filterMid).^2 + (filterY-filterMid).^2) < (radCurv-1).^2 );
            epsGridPad = zeros(Ny+2*(radCurv-1),Nx+2*(radCurv-1));
            epsGridPad(1+(radCurv-1):end-(radCurv-1),1+(radCurv-1):end-(radCurv-1)) = ~obj.epsGrid;
            epsGridOutMask = 1 * ( conv2(1*~epsGridPad,filter,'valid') >= sum(sum(filter)) );
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
            dFdxSpaceM = dFdxSpace.*repmat(epsGridMask,[1 1 size(dFdxSpace,3)]);
            Fnew = dANew*dFdxSpaceM + F0_arr;
            minFnew = min(Fnew,[],3);
            [maxMinFnew,Y] = max(minFnew,[],1);
            [maxMinFnew,x] = max(maxMinFnew);
            y = Y(x);
            dF = maxMinFnew - minF0;
            
            % Consider circle of epsOut
            dFdxSpaceM = dFdxSpace.*repmat(epsGridOutMask,[1 1 size(dFdxSpace,3)]);
            Fnew = -dANew*dFdxSpaceM + F0_arr;
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
            
            if(obj.initialShape)
                % FastShape, add non-circular regions of eps_ and/or epsOut
                dF = sum(dFdxSpace,3);
                dFmax = max(max(dF));
                dFmin = max(max(dF));
                newShapeGrid = 1 * ( (epsGridMask .* (dF>.5*dFmax)) - (epsGridOutMask .* (dF<.5*dFmin)) );
                newShapeGrid = obj.smoothRadiusCurv(newShapeGrid, obj.maskGrid);
                dFGrid = sum(sum( newShapeGrid .* dF ));
            else
                newShapeGrid = [];
                dFGrid = 0;
            end
        end
        
        %% CALCULATE BOUNDARY PERTUBATION
        function [epsGridNew,dFBnd, deltaXNew, dFdxNew, b] = calcBoundaryPerturbation(obj, dFdxBnd, deltaX, F0)
        % Investigate Boundary Perturbation
            dFdxBnd = permute(dFdxBnd,[2,3,1]);
            dFdxBnd = dFdxBnd .* obj.maskGrid;
            if(obj.xBC==1)
                dFdxBnd = .5*(dFdxBnd + fliplr(dFdxBnd));
            end
            if(obj.yBC==1)
                dFdxBnd = .5*(dFdxBnd + flipud(dFdxBnd));
            end
            if(obj.diagBC==1)
                dFdxBnd = .5*(dFdxBnd + dFdxBnd');
            end
            if(obj.quasiNewton==1 && ~isempty(obj.dFdxOld))
                b = (dFdxBnd-obj.dFdxOld)./obj.deltaXOld;
                noDeriv = (b==0) | isnan(b) | isinf(b);
                if( isempty(obj.bOld) )
                    b(noDeriv) = 1;
                else
                    b(noDeriv) = obj.bOld(noDeriv);
                end
                %gradAscent = (b>0)|(obj.deltaXOld==0)|(sign(dFdxBnd)==sign(obj.dFdxOld));
                gradAscent = (b>0)|(sign(dFdxBnd)==sign(obj.dFdxOld));
                
                % Zero-Out velocity to protect minDimension constraint
                dFdxBndTemp = obj.fixMinDimension(obj.epsGrid, dFdxBnd);
                % Zero-Out velocity to protect minPadding constraint
                dFdxBndTemp = obj.fixMinPadding(obj.epsGrid, dFdxBndTemp);
                % Zero-Out velocity to ignore near-optimum boundaries
                dFdxBndTemp = dFdxBndTemp.*gradAscent;
                norm = max(max(abs(dFdxBndTemp))); % Grad-Ascent Norm
                
                deltaX = (~gradAscent).*dFdxBnd./-b + gradAscent.*dFdxBnd*obj.maxMove/norm;
                tooLarge = abs(deltaX)>obj.maxMove;
                deltaX(tooLarge) = sign(deltaX(tooLarge))*obj.maxMove;
                checkSign = all(sign(deltaX)==sign(dFdxBnd));
                if(~checkSign)
                    warning('InverseDesign:FreeForm','quasiNewton failed checkSign');
                end
                if(obj.testFlag)
                    cLim = max(max(max(abs(obj.deltaXOld))));
                    figure(100); imagesc(obj.deltaXOld); caxis([-cLim cLim]); colorbar; title('deltaXOld');
                    cLim = max(max(max(abs(dFdxBnd(:,:,1)))));
                    figure(101); imagesc(obj.dFdxOld(:,:,1)); caxis([-cLim cLim]); colorbar; title('dFdxOld');
                    figure(102); imagesc(dFdxBnd(:,:,1)); caxis([-cLim cLim]); colorbar; title('dFdx');
                    cLim = min(min(min(b)));
                    figure(103); imagesc(b); caxis([cLim 0]); colorbar; title('b');
                    cLim = max(max(max(abs(noDeriv))));
                    figure(104); imagesc(noDeriv); caxis([-cLim cLim]); colorbar; title('noDeriv');
                    cLim = max(max(max(abs(gradAscent)))); 
                    figure(105); imagesc(gradAscent); caxis([-cLim cLim]); colorbar; title('gradAscent');
                    cLim = max(max(max(abs(deltaX))));
                    figure(106); imagesc(deltaX); caxis([-cLim cLim]); colorbar; title('deltaX');
                end
                
                if(obj.testFlag)
                    dFBnd = sum(sum(deltaX.^2.*b.*~gradAscent + dFdxBnd.*deltaX));
                    fprintf(['dFpred Approx2: ' num2str(dFBnd) '\n']);
                end
            else
                b = [];
            end
            [epsGridNew,dFBnd, deltaXNew, dFdxNew] = obj.shapeLoop(deltaX, dFdxBnd, F0);
            if(obj.testFlag)
                fprintf(['dFpred Approx1: ' num2str(dFBnd) '\n']);
            end
            
    end
        
        %% CALCULATE BOUNDARY PERTUBATION - MAIN HELPER
        function [epsGrid,dF,V,Vorig] = shapeLoop(obj, V, dFdxBnd, F0)
            Nx = obj.numX;
            Ny = obj.numY;
            
            epsilon = 1*obj.epsGrid;
            Emask = obj.maskGrid;
            
            dF=0;
            descale = 1;
            numMoves = obj.maxMove;
            
            numMin = length(F0);
            Vorig = dFdxBnd;
            
            V = V .* Emask;
            if(obj.xBC==1)
                V = .5*(V + fliplr(V));
            end
            if(obj.yBC==1)
                V = .5*(V + flipud(V));
            end
            if(obj.diagBC==1)
                V = .5*(V + V');
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
            V = obj.fixMinPadding(epsilon, V);
            
            % Normalize V to normal angle
            % (when theta=pi/4, boundary will naturally move dx*sqrt(2) rather than dx)
            theta = calcTheta(epsilon);
            thetaNorm = abs(cos(theta)) + abs(sin(theta));
            V = V ./thetaNorm .* Emask; % only used for Vmove calc, recalculated inside main loop
            
            % Establish largest velocity (weighted to move by maxMove)
            if(obj.quasiNewton==1 && ~isempty(obj.dFdxOld))
                Vmove = obj.maxMove;
            else
                Vmax = max(max( abs(V.*getBoundaryOut(epsilon) )) );
                V = V/Vmax*obj.maxMove;
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
                %                 theta = calcTheta(epsilon);
                %                 thetaNorm = abs(cos(theta)) + abs(sin(theta));
                %                 Vup = Vup ./thetaNorm .* Emask;
                %                 Vdown = Vdown ./thetaNorm .* Emask;
                
                binIn = ceil(i/2)/(numMoves);
                binOut = ceil(i/2)/(numMoves);
                epsilonRem=getBoundaryIn(epsilon).*(Vdown<(binIn*-Vmove+0.5));
                epsilonAdd=getBoundaryOut(epsilon).*(Vup>(binOut*Vmove-0.5));
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
                    %epsilon = smoothGrid2(epsilon,flag);
                    %epsilon = smoothGrid(epsilon,flag);
                end
                
                epsilonMask = getBoundaryIn(epsilon)+getBoundaryOut(epsilon);
                
                V = obj.shapeLoopVelocity(V).*epsilonMask.*Emask/descale;
                Vorig = obj.shapeLoopVelocity(Vorig).*repmat(epsilonMask,[1 1 numMin]).*repmat(Emask/descale,[1 1 numMin]);
            end
            
            V = round(V);
            
            epsGrid = Emask.*epsilon + (~Emask).*obj.epsGrid;
            
            if(obj.xBC==1)
                epsGrid = 1*(epsGrid | fliplr(epsGrid));
            end
            if(obj.yBC==1)
                epsGrid = 1*(epsGrid | flipud(epsGrid));
            end
            if(obj.diagBC==1)
                epsGrid = 1*(epsGrid | epsGrid');
            end
            
            minF0 = min(F0);
            Fnew = F0 + dF;
            minFnew = min(Fnew);
            dF = minFnew - minF0;
        end
        
        %% CALCULATE BOUNDARY PERTUBATION - SUB HELPER
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
            minDim = round(obj.minDimension/2) + 2; % ie. radius of min dimension
            
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
            
            minDimMask = getBoundaryIn(epsGridIn)+getBoundaryOut(epsGridIn);
            
            velocityMinDim = V.*(minDimMask.*(V<0)+(V>0));
        end
        
        %% SOFT MINIMUM PADDING CONSTRAINT
        % Zeroes any positive boundary velocity impinging on the padding
        % between 2 shapes
        function velocityMinPad = fixMinPadding(obj, epsGrid0, V)
            if(obj.minPadding>0)
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
                minDimMask = getBoundaryIn(epsGridIn)+getBoundaryOut(epsGridIn);
                velocityMinPad = V.*(minDimMask.*(V>0)+(V<0));
            else
                velocityMinPad = V;
            end
        end
        
    end
    
end
