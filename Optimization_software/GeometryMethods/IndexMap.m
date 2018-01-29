classdef IndexMap
    
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
        diagBC;% % 0=nothing, 1=symmetric along diagonal 
        
        %Freeform shape properties
        eps_; % single eps_ value or 'material name'
        epsOut; % single eps_ value or 'material name'
        epsGrid; % freeform eps_ matrix (1 = eps_, 0 = epsOut)
        maskGrid; % 0s indicate non-designable regions
        timeStep;
        
        vel;
        
        % Test Flag
        testFlag = 1;
    end
    
    methods

        function obj = IndexMap(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eps_, epsOut, xBC, yBC, diagBC,timeStep)
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.dx = dx;
            
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx));
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx));
            obj.numZ = floor(1e-4*round(1e4*thickness / dx));
            if(obj.numX<1)
                error('xLengthReal must be > dx');
            end
            if(obj.numX<0)
                error('numY must be > 0');
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
            
            obj.xBC = xBC;
            obj.yBC = yBC;
            obj.diagBC=diagBC;
            obj.timeStep=timeStep;

        end
        
        %% GET GEOMETRY DATA
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % (x0,y0) is the bottom left corner (meters)
        % z0 is the center z coordinate (meters)
        % thickness of extruded planar geometry (meters)
        % dx is the uniform spacing of epsGrid (meters)
        % epsGrid is a binary bitmap of the shapes (1 = eps_, 0 = epsOut)
        % eps_ = dielectric constant or 'Lumerical material name'
        % epsOut = dielectric constant or 'Lumerical material name'
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
        function [obj, dF] = updateShapes(obj,  F0, dFdxBnd, dFdxSpace)
            dFdxSpace0=dFdxSpace;
            figure(345);imagesc(dFdxSpace);
            epsGridOld=obj.epsGrid;
              
            if (obj.yBC==1)
                dFdxSpace=(dFdxSpace+flipud(dFdxSpace))/2;
            end
            
            if (obj.xBC==1)
                dFdxSpace=(dFdxSpace+fliplr(dFdxSpace))/2;
            end
            
            if obj.diagBC==1
                dFdxSpace=(dFdxSpace+dFdxSpace');
            end
            
            obj.epsGrid = obj.epsGrid + (dFdxSpace / max(max(abs(dFdxSpace))) * obj.timeStep);
                               
            dF = sum(sum( (obj.epsGrid-epsGridOld).*dFdxSpace0 ));
            
        end
       
        %% CHANGE STEP-SIZE
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        function obj = changeStepSize(obj, factor)
            obj.timeStep = obj.timeStep * factor;
        end
        
        %% RETURN CONSTRAINT PARAMS
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        function [ub,lb,c] = getConstraintParams(obj)
            % Upper Bound for every parameter change
            ub = obj.maxMove*ones(numel(obj.epsGrid),1);
            % Lower Bound for every parameter change
            lb = -ub;
            % Constraint functions
            c = @(x)( obj.ConstraintFx(x) );
        end
        
        %% MAX AREA CHANGE CONSTRAINT FUNCTION
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        function[c,ceq] = ConstraintFx(obj,x)
            ceq = 0;
            deltaX = abs(x);
            dA = sum(sum(deltaX));
            c = dA - obj.maxArea; % enforcing that this must be less than or equal to 0
        end
        
        %% GET DRAWING
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % return permitivity data for Matlab display to user
        function drawing = getDrawing(obj, x, y, pad)
            drawing = obj.epsGrid;
            %drawing = (obj.epsGrid==1)*obj.epsVec(1) + (obj.epsGrid==0)*obj.epsOutVec(1);
        end
        
        %% SET GEOMETRY
        % takes a binary matrix, 1s = obj.eps_, 0s = obj.epsOut
        % setGeometry(epsGrid, x_grid, y_grid)
        function obj = setGeometry(obj, epsGrid)       
            obj.epsGrid=epsGrid;
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
            obj.maskGrid = maskGrid;
        end
        
       
    end
    
end
