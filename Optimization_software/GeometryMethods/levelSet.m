classdef levelSet
    
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
        
        % Current Shapes
        maxMove; % limits growth of shapes
        maxArea; % limits growth of shapes
        
        
        
        % Constraints      
        radiusCurv;
       
        
        % New shapes
        %level set specifics
        sideMask;
        alpha;
        maxlsIter;
        phi;
        vel;
        iteration;
        sideAnchoring
        
        

        % Test Flag
        testFlag = 1;
    end
    
    methods

        function obj = levelSet(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eps_, epsOut, xBC, yBC, maxArea, radiusCurv,diagBC,maxlsIter,sideAnchoring)
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
            obj.maskGrid=ones(obj.numY,obj.numX);
            obj.xBC = xBC;
            obj.yBC = yBC;
            obj.diagBC=diagBC;

            
            obj.maxArea = floor(1e-4*round(1e4*maxArea / dx^2));
            
            
            obj.radiusCurv = floor(1e-4*round(1e4*radiusCurv / dx));
            
            obj.sideMask=boardercurv(obj.phi,obj.radiusCurv/dx,dx);
            obj.alpha=0.3;
            obj.maxlsIter=maxlsIter;
            obj.iteration=0;
            obj.sideAnchoring=sideAnchoring;
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
        function [obj, dF] = updateShapes(obj, F0, dFdxBnd, dFdxSpace)
            
            dx=obj.dx;
            r0=obj.radiusCurv*dx;
            b=2*r0;
            
            dFdxSpaceraw=dFdxSpace;
            
            % Improvements to make: zero velocity wherever there is no
            % chance any material will ever go
            
            %Give the option of using dFdxBnd
            %New shapes possible....
            
            
            %Include the side masks to the velocity function, and normalize
            %to have the max=1
            %if timestep is tiny error out or propose to relax the
            %curvature constraint to some extent
            
            if (obj.yBC==1)
                dFdxSpace=(dFdxSpace+flipud(dFdxSpace))/2;
            end
            
            if (obj.xBC==1)
                dFdxSpace=(dFdxSpace+fliplr(dFdxSpace))/2;
            end
            
            if obj.diagBC==1
                dFdxSpace=(dFdxSpace+dFdxSpace');
            end
                       
            
            %figure(90);imagesc(obj.phi);
            %%Buffer zone use
            approxSD=messySD(obj.phi,obj.dx);
            %figure(91);imagesc(approxSD);
            bufferzone=(abs(approxSD)<(obj.dx*10));dFdxSpace2=dFdxSpace.*(obj.sideMask==0);
            insidemax=max(max(abs(dFdxSpace2.*bufferzone)));
            outsidemax=max(max(abs(dFdxSpace2)));
            velocity=(dFdxSpace2.*bufferzone)/insidemax+dFdxSpace2.*(1-bufferzone)/outsidemax;
            
            
            %%End buffer zone use  put this back if not using:velocity=dFdxSpace/max(max(abs(dFdxSpace.*(obj.sideMask==0))));
  
            
            obj.vel=velocity;
            velocity=velocity.*(obj.sideMask==0)+obj.sideMask;
            figure(56);imagesc(velocity);colormap('bluewhitered');
            dt=obj.alpha/(2*max(max(abs(velocity)))/dx+4*b/dx^2);
            
            fprintf('dt/dx is %g \n',dt/dx);
            
            finished=0;
            iter=0;
            newphi=obj.phi;
            
            while finished==0 && iter<obj.maxlsIter
                
                
                velDotPhi=vDotPhi(newphi,velocity,dx);
                curv=(abs(curvterm(newphi,dx))>1/r0)*1/r0.*sign(curvterm(newphi,dx));
                gradphi=sqrt(FODiff(newphi,'x','central',dx).^2+FODiff(newphi,'y','central',dx).^2);
                
                newphi=newphi-dt*(velDotPhi-b*gradphi.*curv);
                
                if floor(iter/50)==iter/50
                    newphi=signedDist(newphi,dx);
                    fprintf('.');
                end
                iter=iter+1;
                
                movedArea=sum(sum(abs(sign(newphi)-sign(obj.phi))));                
                if movedArea>obj.maxArea
                    finished=1;
                end
                
            end
            
            if (obj.yBC==1)
                newphi=(newphi+flipud(newphi))/2;
            end
            
            if (obj.xBC==1)
                newphi=(newphi+fliplr(newphi))/2;
            end
            
            if obj.diagBC==1
                newphi=(newphi+newphi');
            end
            
            fprintf('%g level set updates required  \n',iter);
            fprintf('%g pixels changed \n',movedArea);
            
            if movedArea<0.01*obj.maxArea
                error('shape is no longer changing by more than 1% of maxArea, consider increasing maxlsIter')
            end
            
            dF=-sum(sum(((sign(newphi)-sign(obj.phi)).*dFdxSpaceraw)));
            
            obj.phi=newphi;
            
            obj.iteration=obj.iteration+1;
            %%Complete reinitilization 
            if obj.iteration/10==floor(obj.iteration/10)
                obj.phi=messySD(obj.phi,obj.dx);
            end
            %%end reinit
    
            obj.epsGrid=abs((sign(obj.phi)-1)/2);
            
               
            deltaX=0; bndUpdate=1;
            
        end
       
        %% CHANGE STEP-SIZE
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        function obj = changeStepSize(obj, factor)
            obj.maxArea = obj.maxArea * factor;
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
        
 
        
        %% SET GEOMETRY
        % takes a binary matrix, 1s = obj.eps_, 0s = obj.epsOut
        % setGeometry(epsGrid, x_grid, y_grid)
        function obj = setGeometry(obj, phi)       
            obj.phi=messySD(phi,obj.dx);
            if obj.sideAnchoring==1
                obj.sideMask=boardercurv(obj.phi,obj.radiusCurv*obj.dx,obj.dx);
            else
                obj.sideMask=zeros(size(obj.xGrid));
            end
            obj.epsGrid=abs((sign(obj.phi)-1)/2);
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
