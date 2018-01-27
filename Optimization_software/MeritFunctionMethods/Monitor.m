classdef Monitor
    
    %% PROPERTIES
    properties (Constant)
        mu0 = pi*4e-7;
        c = 299792458;
        eps0 = 1 / (Monitor.mu0*Monitor.c^2);
    end
    
    properties
        % To the external grid (minimum x,y,z values resp.)
        x0; y0; z0;
        dx; % Multiplied by numX & numY gives real lengths
        
        % Internal mesh
        numX;
        numY;
        numZ;
        dim;
        
        % Merit properties
        weight;
        exponent;
        
        % Parameters
        
        
        % testFlag
        testFlag;
    end
    
    %% METHODS
    methods
        
        %% CONSTRUCTOR
        function obj = Monitor(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, testFlag)
            obj.weight = weight;
            obj.exponent = exponent;
            if(exponent==0)
                error('Attempted to create a Monitor with exponent = 0.  Not allowed.\n');
            end
            
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.dx = dx;
            
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx)) + 1;
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx)) + 1;
            obj.numZ = floor(1e-4*round(1e4*zLengthReal / dx)) + 1;
            dimVec = find( [obj.numX-1,obj.numY-1,obj.numZ-1] );
            obj.dim = length( dimVec );
            
            obj.testFlag = testFlag;
        end
        
        %% INTERPOLATE FIELDS ON GRID
        % grid cell contains 3 [Nx1] arrays
        % Ensure order is (y,x,z) for 2d or 3d data
        function [eInt,hInt] = interpFields(obj,eField,hField,grid)
            xInt = obj.x0 + obj.dx *(0:obj.numX-1);
            yInt = obj.y0 + obj.dx *(0:obj.numY-1);
            zInt = obj.z0 + obj.dx *(0:obj.numZ-1);
            gridX = grid{1}; gridY = grid{2}; gridZ = grid{3};
            
            exInt = [];
            eyInt = [];
            ezInt = [];
            hxInt = [];
            hyInt = [];
            hzInt = [];
            
            % 0D
            if(obj.dim == 0)
                if(~isempty(eField))
                    exInt = eField{1};
                    eyInt = eField{2};
                    ezInt = eField{3};
                end
                if(~isempty(hField))
                    hxInt = hField{1};
                    hyInt = hField{2};
                    hzInt = hField{3};
                end
                
                % 1D
            elseif(obj.dim == 1)
                if( length(gridX)>1 ) % x axis
                    gridI = gridX; intI = xInt;
                elseif( length(gridY)>1 ) % y axis
                    gridI = gridY; intI = yInt;
                else % z axis
                    gridI = gridY; intI = zInt;
                end
                if(~isempty(eField))
                    exInt = interp1(gridI,eField{1}(:),intI,'pchip');
                    eyInt = interp1(gridI,eField{2}(:),intI,'pchip');
                    ezInt = interp1(gridI,eField{3}(:),intI,'pchip');
                end
                if(~isempty(hField))
                    hxInt = interp1(gridI,hField{1}(:),intI,'pchip');
                    hyInt = interp1(gridI,hField{2}(:),intI,'pchip');
                    hzInt = interp1(gridI,hField{3}(:),intI,'pchip');
                end
                % 2D
            elseif(obj.dim == 2)
                if( length(gridX)==1 ) % y-z plane, size(data) = [y,1,z]
                    x1 = zInt;
                    x2 = yInt;
                    grid1 = gridZ;
                    grid2 = gridY;
                elseif( length(gridY)==1 ) % x-z plane, size(data) = [1,x,z]
                    x1 = zInt;
                    x2 = xInt;
                    grid1 = gridZ;
                    grid2 = gridX;
                elseif( length(gridZ)==1 ) % x-y plane, size(data) = [y,x,1]
                    x1 = xInt;
                    x2 = yInt;
                    grid1 = gridX;
                    grid2 = gridY;
                end
                if(~isempty(eField))
                    exInt = interp2(grid1,grid2',squeeze(eField{1}),x1,x2','spline');
                    eyInt = interp2(grid1,grid2',squeeze(eField{2}),x1,x2','spline');
                    ezInt = interp2(grid1,grid2',squeeze(eField{3}),x1,x2','spline');
                end
                if(~isempty(hField))
                    hxInt = interp2(grid1,grid2',squeeze(hField{1}),x1,x2','spline');
                    hyInt = interp2(grid1,grid2',squeeze(hField{2}),x1,x2','spline');
                    hzInt = interp2(grid1,grid2',squeeze(hField{3}),x1,x2','spline');
                end
                % Do not use 2D Cubic, often introduces a row of NaNs causing an
                % asymmetry in Adjoint dipoles
                
                % 3D
            elseif(obj.dim==3)
                if(~isempty(eField))
                    exInt = interp3(gridX,gridY',gridZ,eField{1},xInt,yInt',zInt,'spline');
                    eyInt = interp3(gridX,gridY',gridZ,eField{2},xInt,yInt',zInt,'spline');
                    ezInt = interp3(gridX,gridY',gridZ,eField{3},xInt,yInt',zInt,'spline');
                end
                if(~isempty(hField))
                    hxInt = interp3(gridX,gridY',gridZ,hField{1},xInt,yInt',zInt,'spline');
                    hyInt = interp3(gridX,gridY',gridZ,hField{2},xInt,yInt',zInt,'spline');
                    hzInt = interp3(gridX,gridY',gridZ,hField{3},xInt,yInt',zInt,'spline');
                end
            end
            
            exInt(isnan(exInt)) = 0;
            eyInt(isnan(eyInt)) = 0;
            ezInt(isnan(ezInt)) = 0;
            hxInt(isnan(hxInt)) = 0;
            hyInt(isnan(hyInt)) = 0;
            hzInt(isnan(hzInt)) = 0;
            eInt = {exInt,eyInt,ezInt};
            hInt = {hxInt,hyInt,hzInt};
        end
        
        %% PERMUTE TO LUMERICAL FORMAT
        % data3d should be in proper [x,y,z] dimensions for Lumerical
        function [data3d] = permuteLum(obj,data)
            data = squeeze(data); % get rid of extra dimensions
            if(obj.dim==1) % 1D
                if( obj.numX>1 ) % x-dimension only
                    data3d = permute(data,[1 2 3]);
                elseif( obj.numY>1 ) % y-dimension only
                    data3d = permute(data,[2 1 3]);
                else % z-dimension only
                    data3d = permute(data,[2 3 1]);
                end
            elseif(obj.dim ==2) % 2D
                if( obj.numX==1 ) % y-z plane, size(data) = [y,1,z]
                    data3d = permute(data, [3 1 2]);
                elseif( obj.numY==1 ) % x-z plane, size(data) = [1,x,z]
                    data3d = permute(data, [1 3 2]);
                else % x-y plane, size(data) = [y,x,1]
                    data3d = permute(data, [2 1 3]);
                end
            else % 3D, size(data) = [y,x,z]
                data3d = permute(data,[2 1 3]);
            end
        end
    end
end
