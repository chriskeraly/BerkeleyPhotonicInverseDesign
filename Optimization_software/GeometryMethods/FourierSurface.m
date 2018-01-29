classdef FourierSurface
    
    %  Notes
    %  Generalize to 2D Shapes (1D Fourier)
    
    
    properties
        %External Parameters
        x0;
        y0; % bottom left corner
        z0; % center
        dx; %dx = dy
        xVec;
        yVec;
        zVec;
        
        %Internal Mesh
        numX;
        numY;
        numZ;
        bndZ; % average height of boundary
        
        %Boundary Conditions 0=nothing, 1 = symmetry about center
        xBC;
        yBC;
        
        %Fourier Surface properties
        epsAbove;
        epsBelow;
        numCoeffx;
        numCoeffy;
        coeff; % size(coeff) =[numCoeffx, numCoeffy, 4]
        % 3d: (:,:,1:4)
        % -> cos(x)cos(y), sin(x)cos(y), cos(x)sin(y), sin(x)sin(y)
        % 2d: (:,:,1:2)
        % -> cos(x), sin(x)
        surface; % size(surface) =[numX, numY]
        
        %Data = cell(numMon, numFreq, numUD) [numY, numX]
        E;
        EA;
        epsBgnd;
        epsBelowVec;
        epsAboveVec;
        
        % Derivative matrix
        % size = [numMon, numFreq, numUD, numY, numX]
        dFdx;
        dFdx_bnd;
        
        % These should be set by the Optimizer object
        % Note that numMon is the number of monitors SIMULATIONS, and that
        % numFreq is the total number of frequencies
        numMon;
        numFreq;
        numUD;
        
        % Current Shapes
        maxMove; % limits growth of shapes
        maxArea; % limits growth of shapes
        eraseSize; % size of ignored field in velocity calcualtion
        velPadding; % size of region in velocity calculation
        
        % Constraints
        fixedAvgHeight;
        
        % Test flag (displaying output if 1)
        testFlag;
        
    end
    
    methods
        function obj = FourierSurface(x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, numCoeffx, numCoeffy, boundaryZReal, epsAbove, epsBelow, xBC, yBC, maxMove, maxArea, eraseSize, velPadding, fixedAvgHeight)
            obj.x0 = x0;
            obj.y0 = y0;
            obj.z0 = z0;
            obj.dx = dx;
            
            obj.numX = floor(1e-4*round(1e4*xLengthReal / dx))+1;
            obj.numY = floor(1e-4*round(1e4*yLengthReal / dx))+1;
            obj.numZ = floor(1e-4*round(1e4*zLengthReal / dx))+1;
            obj.bndZ = floor(1e-4*round(1e4*boundaryZReal / dx));
            obj.xVec = x0+dx*(0:obj.numX-1);
            obj.yVec = y0+dx*(0:obj.numY-1);
            obj.zVec = z0+dx*(0:obj.numZ-1);
            
            obj.numCoeffx = numCoeffx;
            obj.numCoeffy = numCoeffy;
            obj.epsAbove = epsAbove;
            obj.epsBelow = epsBelow;
            
            obj.xBC = xBC;
            obj.yBC = yBC;
            
            obj.eraseSize = floor(1e-4*round(1e4*eraseSize / dx));
            obj.velPadding = floor(1e-4*round(1e4*velPadding / dx));
            obj.maxMove = floor(1e-4*round(1e4*maxMove / dx));
            obj.maxArea = floor(1e-4*round(1e4*maxArea / dx));
            
            obj.fixedAvgHeight = fixedAvgHeight;
            
            obj = obj.setCoeff(zeros(numCoeffx, numCoeffy, 4));
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
        function obj = setCoeff(obj, coeff)
            obj.coeff = coeff / obj.dx;
            obj.surface = obj.calcSurface(obj.coeff);
        end
        
        %% UPDATE FIELD DATA
        % updateFieldData(data, x_grid, y_grid, monIndex, freqInd, userIndex, EAflag)
        function obj = updateFieldData(obj, data, x_vec, y_vec, z_vec, monIndex, freqInd, userIndex, EAflag)
            for i = 1:size(data,1)
                if(obj.numY==1)
                    % incoming y-data is our z-data (ie. in 'thick' direction)
                    % Swapping Ey and Ez and Swapping y-axis and z-axis for each Ex/Ey/Ez
                    data_i{1} = permute(data{i,1}, [3 2 1]);
                    data_i{2} = permute(data{i,3}, [3 2 1]);
                    data_i{3} = permute(data{i,2}, [3 2 1]);
                    Ei = obj.interpolateData(data_i, x_vec, z_vec, y_vec); % Also swap y_vec/z_vec
                else
                    Ei = obj.interpolateData(data(i,:), x_vec, y_vec, z_vec);
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
                    data_i{1} = permute(data{i,1}, [3 2 1]);
                    data_i{2} = permute(data{i,3}, [3 2 1]);
                    data_i{3} = permute(data{i,2}, [3 2 1]);
                    obj.epsBgnd(freqInd(i)) = {obj.interpolateData(data_i, x_vec, z_vec, y_vec)}; % Also swap y_vec/z_vec
                else
                    obj.epsBgnd(freqInd(i)) = {obj.interpolateData(data(i,:), x_vec, y_vec, z_vec)};
                end
                
                epsB = [obj.epsBgnd{freqInd(i)}{1}(1,1),obj.epsBgnd{freqInd(i)}{2}(1,1,1),obj.epsBgnd{freqInd(i)}{3}(1,1,1)];
                epsA = [obj.epsBgnd{freqInd(i)}{1}(1,1,end),obj.epsBgnd{freqInd(i)}{2}(1,1,end),obj.epsBgnd{freqInd(i)}{3}(1,1,end)];
                epsB = mean(epsB(epsB~=0));
                epsA = mean(epsA(epsA~=0));
                obj.epsBelowVec(freqInd(i)) = epsB;
                obj.epsAboveVec(freqInd(i)) = epsA;
                
                if( isnan(epsB) || isnan(epsA) || (epsB==0) || (epsA==0) )
                    error('BAD epsilon data from Lumerical. Ensure that the Index monitor encompasses the entire Velocity region');
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
        
        %% CALCULATE SURFACE HEIGHT
        function[surface] = calcSurface(obj, coeff)
            surface = zeros(obj.numX,obj.numY);
            [y,x] = meshgrid(0:obj.numY-1, 0:obj.numX-1);
            
            for mh=1:obj.numCoeffx
                m = mh - 1;
                for nh=1:obj.numCoeffy
                    n = nh - 1;
                    surface = surface +  coeff(mh,nh,1)*cos(m*pi*x/(obj.numX/2)).*cos(n*pi*y/(obj.numY/2));
                    surface = surface +  coeff(mh,nh,2)*sin(m*pi*x/(obj.numX/2)).*cos(n*pi*y/(obj.numY/2));
                    surface = surface +  coeff(mh,nh,3)*cos(m*pi*x/(obj.numX/2)).*sin(n*pi*y/(obj.numY/2));
                    surface = surface +  coeff(mh,nh,4)*sin(m*pi*x/(obj.numX/2)).*sin(n*pi*y/(obj.numY/2));
                end
            end
            %surface = round(surface);
        end
        
        %% INTERPOLATE FIELDS ONTO GEOMETRY MESH
        % DO NOT AVERAGE FIELDS IN Z, give entire array to the shape
        % Note that this is also used for eps_ = {epsx,epsy,epsz}
        % interpolateData(field, x_grid, y_grid)
        function [E] = interpolateData(obj, field, x_vec, y_vec, z_vec)
            xInt = obj.xVec;
            yInt = obj.yVec;
            zInt = obj.zVec;
            Ex = zeros(obj.numX,obj.numY,obj.numZ);
            Ey = zeros(obj.numX,obj.numY,obj.numZ);
            Ez = zeros(obj.numX,obj.numY,obj.numZ);
            
            %zIndMax = max(max(obj.surface)) + obj.eraseSize + obj.velPadding + 2;
            %zIndMin = min(min(obj.surface)) - obj.eraseSize - obj.velPadding - 2;
            %zIndMax = min(zIndMax,obj.numZ);
            %zIndMin = max(zIndMin,1);
            %zInt = zInt(zIndMin:zIndMax);
            
            if(obj.numY==1)
                Exi = interp2(z_vec, x_vec', squeeze(field{1}), zInt, xInt', 'spline');
                Eyi = interp2(z_vec, x_vec', squeeze(field{2}), zInt, xInt', 'spline');
                Ezi = interp2(z_vec, x_vec', squeeze(field{3}), zInt, xInt', 'spline');
                Ex = permute(Exi, [1 3 2]);
                Ey = permute(Eyi, [1 3 2]);
                Ez = permute(Ezi, [1 3 2]);
            else
                Ex = interp3(y_vec, x_vec', z_vec, field{1}, yInt, xInt', zInt, 'spline');
                Ey = interp3(y_vec, x_vec', z_vec, field{2}, yInt, xInt', zInt, 'spline');
                Ez = interp3(y_vec, x_vec', z_vec, field{3}, yInt, xInt', zInt, 'spline');
            end
            
            Ex(isnan(Ex)) = 0;
            Ey(isnan(Ey)) = 0;
            Ez(isnan(Ez)) = 0;
            
            E = {Ex,Ey,Ez};
        end

        
        %% UPDATE SHAPES
        % updateShapes()
        function [obj, deltaX, dF, bndUpdate] = updateShapes(obj, deltaX, F0, dFdxOpt)
            if(obj.numY == 1)
                figure(32); plot(obj.calcSurface(deltaX) * obj.dx);
                title('Boundary Movement - Real Space','FontWeight','bold','FontSize',12,'FontName','Arial');
                xlabel('x (m)','FontWeight','bold','FontSize',12,'FontName','Arial');
                ylabel('deltaZ (m)','FontWeight','bold','FontSize',12,'FontName','Arial');
                
                figure(33); plot(deltaX(:) * obj.dx);
                title('Boundary Movement - Fourier Space','FontWeight','bold','FontSize',12,'FontName','Arial');
                xlabel('Cx','FontWeight','bold','FontSize',12,'FontName','Arial');
                ylabel('deltaZ (m)','FontWeight','bold','FontSize',12,'FontName','Arial');
            else
                figure(32); imagesc(obj.calcSurface(deltaX).' * obj.dx);
                title('Boundary Movement - Real Space','FontWeight','bold','FontSize',12,'FontName','Arial');
                xlabel('x (m)','FontWeight','bold','FontSize',12,'FontName','Arial');
                ylabel('y (m)','FontWeight','bold','FontSize',12,'FontName','Arial');
            end
            
            dMax = max(max(max(deltaX)));
            dMin = min(min(min(deltaX)));
            dAbsMax = max(max(max(abs(deltaX))));
            deltaX = deltaX / dAbsMax * obj.maxMove;
            obj.coeff = obj.coeff + deltaX;
            obj.surface = obj.calcSurface(obj.coeff);
                        
            dFdxOpt = permute(dFdxOpt, [2 3 4 1]); % Move min-FOM dimension to last
            deltaX = repmat(deltaX, [1 1 1 size(dFdxOpt,4)]);
            dF = sum(sum(sum(sum(deltaX.*dFdxOpt))));
            
            
            bndUpdate = 1;
            %fprintf('  Moving Boundary');
        end
        
        %% CALCULATE DERIVATIVES
        % Note: if numMonSims > 1, minDim = 1
        % dFdx has size [numMon, numFreq, numUD, numY, numX]
        function[obj,lb,ub,c] = calcDerivs(obj, F0, skipInd)
            
            % Construct padding matrices
            [padTop, padBot] = obj.constructPadding;
            pad = padTop + padBot;
            
            figure(30); clf;
            figure(31); clf;
            
            vLim=0;
            obj.dFdx = zeros(obj.numMon,obj.numFreq, obj.numUD, obj.numCoeffx, obj.numCoeffy, 4);
            obj.dFdx_bnd = zeros(obj.numMon,obj.numFreq, obj.numUD, obj.numX, obj.numY);
            if(any(any(any(pad)))) % only calc dFdx if there exists a boundary
                for m = 1:obj.numMon
                    for f = 1:obj.numFreq
                        for u = 1:obj.numUD
                            E = obj.E{1,f,u};
                            EA = obj.EA{m,f,u};
                            if(~isempty(EA))
                                [E_perp_x,E_perp_y,E_perp_z,E_par_x,E_par_y,E_par_z] = obj.calcParPerpField(E, pad);
                                [EA_perp_x,EA_perp_y,EA_perp_z,EA_par_x,EA_par_y,EA_par_z] = obj.calcParPerpField(EA, pad);
                                
                                D_perp_x=padTop.*obj.epsAboveVec(f).*E_perp_x + padBot.*obj.epsBelowVec(f).*E_perp_x;
                                D_perp_y=padTop.*obj.epsAboveVec(f).*E_perp_y + padBot.*obj.epsBelowVec(f).*E_perp_y;
                                D_perp_z=padTop.*obj.epsAboveVec(f).*E_perp_z + padBot.*obj.epsBelowVec(f).*E_perp_z;
                                
                                DA_perp_x=padTop.*obj.epsAboveVec(f).*EA_perp_x + padBot.*obj.epsBelowVec(f).*EA_perp_x;
                                DA_perp_y=padTop.*obj.epsAboveVec(f).*EA_perp_y + padBot.*obj.epsBelowVec(f).*EA_perp_y;
                                DA_perp_z=padTop.*obj.epsAboveVec(f).*EA_perp_z + padBot.*obj.epsBelowVec(f).*EA_perp_z;
                                
%                                     figure(11); imagesc(angle(squeeze(E_par_x)).'); axis equal;
%                                     figure(12); imagesc(angle(squeeze(E_par_y)).'); axis equal;
%                                     figure(13); imagesc(angle(squeeze(E_par_z)).'); axis equal;
%                                     figure(14); imagesc(angle(squeeze(D_perp_x)).'); axis equal;
%                                     figure(15); imagesc(angle(squeeze(D_perp_y)).'); axis equal;
%                                     figure(16); imagesc(angle(squeeze(D_perp_z)).'); axis equal;
                                
                                dFdx_mfu_p = 2*real( (obj.epsBelowVec(f)-obj.epsAboveVec(f)).*(E_par_x.*EA_par_x + E_par_y.*EA_par_y + E_par_z.*EA_par_z) );
                                dFdx_mfu_n = 2*real( (1/obj.epsAboveVec(f)-1/obj.epsBelowVec(f)).*(D_perp_x.*DA_perp_x + D_perp_y.*DA_perp_y + D_perp_z.*DA_perp_z) );
                                dFdx_mfu_pad = dFdx_mfu_p + dFdx_mfu_n;
                                
                                if(obj.testFlag)
                                    vLim = max(max(max( abs(dFdx_mfu_pad) )));
                                    figure(17); imagesc(squeeze(dFdx_mfu_p).'); caxis([-vLim vLim]); colormap('bluewhitered'); axis equal;
                                    figure(18); imagesc(squeeze(dFdx_mfu_n).'); caxis([-vLim vLim]); colormap('bluewhitered'); axis equal;
                                    figure(19); imagesc(squeeze(dFdx_mfu_pad).'); caxis([-vLim vLim]); colormap('bluewhitered'); axis equal;
                                end
                                
                                
                                % Interpolate dFdx_mfu only on boundary
                                dFdx_mfu_bnd = obj.interpolateVelocity(dFdx_mfu_pad, pad);
                                dFdx_mfu_bnd = dFdx_mfu_bnd * obj.dx^(2 + 1*(obj.numY>1));
                                
                                % Zero-Out Velocity if surface will punch-through
                                mask = (obj.surface < (obj.numZ-1.5*obj.maxArea-1)) & (obj.surface > (1.5*obj.maxArea+1));
                                dFdx_mfu_bnd = dFdx_mfu_bnd .* mask;
                                
                                if(obj.numY==1)
                                    figure(30); hold all;
                                    title('Velocity: Real Space','FontWeight','bold','FontSize',12,'FontName','Arial');
                                    xlabel('x (m)');
                                    
                                    plot(squeeze(dFdx_mfu_bnd).');
                                    ylabel('dFdx (FOM/m)');
                                else
                                    figure(30); hold all;
                                    title('Velocity: Real Space','FontWeight','bold','FontSize',12,'FontName','Arial');
                                    xlabel('x (m)'); ylabel('y (m)');
                                    
                                    imagesc(dFdx_mfu_bnd.');
                                end
                                
                                % Tranform to fourier space
                                dFdx_mfu = obj.fft(dFdx_mfu_bnd);
                                if obj.fixedAvgHeight
                                    dFdx_mfu_bnd = dFdx_mfu_bnd - dFdx_mfu(1,1,1);
                                    dFdx_mfu(1,1,1) = 0;
                                end
                                
                                if(obj.xBC==1) % symmetry about x-axis
                                    dFdx_mfu(:,:,2) = 0;
                                    dFdx_mfu(:,:,4) = 0;
                                end
                                
                                if(obj.yBC==1) % symmetry about y-axis
                                    dFdx_mfu(:,:,1) = 0;
                                    dFdx_mfu(:,:,3) = 0;
                                end
                                
                                if(obj.numY==1)
                                    figure(31); hold all;
                                    title('Velocity: Fourier Space','FontWeight','bold','FontSize',12,'FontName','Arial');
                                    xlabel('Cx');
                                    plot(dFdx_mfu(:));
                                    ylabel('dFdC (FOM/m)');
                                end
                                
                                obj.dFdx_bnd(m,f,u,:,:) = dFdx_mfu_bnd;
                                obj.dFdx(m,f,u,:,:,:) = dFdx_mfu;
                            end
                        end
                    end
                end
            end
                        
            ub = obj.maxMove*ones(obj.numCoeffx*obj.numCoeffy*4,1);
            lb = -ub;
            c = @(x)( obj.ConstraintFx(x) );
        end
        
        %% RETURN A STRING CONTAINING EQUATION FOR SURFACE
        % xVar and yVar are strings containing the names of the x and y
        % variables, respectively
        function equationString = surfaceEquation(obj, xVar, yVar, coeff, xLength, yLength)
            switch nargin
                case 1
                    xVar = 'x';
                    yVar = 'y';
                    coeff = obj.coeff * obj.dx * 1e6;
                    xLength = (obj.numX-1) * obj.dx * 1e6;
                    yLength = (obj.numY-1) * obj.dx * 1e6;
                case 3
                    coeff = obj.coeff * obj.dx * 1e6;
                    xLength = (obj.numX-1) * obj.dx * 1e6;
                    yLength = (obj.numY-1) * obj.dx * 1e6;
                case 4
                    xLength = (obj.numX-1) * obj.dx * 1e6;
                    yLength = (obj.numY-1) * obj.dx * 1e6;
            end
            
            equationString = ' ';
            for mh=1:obj.numCoeffx
                m = mh - 1;
                for nh=1:obj.numCoeffy
                    n = nh - 1;
                    if (n==0) && (m==0)
                        equationString = num2str(coeff(1,1,1));
                    else
                        if coeff(mh,nh,1) ~= 0
                            if m==0
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,1)), '*cos(', num2str(n),'*pi*', yVar,'/(', num2str(yLength),'/2))');
                            elseif n==0
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,1)), '*cos(', num2str(m),'*pi*', xVar, '/(', num2str(xLength), '/2))');
                            else
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,1)), '*cos(', num2str(m),'*pi*', xVar, '/(', num2str(xLength), '/2))*cos(', num2str(n),'*pi*', yVar,'/(', num2str(yLength),'/2))');
                            end
                        end
                        if coeff(mh,nh,2) ~= 0
                            if n==0
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,2)), '*sin(', num2str(m),'*pi*', xVar, '/(', num2str(xLength), '/2))');
                            elseif m ~= 0
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,2)), '*sin(', num2str(m),'*pi*', xVar, '/(', num2str(xLength), '/2))*cos(', num2str(n),'*pi*', yVar,'/(', num2str(yLength),'/2))');
                            end
                        end
                        if coeff(mh,nh,3) ~= 0
                            if m==0
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,3)), '*sin(', num2str(n),'*pi*', yVar,'/(', num2str(yLength),'/2))');
                            elseif n ~= 0
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,3)), '*cos(', num2str(m),'*pi*', xVar, '/(', num2str(xLength), '/2))*sin(', num2str(n),'*pi*', yVar,'/(', num2str(yLength),'/2))');
                            end
                        end
                        if coeff(mh,nh,4) ~= 0
                            if (m ~= 0) && (n ~=0)
                                equationString = strcat(equationString, '+', num2str(coeff(mh,nh,4)), '*sin(', num2str(m),'*pi*', xVar, '/(', num2str(xLength), '/2))*sin(', num2str(n),'*pi*', yVar,'/(', num2str(yLength),'/2))');
                            end
                        end
                    end
                end
            end
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
            dA = pi*obj.yr*deltaX(1) + pi*obj.xr*deltaX(2) ...
                + 4*obj.yr*deltaX(3) + 4*obj.yr*deltaX(4);
        end
        
        %% CONSTRUCT PADDING MATRICES
        function [padTop, padBot] = constructPadding(obj)
            padTop = zeros(obj.numX,obj.numY,obj.numZ);
            padBot = padTop;
            
            for i=1:obj.numX
                for j=1:obj.numY
                    padTop_zInd = obj.surface(i,j)+obj.eraseSize : obj.surface(i,j)+obj.eraseSize+obj.velPadding;
                    padBot_zInd = obj.surface(i,j)-obj.eraseSize-obj.velPadding : obj.surface(i,j)-obj.eraseSize;
                    
                    padTop_zInd(padTop_zInd>obj.numZ) = obj.numZ;
                    padTop_zInd(padTop_zInd<1) = 1;
                    padBot_zInd(padBot_zInd>obj.numZ) = obj.numZ;
                    padBot_zInd(padBot_zInd<1) = 1;
                    
                    padTop(i,j,round(padTop_zInd)) = 1;
                    padBot(i,j,round(padBot_zInd)) = 1;
                end
            end
        end
        
        %% INTERPOLATE VELOCITY
        function dFdx_bnd = interpolateVelocity(obj, dFdx_pad, pad)
            [ygrid, xgrid] = meshgrid(1:obj.numY, 1:obj.numX);
            
            % Use padded derivs to interpolate deriv along actual boundary
            
            if(obj.numY==1)
                pad_sqz = squeeze(pad);
                dFdx_pad_sqz = squeeze(dFdx_pad);
                [padIndicies] = find(pad_sqz(:,:,:));
                [padz,padx] = ind2sub(size(pad),find(pad_sqz(:,:)));
                velTri = TriScatteredInterp(padx, padz, dFdx_pad_sqz(padIndicies));
                dFdx_bnd = velTri(obj.surface, xgrid);
                dFdx_bnd = permute(dFdx_bnd, [1 3 2]);
            else
                [padIndicies] = find(pad(:,:,:));
                [pady, padx, padz] = ind2sub(size(pad),find(pad(:,:,:)));
                velTri = TriScatteredInterp(padx, pady, padz, dFdx_pad(padIndicies));
                dFdx_bnd = velTri(ygrid, xgrid, obj.surface);
            end
            dFdx_bnd(isnan(dFdx_bnd)) = 0;
        end
        
        %% FOURIER TRANSFORM VELOCITY
        function dFdx = fft(obj, dFdx_bnd)
            dFdx = zeros(obj.numCoeffx, obj.numCoeffy, 4);
            [y, x] = meshgrid(0:obj.numY-1, 0:obj.numX-1);

            for mh=1:obj.numCoeffx
                m = mh - 1;
                for nh=1:obj.numCoeffy
                    n = nh - 1;
                    
                    if( (mh==1) && (nh==1) )
                        factor = 1;
                    elseif( (mh==1) || (nh==1) )
                        factor = 2;
                    else
                        factor = 4;
                    end
                    
                    dFdx(mh, nh, 1) = factor * sum(sum( dFdx_bnd .* cos(m*pi*x/(obj.numX/2)) .* cos(n*pi*y/(obj.numY/2)) ));
                    dFdx(mh, nh, 2) = factor * sum(sum( dFdx_bnd .* sin(m*pi*x/(obj.numX/2)) .* cos(n*pi*y/(obj.numY/2)) ));
                    dFdx(mh, nh, 3) = factor * sum(sum( dFdx_bnd .* cos(m*pi*x/(obj.numX/2)) .* sin(n*pi*y/(obj.numY/2)) ));
                    dFdx(mh, nh, 4) = factor * sum(sum( dFdx_bnd .* sin(m*pi*x/(obj.numX/2)) .* sin(n*pi*y/(obj.numY/2)) ));
                end
            end
            
            dFdx = dFdx / (obj.numX*obj.numY) * (obj.numCoeffx * obj.numCoeffy)^2;
        end
        
        %% CALCULATE PAR/PERP FIELD COMPONENTS
        % findFieldNP(E, EA, pad, theta, epsPad)
        function [E_perp_x,E_perp_y,E_perp_z,E_par_x,E_par_y,E_par_z] = calcParPerpField(obj, E, pad)
            Ex = E{1};
            Ey = E{2};
            Ez = E{3};
            
            Px=0*Ex;
            E_par_x=Px;
            E_perp_x=Px;
            E_par_y=Px;
            E_perp_y=Px;
            E_par_z=Px;
            E_perp_z=Px;
            
            normal_vec=zeros(obj.numX,obj.numY,3);
            
            for i=1:obj.numX
                for j=1:obj.numY
                    xDev = -obj.xDerivativeSurface(i,j);
                    yDev = -obj.yDerivativeSurface(i,j);
                    v1=[1 0 xDev];
                    v2=[0 1 yDev];
                    n=cross(v1,v2)/norm(cross(v1,v2));
                    normal_vec(i,j,:)=n;
                    E_vec=[squeeze(Ex(i,j,:)),squeeze(Ey(i,j,:)),squeeze(Ez(i,j,:))];
                    P_vec=zeros(size(E_vec));
                    E_perp=P_vec;
                    E_par=P_vec;
                    for k=1:length(Ex(i,j,:))
                        E_perp(k,:)=(sum(E_vec(k,:).*n))*n;
                        E_par(k,:)=E_vec(k,:)-E_perp(k,:);
                    end
                    E_par_x(i,j,:)=E_par(:,1);
                    E_perp_x(i,j,:)=E_perp(:,1);
                    E_par_y(i,j,:)=E_par(:,2);
                    E_perp_y(i,j,:)=E_perp(:,2);
                    E_par_z(i,j,:)=E_par(:,3);
                    E_perp_z(i,j,:)=E_perp(:,3);
                    %E_perp_tot(i,j,:)=E_perp(:,1)+E_perp(:,2)+E_perp(:,3);
                end
            end
            %E_par_tot=squeeze((E_par_x(:,10,:).*conj(E_par_x(:,10,:))+E_par_y(:,10,:).*conj(E_par_y(:,10,:))+E_par_z(:,10,:).*conj(E_par_z(:,10,:))));
            %E_per_tot=squeeze((E_perp_x(:,10,:).*conj(E_perp_x(:,10,:))+E_perp_y(:,10,:).*conj(E_perp_y(:,10,:))+E_perp_z(:,10,:).*conj(E_perp_z(:,10,:))));
            
            E_par_x = E_par_x .* pad;
            E_par_y = E_par_y .* pad;
            E_par_z = E_par_z .* pad;
            E_perp_x = E_perp_x .* pad;
            E_perp_y = E_perp_y .* pad;
            E_perp_z = E_perp_z .* pad;
        end
        
        %% X-DERIVATIVE SURFACE
        function value = xDerivativeSurface(obj, x, y)
            value = 0;
            x=x-1;
            y=y-1;
            
            for mh=1:obj.numCoeffx
                m=mh-1;
                for nh=1:obj.numCoeffy
                    n=nh-1;
                    value = value -  obj.coeff(mh,nh,1)*(m*pi/(obj.numX/2))*sin(m*pi*x/(obj.numX/2))*cos(n*pi*y/(obj.numY/2));
                    value = value +  obj.coeff(mh,nh,2)*(m*pi/(obj.numX/2))*cos(m*pi*x/(obj.numX/2))*cos(n*pi*y/(obj.numY/2));
                    value = value -  obj.coeff(mh,nh,3)*(m*pi/(obj.numX/2))*sin(m*pi*x/(obj.numX/2))*sin(n*pi*y/(obj.numY/2));
                    value = value +  obj.coeff(mh,nh,4)*(m*pi/(obj.numX/2))*cos(m*pi*x/(obj.numX/2))*sin(n*pi*y/(obj.numY/2));
                end
            end
            value = -value;
        end
        
        %% Y-DERIVATIVE SURFACE
        function[value] = yDerivativeSurface(obj, x, y)
            value = 0;
            x=x-1;
            y=y-1;
            
            for mh=1:obj.numCoeffx
                m=mh-1;
                for nh=1:obj.numCoeffy
                    n=nh-1;
                    value = value - obj.coeff(mh,nh,1)*cos(m*pi*x/(obj.numX/2))*(n*pi/(obj.numY/2))*sin(n*pi*y/(obj.numY/2));
                    value = value - obj.coeff(mh,nh,2)*sin(m*pi*x/(obj.numX/2))*(n*pi/(obj.numY/2))*sin(n*pi*y/(obj.numY/2));
                    value = value + obj.coeff(mh,nh,3)*cos(m*pi*x/(obj.numX/2))*(n*pi/(obj.numY/2))*cos(n*pi*y/(obj.numY/2));
                    value = value + obj.coeff(mh,nh,4)*sin(m*pi*x/(obj.numX/2))*(n*pi/(obj.numY/2))*cos(n*pi*y/(obj.numY/2));
                end
            end
            value = -value;
        end

        %% returnData()
        % return nk data for import to Lumerical
        function data = returnData(obj)
            if( ischar(obj.epsBelow) )
                geo_n = 0;
                geo_mat = obj.epsBelow;
            else
                geo_n = sqrt(obj.epsBelow);
                geo_mat = '';
            end
            
            if( ischar(obj.epsAbove) )
                geo_nClad = 0;
                geo_matClad = obj.epsAbove;
            else
                geo_nClad = sqrt(obj.epsAbove);
                geo_matClad = '';
            end
            
            geo_x = obj.x0;
            geo_xspan = (obj.numX-1) * obj.dx;
            if(obj.numY==1)
                geo_y = obj.z0;
                geo_yspan = (obj.numZ-1) * obj.dx;
                geo_z = -2*obj.dx;
                geo_zspan = 4*obj.dx;
            else
                geo_y = obj.y0;
                geo_yspan = (obj.numY-1) * obj.dx;
                geo_z = obj.z0;
                geo_zspan = (obj.numZ-1) * obj.dx;
            end
            
            geo_surface = obj.surfaceEquation('u','v');
            data = struct('geo_nClad',geo_nClad,'geo_matClad',geo_matClad,'geo_mat',geo_mat,'geo_n',geo_n,'geo_x',geo_x,'geo_xspan',geo_xspan,'geo_y',geo_y,'geo_yspan',geo_yspan,'geo_z',geo_z,'geo_zspan',geo_zspan,'geo_surface',geo_surface);
        end
        
        %% getDrawing()
        % return permitivity data for Matlab display to user
        function drawing = getDrawing(obj, x, y, z)
            drawing = obj.surface * obj.dx;
        end
        
    end
    
end
