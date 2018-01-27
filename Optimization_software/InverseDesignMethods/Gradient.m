classdef Gradient
    
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
        
        %Data = cell(numMon, numFreq, numUD) [numY, numX]
        E;
        EA;
        epsBgnd;
        epsVec; % vector of actual eps values
        epsOutVec; % vector of actual eps values
        
        % Derivative matrix
        % size = [numMon, numFreq, numUD, numY, numX]
        dFdxBnd;
        
        dFdxSpace;
        newShape;
        newShapeGrid;
        dFGrid;
        
        % Boundary Field Approximation
        eraseSize; % size of ignored field in velocity calculation
        velPadding; % size of region in velocity calculation
        
        % These should be set by the Optimizer object
        % Note that numMon is the number of monitors SIMULATIONS, and that
        % numFreq is the total number of frequencies
        numMon;
        numFreq;
        numUD;
        
        % Test flag (displaying output if 1)
        testFlag;
    end
    
    methods
        function obj = Gradient(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eraseSize, velPadding)
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
            if(obj.numY<0)
                error('y_span must be > 0');
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
            
            obj.eraseSize = floor(1e-4*round(1e4*eraseSize / dx));
            obj.velPadding = floor(1e-4*round(1e4*velPadding / dx));
            if(obj.eraseSize<0)
                error('eraseSize must be > 0');
            end
            if(obj.velPadding<1)
                error('velPadding must be > 1');
            end
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
        
        %% INIT DATA
        % initData()
        function[obj] = initData(obj)
            obj.E = cell(1,obj.numFreq, obj.numUD);
            obj.EA = cell(obj.numMon,obj.numFreq, obj.numUD);
            obj.epsBgnd = cell(obj.numFreq, 1);
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
        
        %% INTERPOLATE FIELDS ONTO GEOMETRY MESH
        % DO NOT AVERAGE FIELDS IN Z, give entire array to the shape
        % Note that this is also used for eps = {epsx,epsy,epsz}
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
                
        %% CALCULATE DERIVATIVES
        % Note: if numMonSims > 1, minDim = 1
        % dFdx has size [numMon, numFreq, numUD, numY, numX]
        function[obj] = calcDerivs(obj, skipInd, freq, epsGrid)
            obj.dFdxBnd = zeros(obj.numMon, obj.numFreq, obj.numUD, obj.numY, obj.numX);
            
            % Find derivatives at all boundaries
            obj.dFdxBnd = obj.calcBoundaryDerivs(epsGrid);
            
            % Find derivatives for new shape creation
            [obj.dFdxSpace] = obj.calcShapeDerivs(epsGrid, skipInd);
        end
        
        %% CALCULATE NEW SHAPE DERIVATIVES
        % newShape = [x, y, radius, dFdA]
        function dFdxSpace = calcShapeDerivs(obj, epsGrid, skipInd)
            erasePad = obj.getErasePad(epsGrid);
            
            [padIndices] = find(~erasePad);
            [pady, padx] = find(~erasePad(:,:,1));
            
            [xgrid, ygrid] = meshgrid(1:obj.numX, 1:obj.numY);
            
            dFdxSpace = zeros(obj.numY, obj.numX);
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
                            dFdxSpace(:,:,c) = dFdxShape_mfu;
                            
                            if(obj.numY>1)
                                figure(30+f); imagesc(dFdxShape_mfu); colormap(bluewhitered); axis equal; title('dFdx: Sparse Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                            else
                                figure(30+f); plot(dFdxShape_mfu); title('dFdx: Sparse Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                            end
                        end
                    end
                end
            end
            
        end
        
        %% CALCULATE BOUNDARY DERIVATIVES
        function dFdx = calcBoundaryDerivs(obj, epsGrid)
            padOut = obj.getPadOut(epsGrid);
            padIn = obj.getPadIn(epsGrid);
            pad = padOut + padIn;
            theta = obj.getThetaOut(epsGrid) + obj.getThetaIn(epsGrid);
            
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
                                
                                dFdx_mfu = obj.calcV(epsGrid, Ep, EAp, Dn, DAn, Ez, EAz, pad, obj.epsVec(f), obj.epsOutVec(f));
                                dFdx(m,f,u,:,:) = dFdx_mfu;
                                
                                vLim = max(max( abs(dFdx_mfu) ));
                                if(obj.numY>1)
                                    figure(40+f); imagesc(dFdx_mfu); axis equal; caxis([-vLim vLim]); colormap(bluewhitered); title('dFdx: Boundary Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                                else
                                    figure(40+f); plot(dFdx_mfu); title('dFdx: Boundary Perturbation','FontWeight','bold','FontSize',12,'FontName','Arial');
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
        %% CALCULATE BOUNDARY DERIVATIVES - HELPER
        % Averages dFdx_mfu about the extruded dimension (z-axis)
        function dFdx_mfu = calcV(obj, epsGrid, Ep, EAp, Dn, DAn, Ez, EAz, pad, eps, epsOut)
            dFdxPad_p = real( (eps - epsOut) .* Ep .* EAp ); % par component
            if(0*real(obj.epsVec(1))<0)
                dFdxPad_pz = real( (eps - epsOut)/exp(0*1i*angle(eps)+angle(epsOut)) .* Ez .* EAz ); % z component (metal thin-film)
            else
                dFdxPad_pz = real( (eps - epsOut) .* Ez .* EAz ); % par-z component
            end
            dFdxPad_n = real( Dn .* DAn .* (1./epsOut - 1./eps) ); % perp component
            %absorption = 1/7.57e-10 * 0.5* (2*pi*3e8/836e-9) * 8.85e-12 * imag(obj.epsVec(1)) * (abs(Ep).^2 + abs(Ez).^2 + abs((Dn/obj.epsVec(1))).^2);
            
            dFdxPad = dFdxPad_p + dFdxPad_pz + dFdxPad_n;% + .05*absorption;
            
            if(0 && obj.testFlag)
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
                
                clim = max(max(max(abs(dFdxPad_p(:,:,1))))) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0));
                figure(50); imagesc(dFdxPad_p(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); caxis([-clim clim]); colormap(bluewhitered); colorbar;
                clim = max(max(max(abs(dFdxPad_pz(:,:,1))))) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0));
                figure(51); imagesc(dFdxPad_pz(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); caxis([-clim clim]); colormap(bluewhitered); colorbar;
                clim = max(max(max(abs(dFdxPad_n(:,:,1))))) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0));
                figure(52); imagesc(dFdxPad_n(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); caxis([-clim clim]); colormap(bluewhitered); colorbar;
                %clim = max(max(max(abs(absorption(:,:,1))))) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0));
                %figure(53); imagesc(absorption(:,:,1)* obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); caxis([-clim clim]); colormap(bluewhitered); colorbar;
                clim = max(max(max(abs(dFdxPad(:,:,1))))) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0));
                figure(54); imagesc(dFdxPad(:,:,1) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0))); caxis([-clim clim]); colormap(bluewhitered); colorbar;
            
            end
            
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
                dFdx_mfu(:,:,i) = dFdx_z .* obj.getBoundary(epsGrid) * obj.dx^(1 + 1*(obj.numY>1) + 1*(obj.numZ~=0)) * numZ;
            end
            dFdx_mfu = mean(dFdx_mfu,3);
            
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
            if(0*real(obj.epsVec(1))<0)
                Ez = Ez .* exp(1i*angle(epsPad)) .* (epsPad~=0);
            else
                Ez = Ez .* pad;
            end
            
            DAn = (EAx.*cos(theta) + EAy.*sin(theta)) .* epsPad;
            EAp = (EAx.*sin(theta) - EAy.*cos(theta)) .* pad;
            if(0*real(obj.epsVec(1))<0)
                EAz = EAz .* exp(1i*angle(epsPad)) .* (epsPad~=0);
            else
                EAz = EAz .* pad;
            end
            
            
        end

        %% getThetaOut()
        function [thetaOut] = getThetaOut(obj, epsGrid)
            Nx = obj.numX;
            Ny = obj.numY;
            
            thetaOut = zeros(obj.numY, obj.numX);
            gridNew = epsGrid;
            thetaNew = calcTheta(gridNew);
            
            num = obj.eraseSize + obj.velPadding;
            
            for i=1:num
                thetaPad = zeros(Ny+2, Nx+2);
                boundary = getBoundaryOut(gridNew) + getBoundaryIn(gridNew);
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
                
                gridOut = getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                
                if(i > obj.eraseSize)
                    thetaOut = thetaOut + thetaNew.*(1*gridOut);
                end
            end
        end
        
        %% getThetaIn()
        function [thetaIn] = getThetaIn(obj, epsGrid)
            Nx = obj.numX;
            Ny = obj.numY;
            
            thetaIn = zeros(obj.numY, obj.numX);
            gridNew = epsGrid;
            thetaNew = calcTheta(gridNew);
            
            num = obj.eraseSize + obj.velPadding;
            
            for i=1:num
                thetaPad = zeros(Ny+2, Nx+2);
                boundary = getBoundaryOut(gridNew) + getBoundaryIn(gridNew);
                thetaPad(2:Ny+1,2:Nx+1) = boundary.*exp(1i*thetaNew); % exponetial, so angles can be averaged properly
                
                thetaPad2=thetaPad(2:Ny+1,2:Nx+1)+thetaPad(1:Ny,1:Nx)...
                    +thetaPad(1:Ny,3:Nx+2)+thetaPad(3:Ny+2,3:Nx+2)...
                    +thetaPad(3:Ny+2,1:Nx)+thetaPad(1:Ny,2:Nx+1)...
                    +thetaPad(2:Ny+1,1:Nx)+thetaPad(3:Ny+2,2:Nx+1)...
                    +thetaPad(2:Ny+1,3:Nx+2);
                thetaPad2 = 1e-4*round(1e4*thetaPad2); % to remove miniscule numbers which should be zero
                thetaNew = angle(thetaPad2);
                
                gridIn = getBoundaryIn(gridNew);
                gridNew = gridNew & ~gridIn;
                
                if(i > obj.eraseSize)
                    thetaIn = thetaIn + thetaNew.*(1*gridIn);
                end
            end
        end
        
        %% getErasePad()
        function erasePad = getErasePad(obj, epsGrid)
            erasePad = zeros(obj.numY, obj.numX);
            gridNew = epsGrid;
            for i=1:obj.eraseSize
                gridOut = getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                erasePad = erasePad | gridOut;
            end
            gridNew = epsGrid;
            for i=1:obj.eraseSize
                gridIn = getBoundaryIn(gridNew);
                gridNew = ~gridIn & gridNew;
                erasePad = erasePad | gridIn;
            end
        end
        
        %% getPadOut()
        function padOut = getPadOut(obj, epsGrid)
            num = obj.eraseSize + obj.velPadding;
            padOut = zeros(obj.numY, obj.numX);
            gridNew = epsGrid;
            for i=1:num
                gridOut = getBoundaryOut(gridNew);
                gridNew = gridOut | gridNew;
                if(i > obj.eraseSize)
                    padOut = padOut | gridOut;
                end
            end
        end
        
        %% getPadIn()
        function padIn = getPadIn(obj, epsGrid)
            num = obj.eraseSize + obj.velPadding;
            padIn = zeros(obj.numY, obj.numX);
            gridNew = epsGrid;
            for i=1:num
                gridIn = getBoundaryIn(gridNew);
                gridNew = ~gridIn & gridNew;
                if(i > obj.eraseSize)
                    padIn = padIn | gridIn;
                end
            end
        end
        
        %% getBoundary()
        function epsGridBound = getBoundary(obj, epsGrid)
            epsGridBound = getBoundaryIn(epsGrid) +getBoundaryOut(epsGrid);
        end
        
    end
    
end
