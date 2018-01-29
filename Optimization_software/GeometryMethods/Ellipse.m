classdef Ellipse
    % Ellipse class
    % TODO: set normalization for dF (thickness?)
    %       rename epsBgnd
    %       Note: rename dFdx0... derivs to not include area
    %       Note: maybe remove lineStep
    %       NEED TO FIX TREATMENT OF BACKGROUND EPS
    
    %% PROPERTIES
    properties
        epsVec; % vector of actual epsClad values
        eps_; % single eps_ value or 'material name'
        dynamic; % 0=static
        meshOrder; % 1=highest priority
        
        %Shape Definition on internal mesh
        x0;
        y0;
        xr;
        yr;
        rules;
        
        %Boundary Definition on internal mesh
        xPad; 
        yPad;
        VPad;
        xLine;
        yLine;
        nxLine;
        nyLine;
        VLine;
        VErr;
        
        % Normals and tangents on external mesh
        nx;
        ny;
        tx;
        ty;
        
        %New Uncommitted Paramaters
%         x0_new;
%         y0_new;
%         xr_new;
%         yr_new;
        dF;
%         dx0;
%         dy0;
%         dxr;
%         dyr;
        dFdx0;
        dFdy0;
        dFdxr;
        dFdyr;
        dA;
        
        %Flags
        newShape; %1 = new shape created by geometry.createShapes()
        testFlag;
        % add bounds on all parameters (minRadius, maxRadius, etc.)
        % add actual boundary limitation

    end
    
    methods
        %% CONSTRUCTOR
        % params = {x0, y0, xr, yr}
        function obj = Ellipse(eps_, dynamic, meshOrder, rules, params, testFlag)
            obj.eps_ = eps_;
            obj.dynamic = dynamic;
            obj.meshOrder = meshOrder;
            obj.x0 = params(1);
            obj.y0 = params(2);
            obj.xr = params(3);
            obj.yr = params(4);
            % wrappers must check that rules are correct
            obj.rules = rules;
            obj.newShape = 1;
			obj.testFlag = testFlag;
			
            if(obj.testFlag)
                if(obj.dynamic)
                    fprintf('%%%% ADDING DYNAMIC ELLIPSE %%%% \n');
                else
                    fprintf('%%%% ADDING STATIC ELLIPSE %%%% \n');
                end
                fprintf('  (x0,y0) = (%g,%g) \n',obj.x0,obj.y0);
                fprintf('  (xr,yr) = (%g,%g) \n',obj.xr,obj.yr);
                if(ischar(obj.eps_))
                    fprintf('  eps_ = %s \n',obj.eps_);
                else
                    fprintf('  eps_ = %g \n',obj.eps_);
                end
            end
        end
        
        %% ACCEPT NEW PARAMETER VALUES
        function obj = accept(obj, deltaX)
            obj.xr = obj.xr + deltaX(1);
            obj.yr = obj.yr + deltaX(2);
            obj.x0 = obj.x0 + deltaX(3);
            obj.y0 = obj.y0 + deltaX(4);
        end
        
        %% CALCULATE DERIVATIVES OF ELLIPSE GIVEN FIELDS
        % Take Arrays of field in Z, average velocity at the end
        % dFdx is of size [numMon x numFreq x numUD x 4]
        function [obj,dFdx] = calcDerivs(obj, E, EA, epsBgnd, eraseSize, velPadding, lineStep, dx)
            
            % Find boundary points, padding points for given shape
            obj = obj.padding(eraseSize, velPadding, lineStep);
            
            % Calculate V on padding, average back to line
            obj = obj.calcV(E, EA, epsBgnd);
%             if(obj.VErr > 1) % set this to some tolerance!!
%                 fprintf('Error too high: %g \n',obj.VErr);
%             end
            
            [numMon,numFreq,numUD,~] = size(obj.VLine);

            % Calculate derivative w/ respect to given parameters
            % find theta around ellipse
            numT = length(obj.xLine);
            dT = 2*pi / numT;
            theta = linspace(0,2*pi-dT,numT);
            cosT = cos(theta);
            sinT = sin(theta);
            
            % Actual derivatives with respect to parameters
            % Factor of dx so that dFdxr in mesh-units
            dFdx = zeros(numMon,numFreq,numUD,4);
            for m=1:numMon
                for f=1:numFreq
                    for u=1:numUD
                        V = reshape(obj.VLine(m,f,u,:),1,[]);
                        dFdx(m,f,u,1) = obj.yr * dx * trapz(theta, V .* cosT.^2) * lineStep * dx; %dFdxr
                        dFdx(m,f,u,2) = obj.xr * dx * trapz(theta, V .* sinT.^2) * lineStep * dx; %dFdyr
                        dFdx(m,f,u,3) = obj.yr * dx * trapz(theta, V .* cosT) * lineStep * dx; %dFdx0
                        dFdx(m,f,u,4) = obj.xr * dx * trapz(theta, V .* sinT) * lineStep * dx; %dFdx0
                    end
                end
            end

            % Enforce any rules regarding parameters / points needing to
            % stay fixed
            obj = enforceRules(obj);
            
%             dFdx = [obj.dFdxr, obj.dFdyr, obj.dFdx0, obj.dFdy0];
        end
        
        %% RETURN CONSTRAINTS
        function[dA] = deltaArea(obj,deltaX)
            deltaX = abs(deltaX);
            dA = pi*obj.yr*deltaX(1) + pi*obj.xr*deltaX(2) ...
                + 4*obj.yr*deltaX(3) + 4*obj.xr*deltaX(4);
        end
        
        %% CALCULATE NORMAL AND TANGENT VECTORS ON EXTERNAL MESH
        function[obj] = calcNormals(obj, xGrid, yGrid)
            % Using level set function phi = (x-x0)^2/xr^2 + (y-y0)^2/yr^2 - 1
            dPhiDx = 2*(xGrid-obj.x0)/obj.xr^2;
            dPhiDy = 2*(yGrid-obj.y0)/obj.yr^2;
            obj.nx = dPhiDx ./ (sqrt( dPhiDx.^2 + dPhiDy.^2 )+1e-6); % ensures no division by 0
            obj.ny = dPhiDy ./ (sqrt( dPhiDx.^2 + dPhiDy.^2 )+1e-6);
            obj.tx = -obj.ny; % Tangent vector in plane
            obj.ty = obj.nx;
        end
        
        %% CALCULATE V ON EXTERNAL MESH
        % E is assumed to have dimensions {1 x numFreq x numUser}{1x3}
        % EA is assumed to have dimensions {numMon x numFreq x numUser}{1x3}
        % For now, it is assumed that eps_ = avg(epsx, epsy, epsz)
        % VLine is a matrix of size [numMon x numFreq x numUser x padSize]
        function[obj, xGrid, yGrid, VCell, epsCell] = calcV(obj, E, EA, epsBgnd)
            
            szPad = size(obj.xPad);
            [xGrid,yGrid] = meshgrid(1:size(E{1,1}{1},2), 1:size(E{1,1}{1},1));
            obj = obj.calcNormals(xGrid, yGrid);
            
            [numMon, numFreq, numUser] = size(EA);
            VCell = cell(numFreq,numUser);
            epsCell = cell(numFreq, 1);
            obj.VLine = zeros(numMon,numFreq,numUser,szPad(2));
            
            numZ = size(epsBgnd{1}{1},3);
            tx = repmat(obj.tx,[1 1 numZ]);
            ty = repmat(obj.ty,[1 1 numZ]);
            nx = repmat(obj.nx,[1 1 numZ]);
            ny = repmat(obj.ny,[1 1 numZ]);
            
            for m = 1:numMon
                for f = 1:numFreq
                    epsi = epsBgnd{f}{1} + epsBgnd{f}{2} + epsBgnd{f}{3};
                    numEps = 1*(max(max(max(abs(epsBgnd{f}{1})))) > 0) + 1*( max(max(max(abs(epsBgnd{f}{2})))) > 0) + 1*(max(max(max(abs(epsBgnd{f}{3})))) > 0);
                    epsi = epsi / numEps;
                    for zInd = 1:size(epsi,3)
                        FepsR = TriScatteredInterp(reshape(xGrid,[],1),reshape(yGrid,[],1),reshape(real(epsi(:,:,zInd)),[],1));
                        FepsI = TriScatteredInterp(reshape(xGrid,[],1),reshape(yGrid,[],1),reshape(imag(epsi(:,:,zInd)),[],1));
                        epsPad = FepsR(obj.xPad, obj.yPad) + 1i*FepsI(obj.xPad, obj.yPad);
                        epsInside = repmat( epsPad(szPad(1)/2,:), szPad(1), 1);
                        epsOutside = repmat( epsPad(szPad(1)/2+1,:), szPad(1), 1);
                        deltaEps(:,:,zInd) = epsInside - epsOutside;
                        deltaEpsInv(:,:,zInd) = 1./epsOutside - 1./epsInside;
                    end
                    deltaEps(isnan(deltaEps)) = 0;
                    deltaEpsInv(isnan(deltaEpsInv)) = 0;
                    
                    for u = 1:numUser % numUser
                        if( ~isempty(EA{m,f,u}) )
                            Epar1 = E{1,f,u}{1} .* tx + E{1,f,u}{2} .* ty;
                            Epar2 = E{1,f,u}{3};
                            Dperp = epsBgnd{f}{1} .* E{1,f,u}{1} .* nx + epsBgnd{f}{2} .* E{1,f,u}{2} .* ny;
                            EApar1 = EA{m,f,u}{1} .* tx + EA{m,f,u}{2} .* ty;
                            EApar2 = EA{m,f,u}{3};
                            DAperp = epsBgnd{f}{1} .* EA{m,f,u}{1} .* nx + epsBgnd{f}{2} .*EA{m,f,u}{2} .* ny;
                            
                            Vpar = Epar1 .* EApar1 + Epar2 .* EApar2;
                            Vperp = Dperp .* DAperp;
                            
                            for zInd = 1:size(epsi,3)
                                FparR = TriScatteredInterp(reshape(xGrid,[],1),reshape(yGrid,[],1),reshape(real(Vpar(:,:,zInd)),[],1));
                                FparI = TriScatteredInterp(reshape(xGrid,[],1),reshape(yGrid,[],1),reshape(imag(Vpar(:,:,zInd)),[],1));
                                FperpR = TriScatteredInterp(reshape(xGrid,[],1),reshape(yGrid,[],1),reshape(real(Vperp(:,:,zInd)),[],1));
                                FperpI = TriScatteredInterp(reshape(xGrid,[],1),reshape(yGrid,[],1),reshape(imag(Vperp(:,:,zInd)),[],1));
                                VparPad(:,:,zInd) = FparR(obj.xPad, obj.yPad) + 1i*FparI(obj.xPad, obj.yPad);
                                VperpPad(:,:,zInd) = FperpR(obj.xPad, obj.yPad) + 1i*FperpI(obj.xPad, obj.yPad);
                            end
                            
                            VparPad(isnan(VparPad)) = 0;
                            VperpPad(isnan(VperpPad)) = 0;
                            
                            Vmfu = 2*real( deltaEps.*VparPad + deltaEpsInv.*VperpPad );
                            Vmfu = mean( mean(Vmfu,1), 3); % 1 = mean over padded region, 3 = mean over z
                            VCell(m,f,u) = {Vmfu};
                            epsCell(m,f,u) = {epsInside};
                            obj.VLine(m,f,u,:) = reshape(Vmfu,1,1,1,[]);
                        end
                    end
                end
            end
%             obj.VErr = max( std(obj.VPad, 1, 1) ./ mean((mean(abs(obj.VPad)))) );
        end
        
        %% PAD THE ELLIPSE TO DELETE VELOCITY DATA AND INTERPOLATE
        function obj = padding(obj, eraseSize, velPadding, lineStep)
            % Units for velPadding, eraseSize, and lineStep in internal mesh cells
            % Perimeter approximation from Ramanujan
            perim = pi*( 3*(obj.xr+obj.yr) - sqrt( (3*obj.xr+obj.yr)*(obj.xr+3*obj.yr) ) );
            numPts = round( perim / lineStep );
            
            % xLine, yLine have dimensions [1,numPts]
            dTheta = 2*pi / numPts;
            theta = linspace(0,2*pi-dTheta,numPts);
            cosT = cos(theta); 
            sinT = sin(theta);
            obj.xLine = obj.x0 + obj.xr*cosT;
            obj.yLine = obj.y0 + obj.yr*sinT;
            obj.nxLine = cosT .* ( cosT.^2 + (obj.xr/obj.yr)^2 * sinT.^2 ).^(-1/2);
            obj.nyLine = sinT .* ( (obj.yr/obj.xr)^2 * cosT.^2 + sinT.^2 ).^(-1/2);
            
            % delta = vector containing distances from the boundary
            %   e.g. if eraseSize = 1 & velPadding = 2
            %   delta = [-2,-1,1,2]
            delta = [(-velPadding-eraseSize+1:-eraseSize),(eraseSize:eraseSize+velPadding-1)].';
            
            % xPad, yPad have dimensions [2*velPadding, numPoints]
            % xPad, yPad in order from farthest inside (row 1) to farthest outside (row 2*velPadding)
            obj.xPad = repmat(obj.xLine,2*velPadding,1) + delta * cosT;
            obj.yPad = repmat(obj.yLine,2*velPadding,1) + delta * sinT;

        end
        
        %% ENFORCE RULES
        function obj = enforceRules(obj)
            % ? Add rules for aspect ratio, fix one radius
            % 7 rules
            % fix x0, y0 
            % fix to a circle
            % fix four possible points
            % fix xr, yr
            % fix xr/yr (aspect ratio)

            % This should be given as a constraint to optimization routine
                
%             % rule 1: fix x0 to a given point
%             if( obj.rules(1) )
%                 obj.dx0 = 0;
%             end
%             
%             % rule 2: fix y0 to a given location
%             if( obj.rules(2) )
%                 obj.dy0 = 0;
%             end
%             
%             % rule 3: fix ellipse to a circle
%             if( obj.rules(3) )
%                 obj.dxr = 0.5*(Dxr+Dyr);
%                 obj.dyr = obj.dxr;
%             end
%             
%             % rule 4: fix x1 (make same as rectangle order)
%             %   where x1 is x0 - xr
%             if( obj.rules(4) )
%                 obj.dx0 = 0.5*(Dx0+Dxr);
%                 obj.dxr = obj.dx0;
%             end
%             
%             % rule 5: fix y1
%             %   where y1 is y0 - yr
%             if( obj.rules(5) )
%                 obj.dy0 = 0.5*(Dy0+Dyr);
%                 obj.dyr = obj.dy0;
%             end
%             
%             % rule 6: fix x2
%             %   where x2 is x0 + xr
%             if( obj.rules(6) )
%                 obj.dx0 = 0.5*(Dx0-Dxr);
%                 obj.dxr = -obj.dx0;
%             end
%             
%             % rule 7: fix y2
%             %   where y2 is y0 + yr
%             if( obj.rules(7) )
%                 obj.dy0 = 0.5*(Dy0-Dyr);
%                 obj.dyr = -obj.dy0;
%             end
%             
%             % rule 8: fix xr
%             if( obj.rules(8) )
%                 obj.dxr = 0;
%             end
%             
%             % rule 9: fix yr
%             if( obj.rules(9) )
%                 obj.dyr = 0;
%             end
%             
%             % rule 10: fix aspect ratio (xr/yr)
%             if( obj.rules(10) )
%                 obj.dxr = 0.5*(Dx0+obj.xr*Dyr/obj.yr);
%                 obj.dyr = obj.yr*Dxr/obj.xr;
%             end
        end
        
        %% CHECK IF ELLIPSE TOO SMALL (IN WHICH CASE MAYBE DELETE)
        function flag = tooSmall(obj,smallSize)
            flag = 0;
            if (obj.xr < smallSize  || obj.yr < smallSize)
                flag = 1;
            end
        end
        
        %% GET DRAWING OF CURRENT ELLIPSE
        % padding in units of internal mesh
        function [mask] = getDrawing(obj, xLengthGeo, yLengthGeo, padding)
            [xGrid, yGrid] = meshgrid(1:xLengthGeo, 1:yLengthGeo);
            mask = 1.*( ((xGrid-obj.x0).^2/(obj.xr+padding)^2 + (yGrid-obj.y0).^2/(obj.yr+padding)^2) <=1 );
        end
        
        %% GET DRAWING OF NEW ELLIPSE
        function [mask] = getNewDrawing(obj, xLengthGeo, yLengthGeo, padding, deltaX)
            xrnew = obj.xr + deltaX(1);
            yrnew = obj.yr + deltaX(2);
            x0new = obj.x0 + deltaX(3);
            y0new = obj.y0 + deltaX(4);
            [xGrid, yGrid] = meshgrid(1:xLengthGeo, 1:yLengthGeo);
            mask = 1.*( ((xGrid-x0new).^2/(xrnew+padding)^2 + (yGrid-y0new).^2/(yrnew+padding)^2) <=1 );
        end
    end
    
end
