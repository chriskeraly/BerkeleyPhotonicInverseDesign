classdef Rectangle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
	% Fix treatment of parallel components!!
    
    properties
        epsVec; % vector of actual epsClad values
        eps_; % single eps_ value or 'material name'
        dynamic; % 0=static
        meshOrder; % 1=highest priority
        
        %Shape Definition on internal mesh
        x0;
        y0;
        xLength;
        yLength;
        rules;
        
        %Boundary Definition on internal mesh
        xLine;
        yLine;
        nxLine;
        nyLine;
        VLine;
        left; top; right; bottom; % indices of the rectangle's sides
        
        %New Uncommitted Paramaters
        x0_new;
        y0_new;
        xLength_new;
        yLength_new;
        dF;
        dx0;
        dy0;
        dxLength;
        dyLength;
        dxr;
        dyr;
        dFdx0;
        dFdy0;
        dFdxLength;
        dFdyLength;
        dA;
        
        % add bounds on all parameters (minRadius, maxRadius, etc.)
        % don't let padding past maxIndices
        
        %Flags
        newShape; %1 = new shape created by geometry.createShapes()
		testFlag;
		
    end
    
    methods
        function obj = Rectangle(eps_, dynamic, meshOrder, rules, params, testFlag)
            obj.eps_ = eps_;
            obj.dynamic = dynamic;
            obj.meshOrder = meshOrder;
            obj.x0 = params(1);
            obj.y0 = params(2);
            obj.xLength = params(3);
            obj.yLength = params(4);
            % wrappers must check that rules are correct
            obj.rules = rules;
            obj.newShape=1;
			obj.testFlag = testFlag;
			
			if(obj.testFlag)
				if(obj.dynamic)
					fprintf('%%%% ADDING DYNAMIC RECTANGLE %%%% \n');
					fprintf('  (x0,y0) = (%g,%g) \n',obj.x0,obj.y0);
					fprintf('  (numX,numY) = (%g,%g) \n',obj.xLength,obj.yLength);
					fprintf('  eps_ = %g \n',obj.eps_);
				else
					fprintf('%%%% ADDING STATIC RECTANGLE %%%% \n');
					fprintf('  (x0,y0) = (%g,%g) \n',obj.x0,obj.y0);
					fprintf('  (numX,numY) = (%g,%g) \n',obj.xLength,obj.yLength);
					fprintf('  eps_ = %g \n',obj.eps_);
				end
			end
        end
        
        %% ACCEPT NEW PARAMETER VALUES
        function obj = accept(obj, deltaX)
            deltaX = obj.enforceRules(deltaX);
            obj.x0 = obj.x0 + deltaX(1);
            obj.y0 = obj.y0 + deltaX(2);
            obj.xLength = obj.xLength + deltaX(3);
            obj.yLength = obj.yLength + deltaX(4);
        end
        
        %% CALCULATE DERIVATIVES OF RECTANGLE GIVEN FIELDS
        % Take Arrays of field in Z, average velocity at the end
        % dFdx is of size [numMon x numFreq x numUD x 4]
        function [obj,dFdx] = calcDerivs(obj, E, EA, epsBgnd, eraseSize, velPadding, lineStep, dx)
            
            % Find boundary points, padding points for given shape
            [obj, xPad, yPad] = obj.padding(eraseSize, velPadding, lineStep);
            
            % Calculate V on padding, average back to line
            [obj.VLine, Verror] = obj.calcV(xPad, yPad, E, EA, epsBgnd);
            if(Verror > 1) % set this to some tolerance!!
                fprintf('Error too high \n');
            end
            
            [numMon,numFreq,numUD,~] = size(obj.VLine);

            % Actual derivatives with respect to parameters
            % Factor of dx so that dFdxr in mesh-units
            dFdx = zeros(numMon,numFreq,numUD,4);
            for m=1:numMon
                for f=1:numFreq
                    for u=1:numUD
                        % Calculate derivative w/ respect to given parameters
                        dFdx(m,f,u,1) = dx^2 * lineStep * sum(obj.VLine(obj.right) - obj.VLine(obj.left));
                        dFdx(m,f,u,2) = dx^2 * lineStep * sum(obj.VLine(obj.top) - obj.VLine(obj.bottom));
                        dFdx(m,f,u,3) = dx^2 * lineStep * sum(obj.VLine(obj.right) + obj.VLine(obj.left));
                        dFdx(m,f,u,4) = dx^2 * lineStep * sum(obj.VLine(obj.top) + obj.VLine(obj.bottom));
                    end
                end
            end
            
            % normalize with area differentials
            dAdx = obj.yLength;
            dAdy = obj.xLength;
            
            obj.dx0 = obj.dFdx0 / dAdx;
            obj.dy0 = obj.dFdy0 / dAdy;
            obj.dxLength = obj.dFdxLength / dAdx;
            obj.dyLength = obj.dFdyLength / dAdy;
            
            % Enforce any rules regarding parameters / points needing to
            % stay fixed
            %obj = enforceRules(obj);
            obj.dxr = 2 * obj.dxLength;
            obj.dyr = 2 * obj.dyLength;
        end

        %% RETURN CONSTRAINTS
        function[dA] = deltaArea(obj,deltaX)
            deltaX = abs(deltaX);
            dA = obj.yLength*deltaX(1) + obj.xLength*deltaX(2) ...
                + obj.yLength*deltaX(3) + obj.xLength*deltaX(4);
        end
        
        %% NORMALIZE DERIVATIVES
        function obj = normDerivs(obj, maxMove, normDParams, dx)
            obj.dx0 = maxMove * obj.dx0 / normDParams;
            obj.dy0 = maxMove * obj.dy0 / normDParams;
            obj.dxLength = (2*maxMove) * obj.dxLength / normDParams;
            obj.dyLength = (2*maxMove) * obj.dyLength / normDParams;
            
            obj.x0_new = obj.x0 + obj.dx0;
            obj.y0_new = obj.y0 + obj.dy0;
            obj.xLength_new = obj.xLength + obj.dxLength;
            obj.yLength_new = obj.yLength + obj.dyLength;
            
            % Calculate dF
            obj.dF = dx*(obj.dx0*obj.dFdx0 + obj.dy0*obj.dFdy0 ...
                    +obj.dxLength*obj.dFdxLength + obj.dyLength*obj.dFdyLength);
                
            % normalize with area differentials
            dAdx = obj.yLength;
            dAdy = obj.xLength;
            
            % Potentially try to figure out actual new shape - curr shape
            obj.dA = dAdx * abs(obj.dx0) + dAdy * abs(obj.dy0) + dAdx * abs(obj.dxLength) + dAdy * abs(obj.dyLength);
             
        end
        
        %% PAD THE RECTANGLE TO DELETE VELOCITY DATA AND INTERPOLATE
        function [obj, xPad, yPad] = padding(obj, eraseSize, velPadding, lineStep)
            margin = eraseSize + velPadding - 1;
            
            % define indices for rectangles's sides
            numLeft = floor((obj.yLength - 2*margin)/lineStep) + 1;
            numTop = floor((obj.xLength - 2*margin)/lineStep) + 1;
            obj.left = 1:numLeft;
            obj.top = numLeft + (1:numTop) ;
            obj.right = numLeft + numTop + (1:numLeft);
            obj.bottom = 2*numLeft + numTop + (1:numTop);
            
            % determine boundary points
            XLine(obj.left) = obj.x0 - .5*obj.xLength;
            XLine(obj.top) = linspace(obj.x0 - .5*obj.xLength + margin,obj.x0 + .5*obj.xLength - margin, numTop);
            XLine(obj.right) = obj.x0 + .5*obj.xLength;
            XLine(obj.bottom) = linspace( obj.x0 + .5*obj.xLength - margin, obj.x0 - .5*obj.xLength + margin, numTop);
            
            YLine(obj.left) = linspace(obj.y0 - .5*obj.yLength + margin, obj.y0 + .5*obj.yLength - margin, numLeft);
            YLine(obj.top) = obj.y0 + .5*obj.yLength;
            YLine(obj.right) = linspace(obj.y0 + .5*obj.yLength - margin, obj.y0 - .5*obj.yLength + margin, numLeft);
            YLine(obj.bottom) = obj.y0 - .5*obj.yLength;

            % determine normal vectors
            NXLine(obj.left) = -1;
            NXLine(obj.top) = 0;
            NXLine(obj.right) = +1;
            NXLine(obj.bottom) = 0;
            
            NYLine(obj.left) = 0;
            NYLine(obj.top) = +1;
            NYLine(obj.right) = 0;
            NYLine(obj.bottom) = -1;
            
            % normalize normal vectors (probably not needed)
            nLineNorm = sqrt(NXLine.^2+NYLine.^2);
            NXLine=NXLine./nLineNorm;
            NYLine=NYLine./nLineNorm;
            
            % determine points for padded boundaries
            xPad = meshgrid(XLine, 1:2*velPadding);
            xPad = xPad - [-(velPadding+eraseSize-1):-eraseSize , eraseSize:eraseSize+velPadding-1]'*NXLine;
            
            yPad = meshgrid(YLine, 1:2*velPadding);
            yPad = yPad - [-(velPadding+eraseSize-1):-eraseSize , eraseSize:eraseSize+velPadding-1]'*NYLine;
        
            % xLine,yLine,nxLine,nyLine will change size as the rectangle changes shape,
            % --> important to overwrite with fresh array
            obj.xLine = XLine;
            obj.yLine = YLine;
            obj.nxLine = NXLine;
            obj.nyLine = NYLine;
        end
        
        %% CALCULATE V
        function [VLine, Verror, ExPad, EyPad, epsPad, EpPad, DnPad] = calcV(obj, xPad, yPad, E, EA, epsBgnd)            
            % size of entire internal mesh
            xLen = size(E{1}{1},2);
            yLen = size(E{1}{1},1);
            szPad = size(xPad);

            [numMon, numFreq, numUser] = size(EA);
            VLine = zeros(numMon,numFreq,numUser,szPad(2));
            
            numPads = size(xPad,1);            
            numZ = size(epsBgnd{1}{1},3);
            nxLinePad = meshgrid(obj.nxLine,1:numPads);
            nyLinePad = meshgrid(obj.nyLine,1:numPads);
            nxLinePad = repmat(nxLinePad,[1 1 numZ]);
            nyLinePad = repmat(nyLinePad,[1 1 numZ]);
            
            for m = 1:numMon
                for f = 1:numFreq
                    epsi = epsBgnd{f}{1} + epsBgnd{f}{2} + epsBgnd{f}{3};
                    numEps = 1*(max(max(max(abs(epsBgnd{f}{1})))) > 0) + 1*( max(max(max(abs(epsBgnd{f}{2})))) > 0) + 1*(max(max(max(abs(epsBgnd{f}{3})))) > 0);
                    epsi = epsi / numEps;
                    for zInd = 1:numZ
                        for i = 1:numPads
                            epsPadR(i,:) = interp2(1:xLen, 1:yLen, real(epsi(:,:,zInd)), xPad(i,:), yPad(i,:));
                            epsPadI(i,:) = interp2(1:xLen, 1:yLen, imag(epsi(:,:,zInd)), xPad(i,:), yPad(i,:));
                        end
                        epsPad(:,:,zInd) = epsPadR + 1i*epsPadI;
                    end
                    epsPad(isnan(epsPad)) = 0;
                    epsInside = repmat( epsPad(end,:,:), szPad(1), 1); %szPad(1)/2+1
                    epsOutside = repmat( epsPad(1,:,:), szPad(1), 1); %szPad(1)/2
                    
                    for u=1:numUser
                        if( ~isempty(EA{m,f,u}) )
                            % interpolate fields and epsBgnd for each padded row
                            for zInd = 1:numZ
                                for i = 1:numPads
                                    ExPad(i,:,zInd) = interp2(1:xLen, 1:yLen, E{m,f,u}{1}(:,:,zInd), xPad(i,:), yPad(i,:));
                                    EyPad(i,:,zInd) = interp2(1:xLen, 1:yLen, E{m,f,u}{2}(:,:,zInd), xPad(i,:), yPad(i,:));
                                    EzPad(i,:,zInd) = interp2(1:xLen, 1:yLen, E{m,f,u}{3}(:,:,zInd), xPad(i,:), yPad(i,:));
                                    EAxPad(i,:,zInd) = interp2(1:xLen, 1:yLen, EA{m,f,u}{1}(:,:,zInd), xPad(i,:), yPad(i,:));
                                    EAyPad(i,:,zInd) = interp2(1:xLen, 1:yLen, EA{m,f,u}{2}(:,:,zInd), xPad(i,:), yPad(i,:));
                                    EAzPad(i,:,zInd) = interp2(1:xLen, 1:yLen, EA{m,f,u}{3}(:,:,zInd), xPad(i,:), yPad(i,:));
                                end
                            end
                            
                            % calculate par and perp components
                            EpPad = nyLinePad.*ExPad - nxLinePad.*EyPad;
                            EApPad = nyLinePad.*EAxPad - nxLinePad.*EAyPad;
                            DnPad = epsPad.*(nxLinePad.*ExPad + nyLinePad.*EyPad);
                            DAnPad = epsPad.*(nxLinePad.*EAxPad + nyLinePad.*EAyPad);
                            
                            % calculate V for each padded row
                            VPad = 2*real( ...
                                (epsInside - epsOutside) .* EpPad .* EApPad ... %par component
                                + (epsInside - epsOutside) .* EzPad .* EAzPad ... %par component
                                + (1./epsOutside - 1./epsInside) .* DnPad .* DAnPad ); %perp component
                            
                            VLine(m,f,u,:) = mean( mean(VPad,1) , 3); % 1 = mean over padded region, 3 = mean over z
                            Verror = 0;%max( std(VPad, 1, 1) ./ mean((mean(abs(VPad)))) );
                        end
                    end
                end
            end
        end
        
        %% ENFORCE RULES
        function deltaX = enforceRules(obj, deltaX)            
            % 7 rules
            % fix x0, y0
            % fix to a square
            % fix four possible sides
            Dx0 = deltaX(1); 
            Dy0 = deltaX(2);
            DxLength = deltaX(3);
            DyLength = deltaX(4);
            
            % rule 1: fix x0
            if(obj.rules(1))
                deltaX(1) = 0;
            end
            
            % rule 2: fix y0
            if(obj.rules(2))
                deltaX(2) = 0;
            end
            
            % rule 3: fix rectangle to a square
            if(obj.rules(3))
                deltaX(3) = 0.5*(DxLength+DyLength);
                deltaX(4) = deltaX(3);
            end
            
            % rule 4: fix left side
            if( obj.rules(4) )
                deltaX(1) = 0.5*(Dx0+DxLength/2);
                deltaX(3) = deltaX(1)*2;
            end
            
            % rule 5: fix bottom side
            if( obj.rules(5) )
                deltaX(2) = 0.5*(Dy0+DyLength/2);
                deltaX(4) = deltaX(2)*2;
            end
            
            % rule 6: fix right side
            if( obj.rules(6) )
                deltaX(1) = 0.5*(Dx0-DxLength/2);
                deltaX(3) = -deltaX(1)*2;
            end
            
            % rule 7: fix top side
            if( obj.rules(7) )
                deltaX(2) = 0.5*(Dy0-DyLength/2);
                deltaX(4) = -deltaX(2)*2;
            end
            
            % !!! Add rules for xLength, yLength, aspect ratio
            
        end
        
        %% CHECK IF RECTANGLE TOO SMALL (IN WHICH CASE MAYBE DELETE)
        function flag = tooSmall(obj, smallSize)
            flag = 0;
            if (obj.xLength < 2*smallSize  || obj.yLength < 2*smallSize)
                flag = 1;
            end
        end
        
        %% GET DRAWING OF CURRENT RECTANGLE
        % padding in units of internal mesh
        function mask = getDrawing(obj, xLengthGeo, yLengthGeo, padding)
            xRange = max( round(obj.x0-.5*(obj.xLength+2*padding)), 1 ) ...
                : min( round(obj.x0+.5*(obj.xLength+2*padding)), xLengthGeo) ;
            yRange = max( round(obj.y0-.5*(obj.yLength+2*padding)), 1 ) ...
                : min( round(obj.y0+.5*(obj.yLength+2*padding)), yLengthGeo );

            mask = zeros(yLengthGeo,xLengthGeo);
            mask(yRange,xRange) = 1;
        end

        %% GET DRAWING OF NEW RECTANGLE
        function mask = getNewDrawing(obj, xLengthGeo, yLengthGeo, padding, deltaX)
            x0 = obj.x0 + deltaX(1);
            y0 = obj.y0 + deltaX(2);
            xLength = obj.xLength + deltaX(3);
            yLength = obj.yLength + deltaX(4);
            
            xRange = round( max(1,x0-.5*(xLength+2*padding)) : min(xLengthGeo, x0+.5*(xLength+2*padding)) );
            yRange = round( max(1,y0-.5*(yLength+2*padding)) : min(yLengthGeo, y0+.5*(yLength+2*padding)) );
            
            mask = zeros(yLengthGeo,xLengthGeo);
            mask(yRange,xRange) = 1;
        end
        
    end
    
end
