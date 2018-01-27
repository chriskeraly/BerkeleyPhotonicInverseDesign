classdef Optimizer
    % Controls an optimization from start to finish
    
    %% PROPERTIES
    properties
        % Basic flags to simulation engine
        simFlag; % 1 = forward simulations, 2 = adjoint simulations
        notFinished;
        
        % basic optimization properties
        iter = 0; % current iteration (default = 0)
        numIter;
        alpha; % parameter on the order of 5%, allowable decrease in merit function
        testFlag;
        numUserSims; % = 0 if want to use already-defined geometry
        tempStep = 0.2; % temporary step size factor (on order of 0.2)
        longStep = 0.9; % long-term step size factor (on order of 0.9)
        prevFOM; % Keeps track of FOM that needs to improve
        prevGeo;
        bndUpdate; % 1 = moved boundary, 0 = new shape
        reRun = 0;
        reRunLimit = 2;
        
        % frequencies
        %   Should be an array of size [Nx3], where each row is of the form
        %   [f1 f2 nf], i.e. startFreq, endFreq, numFreq's
        %   N = number of different-frequency simulations to run
        freq;
        
        % (Current) geometry
        geo;
        
        % (Current) gradient calculation
        grad;
        
        % Merit function (May want to change the data structure)
        mf;
        dFdx;
        deltaXorig;
        deltaForig;
        deltaX;
        deltaFpred;
        
        % Binary variable dictating whether to try to approximate many
        % adjoint frequencies in a single simulation or whether to just do
        % them separately
        approxFreq;
        
        % Cell arrays that contain the relevant data for each simulation to
        %       run.  Size of arrays corresponds to number of simulations.
        % forwardSolves is of size [1xMxN], where M is the number of
        %       frequencies to do and N is the number of userSims to
        %       do.  Should never have to do more than 1 sim for monitors
        % adjointSolves is of size [OxPxQ], O is the number of monitor sims
        %       to do, and P and Q are the number of userSims and monSims
        %       to do (will be different from M and N if approxFreq = 1).
        % These are definable when the merit function is known
        forwardSolves;
        adjointSolves;
        skipInd;
        
        % Data to Save
        % Required: FHist, FOMHist
        % possibilities: EHist, EAHist, GradHist, GeoHist
        % Note that if GeoHist is saved, the field data within the geo
        %       object is NOT!
        data;
        
        % Plots to make
        % [FOMHist,FHist,Geometry,dFHist]
        plotFlags;
        
    end
    
    %% METHODS
    methods
        %% CONSTRUCTOR
        function[obj] = Optimizer(geo,grad,mf,freq,approxFreq,numUserSims,dataList,numIter,alpha,testFlag)
                  
            if(~all(freq(:,2)>=freq(:,1)))
                error('freq must be in increasing order');
            end
            if( (freq(1,1)~=freq(1,2)) && (freq(1,3)<=1) )
                error('numFreq must be >1 if using a finite frequency range');
            end
            if( (freq(1,1)==freq(1,2)) && (freq(1,3)>1) )
                error('must specify a finite frequency range if numFreq >1');
            end
            if( ~isa(geo,'FreeForm') && strcmp(mf.operators{end},'min'))
                error('min in Merit Function is only implemented for FreeForm Geometry objects. Use minFast instead.');
            end
            
            % Signify just created Optmizer
            obj.simFlag = 0;
            obj.iter = 0;
            
            % Attach geometry, merit function, etc.
            obj.geo = geo;
            obj.grad = grad;
            obj.mf = mf;
            obj.freq = freq;
            obj.approxFreq = approxFreq;
            obj.numUserSims = numUserSims;
            obj.numIter = numIter;
            obj.alpha = alpha;
            obj.testFlag = testFlag;
            obj.prevFOM = [];
            
            % Now that merit function is known, discern forward/adjoint
            % solves
            obj = obj.calcSolves;
            
            numFreq = sum(freq(:,3));
            obj.mf = obj.mf.setOptData(testFlag,numFreq,max(numUserSims,1));
            
            numMon = size(obj.adjointSolves,1);
            obj.grad = obj.grad.setOptData(testFlag,numMon,numFreq,max(numUserSims,1));
            
            % Create data structure for data to save
            obj.data = struct();
            obj.data.FOMHist = zeros(1,obj.numIter+1);
            obj.data.FHist = zeros(obj.mf.numMonitors, obj.mf.numFreq, obj.mf.numUser,obj.numIter+1);
            obj.data.dFOMactual = zeros(1,obj.numIter);
            obj.data.dFOMpred = zeros(1,obj.numIter);
            obj.data.dFactual = zeros(obj.mf.numMonitors, obj.mf.numFreq, obj.mf.numUser,obj.numIter+1);
            obj.data.dFpred = obj.data.dFactual;
            % At some point include deltaFpred and dFa
            % obj.data.dFa = zeros(obj.mf.numMonitors, obj.mf.numFreq, obj.mf.numUser,obj.numIter);
            for i=1:length(dataList)
                str = dataList{i};
                if( strcmp(str,'GradHist') )
                    obj.data.GradHist = cell(1,obj.numIter);
                elseif( strcmp(str,'EHist') )
                    obj.data.EHist = cell(1,obj.numIter+1);
                elseif( strcmp(str,'EAHist') )
                    obj.data.EAHist = cell(1,obj.numIter);
                elseif( strcmp(str,'GeoHist') )
                    obj.data.GeoHist = cell(1,obj.numIter+1);
                else
                    warning('Invalid data to save option: %s.  Ignoring...\n',str);
                end
            end
            
            % Determine index to skip
            obj.skipInd = 0;
            if(strcmp(obj.mf.operators{3},'min'))
                obj.skipInd = obj.mf.operatorOrder(3);
            end
        end
        
        %% RUN OPTIMIZATION HALF-STEP
        % Return data to simulation engine
        function[obj,simFlag,notFinished] = runOptStep(obj)
            notFinished = 1;
            if(obj.simFlag==0) % Just created optimizer
                simFlag = 1;
            elseif(obj.simFlag==1) % Just ran forward simulations               
                obj = obj.calcMerit;
                simFlag = 2;
                if(~isempty(obj.prevFOM))
                    if( (obj.alpha~=1) && (obj.mf.FOM < (1-obj.alpha)*obj.prevFOM) )
                        obj.reRun = obj.reRun + 1;
                        if(obj.reRun <= obj.reRunLimit)
                            if(obj.bndUpdate)
                                obj.deltaX = obj.tempStep * obj.deltaX;
                                obj.prevGeo = obj.prevGeo.changeStepSize(.5);
                            elseif(isa(obj.geo,'FreeForm'))
                                obj.prevGeo.newShapeRad = .5*obj.prevGeo.newShapeRad;
                            end
                            obj = obj.updateGeometry;
                            simFlag = 1; % Re-run a forward simulation
                            fprintf('Reducing step size, re-running simulation \n');
                        end
                    elseif( (obj.mf.FOM == obj.prevFOM))
                        obj.reRun = obj.reRun + 1;
                        if(obj.reRun <= obj.reRunLimit)
                            if(obj.bndUpdate)
                                obj.deltaX = obj.deltaX / obj.tempStep;
                                obj.prevGeo = obj.prevGeo.changeStepSize(1.5);
                            elseif(isa(obj.geo,'FreeForm'))
                                obj.prevGeo.newShapeRad = obj.prevGeo.newShapeRad*1.5;
                            end
                            obj = obj.updateGeometry;
                            simFlag = 1; % Re-run a forward simulation
                            fprintf('Increasing step size, re-running simulation \n');
                        end
                    end
                end
                if(simFlag==2)
                    obj.reRun = 0;
                    if(obj.prevFOM == obj.mf.FOM)
                        % End optimization if FOM has not changed
                        notFinished = 0;
                        simFlag = 1;
                    end
					if(isempty(obj.prevFOM) || (obj.prevFOM < obj.mf.FOM))
                        % Only update prevFOM if FOM has improved
                        % Allows for proper step size reductions
                        obj.prevFOM = obj.mf.FOM;
                    end
                    
                    obj.prevGeo = obj.geo;
                    obj = obj.saveOptData;
                    obj = obj.generatePlots;
                    obj.iter = obj.iter + 1;
                    if(obj.iter > obj.numIter)
                        notFinished = 0;
                        simFlag = 1;
                    end
                end
            elseif(obj.simFlag==2) % Just ran adjoint simulations
                simFlag = 1;
                obj = obj.calcDerivs;
                obj.prevGeo = obj.geo; % Make sure prevGeo also gets derivs
                obj = obj.updateGeometry;
                obj.grad = obj.grad.initData;
            end
            obj.simFlag = simFlag;
            obj.notFinished = notFinished;
        end
        
        %% RESUME A PREVIOUS OPTIMIZATION
        % iterNo currently defined as the iterNo+1 geometry (i.e. iterNo = 10
        % refers to the 11th geometry, which was created by the 10th
        % iteration)
        function[obj] = resume(obj,iterNo,numNewIter)
            obj.simFlag = 0;
            obj.numIter = iterNo + numNewIter;
            obj.iter = iterNo;
            obj.geo = obj.data.GeoHist{iterNo};
            obj.grad = obj.grad.initData;
            obj.mf = obj.mf.initData;
            obj.prevFOM = [];
            
            % Store current values of data
            FOMHistS = obj.data.FOMHist;
            FHistS = obj.data.FHist;
            dFOMa = obj.data.dFOMactual;
            dFOMp = obj.data.dFOMpred;
            dFa = obj.data.dFactual;
            dFp = obj.data.dFpred;
            obj.data.FOMHist = [FOMHistS(1:iterNo), zeros(1,numNewIter+1)];
            obj.data.FHist = zeros(obj.mf.numMonitors, obj.mf.numFreq, obj.mf.numUser,obj.numIter+1);
            obj.data.FHist(:,:,:,1:iterNo) = FHistS(:,:,:,1:iterNo);
            obj.data.dFOMactual = [dFOMa(1:iterNo), zeros(1,numNewIter+1)];
            obj.data.dFOMpred = [dFOMp(1:iterNo), zeros(1,numNewIter+1)];
            obj.data.dFactual = zeros(obj.mf.numMonitors, obj.mf.numFreq, obj.mf.numUser,obj.numIter+1);
            obj.data.dFactual(:,:,:,1:iterNo) = dFa(:,:,:,1:iterNo);
            obj.data.dFpred = zeros(obj.mf.numMonitors, obj.mf.numFreq, obj.mf.numUser,obj.numIter+1);
            obj.data.dFpred(:,:,:,1:iterNo) = dFp(:,:,:,1:iterNo);
            % At some point include deltaFpred and dFa
            
            if(isfield(obj.data,'EHist'))
                EHistS = obj.data.EHist;
                obj.data.EHist = cell(1,obj.numIter+1);
                obj.data.EHist(1:iterNo) = EHistS(1:iterNo);
            end
            if(isfield(obj.data,'GradHist'))
                GradHistS = obj.data.GradHist;
                obj.data.GradHist = cell(2,obj.numIter);
                obj.data.GradHist(1:iterNo) = GradHistS(1:iterNo);
            end
            if(isfield(obj.data,'EAHist'))
                EAHistS = obj.data.EAHist;
                obj.data.EAHist = cell(1,obj.numIter);
                obj.data.EAHist(1:iterNo) = EAHistS(1:iterNo);
            end
            if(isfield(obj.data,'GeoHist'))
                GeoHistS = obj.data.GeoHist;
                obj.data.GeoHist = cell(1,obj.numIter+1);
                obj.data.GeoHist(1:iterNo) = GeoHistS(1:iterNo);
            end

        end
                
        %% CALCULATE WHAT FORWARD AND ADJOINT SIMULATIONS TO RUN
        function[obj] = calcSolves(obj)
            
            % First figure out forward solves
            numFreqSims = size(obj.freq,1);
            
            obj.forwardSolves = cell(1,numFreqSims,max(obj.numUserSims,1));

            for j = 1:numFreqSims
                freqInd = sum(obj.freq(1:j-1,3))+1:sum(obj.freq(1:j,3)); % indices among total number of frequencies
                freqVec = linspace(obj.freq(j,1),obj.freq(j,2),obj.freq(j,3)); % actual frequency values to simulate with
                for k = 1:max(obj.numUserSims,1)
                    obj.forwardSolves{1,j,k} = struct(...
                        'freqVec',freqVec,...
                        'runUserSim',1*(obj.numUserSims>0),...
                        'freqInd',freqInd);
                end
            end
            
            % Now figure out adjoint solves
            if(obj.approxFreq)
                numFreqSims = size(obj.freq,1);
            else
                numFreqSims = sum(obj.freq(:,3));
            end
            
            % Figure out whether multiple sims are needed for different monitors 
            numMonSims = 1;
            monInd = 1:obj.mf.numMonitors;
            if( strcmp(obj.mf.operators{3},'min') )
                if( obj.mf.operatorOrder(3) == 1 )
                    numMonSims = obj.mf.numMonitors;
                    monInd = (1:obj.mf.numMonitors).';
                end
            end
            obj.adjointSolves = cell(numMonSims,numFreqSims,max(obj.numUserSims,1));
            
            for i = 1:numMonSims
                for j = 1:numFreqSims
                    if(obj.approxFreq)
                        freqInd = sum(obj.freq(1:j-1,3))+1:sum(obj.freq(1:j,3));
                        freqVec = linspace(obj.freq(j,1),obj.freq(j,2),obj.freq(j,3));
                    else
                        freqInd = j;
                        a = cumsum(obj.freq(:,3),1);
                        [freqSubInd,freqRow] = min(j-a(j-a>0));
                        if(isempty(freqSubInd))
                            freqRow = 0;
                            freqSubInd = j;
                        end
                        freqRow = freqRow + 1;
                        fmin = obj.freq(freqRow,1);
                        fmax = obj.freq(freqRow,2);
                        Nf = obj.freq(freqRow,3);
                        if(Nf>1)
                            freqVec = fmin + (freqSubInd-1)*(fmax-fmin)/(Nf-1);
                        else
                            freqVec = fmin;
                        end
                    end
                    
                    for k = 1:max(obj.numUserSims,1)
                        obj.adjointSolves{i,j,k} = struct( ...
                            'freqInd',freqInd, ...
                            'freqVec',freqVec, ...
                            'monInd',monInd(i,:), ...
                            'runUserSim',1*(obj.numUserSims>0) );
                    end
                end
            end
        end
        
        %% UPDATE MERIT FUNCTION FIELD DATA
        function[obj] = updateMeritData(obj, merit_E,merit_H,merit_pos,merit_eps,merit_eps_pos,normParam,freqInd, userSimNo)
            obj.mf = obj.mf.updateData(merit_E,merit_H,merit_pos,merit_eps,merit_eps_pos,normParam,freqInd, userSimNo);
        end
        
        %% UPDATE GEOMETRY FIELD DATA
        function[obj] = updateGeoData(obj, E, pos_E, eps, pos_eps, monSimNo, freqInd, userSimNo, EAflag)
            obj.grad = obj.grad.updateFieldData(E, pos_E{1}, pos_E{2}, pos_E{3}, monSimNo, freqInd, userSimNo, EAflag);
            if(EAflag == 0)
                obj.grad = obj.grad.updateEpsData(eps, pos_eps{1}, pos_eps{2}, pos_eps{3}, freqInd);
            end
        end
        
        %% CALCULATE OVERALL MERIT FUNCTION
        function[obj] = calcMerit(obj)
            obj.mf = obj.mf.calcMeritAndDipoles;
        end
        
        %% CALCULATE DERIVATIVES AND DESIRED UPDATES
        function[obj] = calcDerivs(obj)
            [epsGrid, ~, ~, ~, ~, ~, ~, ~] = obj.geo.getGeometry;
            
            obj.grad = obj.grad.calcDerivs(obj.skipInd,obj.freq,epsGrid);
            dFdxi = obj.grad.dFdxBnd;
            
            if(isa(obj.geo,'FreeForm'))
                [lb,ub,c] = obj.geo.getConstraintParams;
            else
                lb=-1;
                ub=1;
                c=0;
            end
            
            for i=1:3
                if(i~=obj.skipInd), dFdxi = sum(dFdxi,i); end
            end
            
            if(obj.skipInd==1)
                perm = [1 4 5 6 2 3];
            elseif(obj.skipInd==2)
                perm = [2 4 5 6 1 3];
            else
                perm = [3 4 5 6 1 2];
            end
            
            obj.dFdx = permute(dFdxi,perm); % get rid of other mfu dimensions
                        
            numMin = size(obj.dFdx,1);
            dim1 = size(obj.dFdx,2);
            dim2 = size(obj.dFdx,3);
            dim3 = size(obj.dFdx,4);
            
            if(all(all(all(~obj.dFdx))))
                obj.deltaX = zeros(dim1,dim2);
            else
                if(numMin > 1)
                    if((dim2 > 1) || (dim3 > 1))
                        [xInd,yInd] = find(permute(obj.dFdx(1,:,:,:), [2 3 4 1]));
                        ind = sub2ind([dim1,dim2,dim3], xInd, yInd);
                        lb = lb(ind);
                        ub = ub(ind);
                        for i=1:numMin
                            [~,~,dFdx_Comp(i,:)] = find(permute(obj.dFdx(i,:,:,:), [2 3 4 1]));
                        end
                    else
                        dFdx_Comp = obj.dFdx;
                    end
                    
                    %figure(203); plot(1:length(dFdx_Comp(1,:)), dFdx_Comp(1,:), 1:length(dFdx_Comp(1,:)), dFdx_Comp(2,:));
                    
                    tic;
                    options = optimset('MaxIter',100,'MaxFunEvals',1e5,'Display','notify','PlotFcns',{@optimplotfval,@optimplotconstrviolation,@optimplotstepsize});
                    [deltaX,obj.deltaFpred] = minimaxOpt(obj.mf.Fmin,dFdx_Comp,lb,ub,c,options);
                    toc
                    
                    %figure(204); plot(deltaX);
                    
                    if((dim2 > 1) || (dim3 > 1))
                        obj.deltaX = zeros(dim1,dim2,dim3);
                        obj.deltaX(ind) = deltaX;
                    else
                        obj.deltaX = deltaX;
                    end
                    
                else
                    % more sophisticated update?
                    obj.dFdx = obj.dFdx;
                    obj.deltaX = obj.dFdx;
                    if( (~isempty(obj.deltaX))&&(max(max(abs(obj.deltaX)))~=0) )
                        obj.deltaX = obj.deltaX * max(max(abs(ub)),max(abs(lb)))/max(max(abs(obj.deltaX)));
                    end
                    obj.deltaFpred = sum(sum(obj.dFdx.*obj.deltaX));
                    obj.deltaX = permute(obj.deltaX, [2 3 4 1]); % clever way to get rid of first dimension (min parameter dimension)
                    % A simple squeeze(deltaX) could mess up the dimensionality of the rest of the deltaX matrix
                end
            end
        end
            
        %% UPDATE GEOMETRY
        function[obj] = updateGeometry(obj)
            if(isa(obj.geo,'FreeForm'))
                [obj.geo, obj.deltaX, obj.deltaFpred, obj.bndUpdate] = ...
                    obj.prevGeo.updateShapes(obj.deltaX, obj.mf.Fmin, obj.dFdx, obj.grad.dFdxSpace);
            else
                [obj.geo, obj.deltaFpred] = ...
                    obj.prevGeo.updateShapes(obj.mf.Fmin, obj.dFdx, obj.grad.dFdxSpace);
                obj.deltaX = zeros(size(obj.dFdx));
                obj.bndUpdate = 1;
            end
            obj.data.dFOMpred(obj.iter) = obj.deltaFpred;
        end
            
        %% GET GEOMETRY DATA (For Import to Lumerical)
        function[geoData] = exportGeometry(obj)
            
            if(isa(obj.geo,'FourierSurface') || isa(obj.geo,'Polygons'))
                geoData = obj.geo.returnData;
            else
                [epsGrid, eps, epsOut, x0, y0, z0, dx, thickness] = obj.geo.getGeometry;
                
                [numY, numX] = size(epsGrid);
                numZ = floor(1e-4*round(1e4*thickness / dx));
                [xGrid0, yGrid0] = meshgrid(x0+dx*(0:numX-1)+dx/2, y0+dx*(0:numY-1)+dx/2);
                if(numZ>0)
                    zVec = z0 + dx*(-numZ/2:numZ/2);
                else
                    zVec = z0 + dx*[-1 1];
                end
                
                % Check if epsGrid is binary valued
                binaryFlag = all(all((epsGrid==0) | (epsGrid==1)));
                
                % FreeForm Binary Import
                if( binaryFlag && (numX>1) && isreal(eps) && isreal(epsOut) )
                    [xGrid, yGrid] = meshgrid(1:numX, 1:numY);
                    cnt = 0;
                    yArr=[];
                    xArr=[];
                    blockSizeArr=[];
                    
                    if(numY==1)
                        bndIn = getBoundaryIn(epsGrid);
                        ind = find(bndIn);
                        x1 = xGrid0(ind(1:2:end-1));
                        x2 = xGrid0(ind(2:2:end));
                        % eps center coordinates and square dimensions
                        xArr = .5*(x2+x1);
                        yArr = yGrid(1,1);
                        blockSizeArr = x2-x1;
                    else
                        for blockSize = 100:-2:0 %min(numX,numY):-2:0
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
                        % eps center coordinates and square dimensions
                        xArr = xGrid0(1,xArr);
                        yArr = yGrid0(yArr,1);
                        blockSizeArr = blockSizeArr*dx;
                    end
                    
                    % epsOut center coordinates and dimensions
                    geo_x = mean([xGrid0(1,1) xGrid0(1,end)]);
                    geo_xspan = (numX)*dx;
                    geo_y = mean([yGrid0(1,1) yGrid0(end,1)]);
                    geo_yspan = (numY)*dx;
                    % center z coordinate and extrusion depth
                    geo_z = z0;
                    geo_zspan = (numZ) * dx;
                    
                    if( ischar(eps) )
                        geo_n = 0;
                        geo_mat = eps;
                    else
                        geo_n = sqrt(eps);
                        geo_mat = '';
                    end
                    if( ischar(epsOut) )
                        geo_nClad = 0;
                        geo_matClad = epsOut;
                    else
                        geo_nClad = sqrt(epsOut);
                        geo_matClad = '';
                    end
                    
                    geo_type = 'binary';
                    geoData = struct('geo_type',geo_type,'xArr',xArr,'yArr',yArr,'blockSizeArr',blockSizeArr,'geo_n',geo_n,'geo_mat',geo_mat,'geo_nClad',geo_nClad,'geo_matClad',geo_matClad,'geo_x',geo_x,'geo_xspan',geo_xspan,'geo_y',geo_y,'geo_yspan',geo_yspan,'geo_z',geo_z,'geo_zspan',geo_zspan);
                    
                    
                    % FreeForm importnk Matrix
                elseif( ~ischar(eps) && ~ischar(epsOut) )
                    n = sqrt(epsGrid*(eps-epsOut)+epsOut);
                    n = n.'; % transpose for Lumerical
                    zVec = [zVec(1),zVec(end)];
                    len = length(zVec);
                    n = repmat(n,[1 1 len]);
                    x = xGrid0(1,:);
                    y = yGrid0(:,1);
                    z = zVec;
                    geo_type = 'importnk';
                    if(numY==1)
                        n = permute(n,[1 3 2]);
                        geoData = struct('geo_type',geo_type,'geo_n',n,'geo_x',x,'geo_y',z,'geo_z',y);
                    else
                        geoData = struct('geo_type',geo_type,'geo_n',n,'geo_x',x,'geo_y',y,'geo_z',z);
                    end
                    
                    % Error
                else
                    error('No export option available for material types. Try using materials from the Lumerical Material Database.');
                end
                
            end
            %             FreeForm_image Import
            %             if( ischar(obj.eps) )
            %                 geo_n = 0;
            %                 geo_mat = obj.eps;
            %             else
            %                 geo_n = sqrt(obj.eps);
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
        
        %% GET WEIGHTS
        function[weights] = getWeights(obj,i,j,k) % note that i,j,k may be vectors!
            weights = obj.mf.w(i,j,k);
        end
        
        %% GET DIPOLES
        function[dipoles] = getDipoles(obj,i,j,k) % note that i,j,k may be vectors!
            dipoles = obj.mf.dip(i,j,k);
        end
        
        %% GENERATE PLOTS
        function[obj] = generatePlots(obj)
            
            figCnt = 0;
            if(obj.plotFlags(1)) % FOMHist
                figCnt = figCnt + 1;
                figure(figCnt);
                hold off;
                plot(1:obj.iter+1,obj.data.FOMHist(1:obj.iter+1));
                title('Figure of Merit History','FontWeight','bold','FontSize',12,'FontName','Arial');
                xlabel('iteration');
            end
            if(obj.plotFlags(2)) % FHist
                sz = size(obj.data.FHist);
                if( (sz(1)+sz(2)+sz(3))>3 )
                    figCnt = figCnt + 1;
                    figure(figCnt);
                    clf;
                    cnt = 0;
                    for i=1:sz(1)
                        for k=1:sz(3)
                            cnt = cnt+1;
                            subplot(sz(1),sz(3),cnt);
                            hold all;
                            for j=1:sz(2)
                                plot(1:obj.iter+1, reshape(obj.data.FHist(i,j,k,1:obj.iter+1),1,[]) );
                            end
                            xlabel('iteration');
                            if(sz(3)>1)
                                title(['Monitor ',num2str(i),': Custom Parameter ',num2str(k)],'FontWeight','bold','FontSize',12,'FontName','Arial');
                            else
                                title(['Monitor ',num2str(i)],'FontWeight','bold','FontSize',12,'FontName','Arial');
                            end
                            hold off;
                        end
                    end
                end
            end
            if(obj.plotFlags(3)) % Geometry
                figCnt = figCnt + 1;
                figure(figCnt);
                [epsGrid, ~, ~, ~, ~, ~, ~, ~] = obj.geo.getGeometry;
                drawing = real(epsGrid * (obj.grad.epsVec(1)-obj.grad.epsOutVec(1)) + obj.grad.epsOutVec(1));
                %drawing = real(obj.geo.getDrawing(0));
                if(sum(size(drawing)>1)>1)
                    imagesc(0:(size(drawing,2)-1),0:(size(drawing,1)-1),drawing); axis equal;
                else
                    plot(0:(length(drawing)-1),drawing(:));
                end
                title('Re[permittivity]','FontWeight','bold','FontSize',12,'FontName','Arial');
            end
            if(obj.plotFlags(4))
                if(obj.iter>0)
                    figCnt = figCnt + 1;
                    figure(figCnt);
                    clf;
                    plot(1:obj.iter,obj.data.dFOMpred(1:obj.iter),1:obj.iter,obj.data.dFOMactual(1:obj.iter));
                    title('Change in FOM','FontWeight','bold','FontSize',12,'FontName','Arial');
                    legend('predicted','actual');
                    xlabel('iteration');
                end
            end
                
%             posOld = posNew;
%             posNew = {geo.xLine, geo.yLine, geo.nxLine, geo.nyLine};
%             velocity = geo.VLine;
%             
%             figCnt = figCnt + 1;
%             figure(figCnt); clf;  set(gcf,'color','white');
%             imagesc(geo.velocity);
%             maxVal = max(max(abs( geo.velocity )));
%             caxis(maxVal*[-1,1]);
%             hold on; plot(geo.xLine,geo.yLine,'kx');
%             title('V_{bgnd}');
%             
%             figCnt = figCnt + 1;
%             figure(figCnt); clf; set(gcf,'color','white'); hold on;
%             plot(posOld{1}, posOld{2}, 'x');
%             plot(posNew{1}, posNew{2}, 'go');
%             xlim([1,geo.numX]); ylim([1,geo.numY]);
%             hold off; axis equal;
%             title('New (circle) and old (x) positions');
%             
%             figCnt = figCnt + 1;
%             figure(figCnt); clf; set(gcf,'color','white'); hold on;
%             plot(posOld{1}, posOld{2}, 'o');
%             quiver(posOld{1}, posOld{2}, posOld{3}.*velocity, posOld{4}.*velocity );
%             xlim([1,geo.numX]); ylim([1,geo.numY]);
%             hold off; axis equal;
%             title('Velocity and position');
%             
%             figCnt = figCnt + 1;
%             figure(figCnt); clf; set(gcf,'color','white'); hold on;
%             plot(posOld{1}, posOld{2}, 'o');
%             quiver(posOld{1}, posOld{2}, posOld{3}, posOld{4} );
%             xlim([1,geo.numX]); ylim([1,geo.numY]);
%             hold off; axis equal;
%             title('Normal vectors');
%             
%             figCnt = figCnt + 1;
%             figure(figCnt); clf; plot(FOMHist);
%             set(gcf,'color','white');
%             title('FOM History');
%             
%             if(iter>1)
%                 figCnt = figCnt + 1;
%                 figure(figCnt); plot(1:iter,dFHist(1:iter),1:iter,dFpredHist(1:iter));
%                 set(gcf,'color','white');
%                 title('deltaFOM History');
%             end
            
            % 	numDip = length(mf.dip{1}{1});
            % 	figCnt = figCnt + 1;
            % 	figure(figCnt); clf; plot(1:numDip,abs(pDip{1}{1}),1:numDip,abs(mDip{1}{1}));
            % 	set(gcf,'color','white');
            % 	title('x-dipole amplitudes');
            % 	figCnt = figCnt + 1;
            % 	figure(figCnt); clf; plot(1:numDip,abs(pDip{1}{2}),1:numDip,abs(mDip{1}{2}));
            % 	set(gcf,'color','white');
            % 	title('y-dipole amplitudes');
            % 	figCnt = figCnt + 1;
            % 	figure(figCnt); clf; plot(1:numDip,abs(pDip{1}{3}),1:numDip,abs(mDip{1}{3}));
            % 	set(gcf,'color','white');
            % 	title('z-dipole amplitudes');
        end
        
        %% SAVE DATA
        function[obj] = saveOptData(obj)
            obj.data.FOMHist(obj.iter+1) = obj.mf.FOM;
            obj.data.FHist(:,:,:,obj.iter+1) = obj.mf.F;
            if(obj.iter>0)
                obj.data.dFOMactual(obj.iter) = obj.data.FOMHist(obj.iter+1)-obj.data.FOMHist(obj.iter);
                obj.data.dFactual(:,:,:,obj.iter) = obj.data.FHist(:,:,:,obj.iter+1)-obj.data.FHist(:,:,:,obj.iter);
            end
            if(isfield(obj.data,'EHist'))
                obj.data.EHist{obj.iter+1} = obj.grad.E;
            end
            if(isfield(obj.data,'EAHist'))
                obj.data.EAHist{obj.iter+1} = obj.grad.EA;
            end
            if(isfield(obj.data,'GeoHist'))
                geoS = obj.geo;
                if(isa(geoS,'FreeForm'))
                    %geoS.newShapeGrid=[];
                end
                obj.data.GeoHist{obj.iter+1} = geoS;
            end
            if(isfield(obj.data,'GradHist'))
                if(obj.iter>0)
                    obj.data.GradHist(:,obj.iter) = {obj.grad.dFdxBnd};
                end
            end
        end
        
        %% SET PLOTTING FLAGS
        function[obj] = setPlotFlags(obj,plotFlags)
            obj.plotFlags = plotFlags;
        end
        
    end

end
       
