classdef MeritFunction
    %   Weights should be properties of actual MeritFunction classes, not
    %   this one!!
    %   Changed accordingly for Transmission (Needed for ModeMatch,
    %   FieldEnergy!)
    
    properties
        
        % contains all separate merit functions
        monitors;
        monNames;
        
        % Single Figure of Merit
        FOM;
        
        % Data structures for all merit and dipole data
        F;
        w;
        Fmin;
        dip;
		pos;
      
        % Sizes for important data structures
        % Note that numFreq is number of actual frequencies, not the number
        % of frequency simulations
        numFreq;
        numUser;
        numMonitors;
        
		% testFlag;
		testFlag;
        
        % 1 = monitor, 2 = frequency, 3 = user-defined
        operatorOrder;
        
        % functions implementing sum, prod, and min
        operators;
    end
    
    methods
        %% CONSTRUCTOR
        % operatorOrder and operators defined here
        % Also should specify details about: numFreq, numZ, numU
        function obj = MeritFunction(operatorOrder, operators)
            obj.monitors = {};
            obj.numMonitors = 0;
            obj.FOM = 0;
            obj.monNames = {};
            obj.operatorOrder=operatorOrder
%             obj.operatorOrder(ismember(operatorOrder,'monitor'))=1;
%             obj.operatorOrder(ismember(operatorOrder,'frequency'))=2;
%             obj.operatorOrder(ismember(operatorOrder,'user-defined'))=3;
            obj.operators = operators;
            
            if( length(operatorOrder)~=3 )
                error(' Incorrect operatorOrder length: %g \n',length(operatorOrder));
            end
            if( sum(obj.operatorOrder)~=6 )
                error([' Incorrect operatorOrder elements:' operatorOrder ' \n']);
            end
            if( length(operators)~=3 )
                error(' Incorrect operators length: %g \n',length(operators));
            end
            if(strcmp(operators{1},'min'))% || strcmp(operators{1},'minFast'))
                error(['@min operator appearing at index: 1.  @min is only',...
                    ' allowed to appear once, at index 3.  If fewer',...
                    ' than three operators are needed, pad with',...
                    ' dummy @sum operators.\n  Exiting...\n']);
            end
            if(strcmp(operators{2},'min'))% || strcmp(operators{2},'minFast'))
                error(['@min operator appearing at index: 2.  @min is only',...
                    ' allowed to appear once, at index 3.  If fewer',...
                    ' than three operators are needed, pad with',...
                    ' dummy @sum operators.\n  Exiting...\n']);                
            end
        end
        
        %% GET DATA CONSISTENT WITH OPTIMIZER
        function[obj] = setOptData(obj,testFlag,numFreq,numUser)
            obj.testFlag = testFlag;
            obj.numFreq = numFreq;
            obj.numUser = numUser;
            obj = obj.initData; % Can initialize data now
        end
        
        %% ADD A MERIT FUNCTION TO CURRENT OBJECT
        % newMerit = {type, name, weight, exponent, x0, y0, z0, xLength, yLength, zLength, dx, params}
        function obj = addMerit(obj,newMerit)
            obj.numMonitors = obj.numMonitors + 1;
            obj.monNames{obj.numMonitors} = newMerit{2};
            
            if( strcmp('transmission', newMerit{1}) )
                obj.monitors{obj.numMonitors} = Transmission( newMerit{3}, newMerit{4}, newMerit{5}, newMerit{6}, newMerit{7}, newMerit{8}, newMerit{9}, newMerit{10}, newMerit{11}, newMerit{12}, obj.testFlag );
            elseif( strcmp('modematch', newMerit{1}) )
                obj.monitors{obj.numMonitors} = ModeMatch( newMerit{3}, newMerit{4}, newMerit{5}, newMerit{6}, newMerit{7}, newMerit{8}, newMerit{9}, newMerit{10}, newMerit{11}, newMerit{12}, obj.testFlag );
            elseif( strcmp('fieldenergy', newMerit{1}) )
                obj.monitors{obj.numMonitors} = FieldEnergy( newMerit{3}, newMerit{4}, newMerit{5}, newMerit{6}, newMerit{7}, newMerit{8}, newMerit{9}, newMerit{10}, newMerit{11}, newMerit{12}, obj.testFlag );
            else
                error('Merit function monitor type must be transmission, modematch or fieldenergy.');
            end
        end
        
        %% SET/RESET F, PDIP, AND POS VARIABLES
        function[obj] = initData(obj)
            obj.F = zeros(obj.numMonitors, obj.numFreq, obj.numUser);
            obj.dip = cell(obj.numMonitors, obj.numFreq, obj.numUser);
            obj.pos = cell(obj.numMonitors,1);
        end
        
        %% UPDATE GIVEN NEW SIMULATION
        function[obj] = updateData(obj,merit_E,merit_H,merit_pos,merit_eps,merit_eps_pos, normParam, freqInd, usrIndex)
            for i=1:obj.numMonitors
                if(isa(obj.monitors{i},'FieldEnergy'))
                    [Fi,dipi,posi] = obj.monitors{i}.calcMerit(merit_E{i}, merit_H{i}, merit_pos{i}, merit_eps{i}, merit_eps_pos{i}, normParam, freqInd);
                else
                    [Fi,dipi,posi] = obj.monitors{i}.calcMerit(merit_E{i}, merit_H{i}, merit_pos{i}, normParam, freqInd);
                end
                obj.F(i,freqInd,usrIndex) = Fi;
                obj.dip(i,freqInd,usrIndex) = dipi;
                obj.pos(i) = posi;
            end
        end
        
        %% CALCULATE SINGLE OVERALL MERIT FUNCTION & DIPOLE WEIGHTS
        % F = merit function contracted by ind
        % w = weight matrix to multiply dipole arrays
        % Fmin is a list of the merit function contracted fully EXCEPT for
        %   dimension along which the min happens (if it does)
        function[obj] = calcMeritAndDipoles(obj)
            F = obj.F;
            w = ones(size(F));
            Fmin = [];
            
            indexHist = zeros(1,3);
            F0 = F;
            for i = 1:3
                index = obj.operatorOrder(i); 
                op = obj.operators{i};
                indexHist(index) = 1;
                
                if(obj.testFlag)
                    fprintf('Contracting merit on dimension: %g \n',index);
                    fprintf('Operator: %s%s \n','@',op);
                end
                
                if(strcmp(op,'prod'))
                    F = prod(F,index);
                    repeatMat = [size(F0,1),size(F0,2),size(F0,3)];
                    repeatMat( ~logical(indexHist) ) = 1;
                    
                    nw = repmat(F,repeatMat); % new weights
                    nw = nw./F0;
                    nw(isnan(nw)) = 0;
                    w = w.*nw;
                elseif(strcmp(op,'sum')) % Don't touch any of the dipole weights
                	F = sum(F,index);
                elseif(strcmp(op,'min')) % Don't touch any of the dipole weights
                    Fmin = reshape(F,[],1);
                    [F,ind] = min(F,[],index);
                elseif(strcmp(op,'minFast')) % Keep the dipole weights corresponding to min(F)
                    [F,ind] = min(F,[],index);
                    mask = zeros(size(w));
                    for i2 = 1:size(ind,1)
                        for i3 = 1:size(ind,2)
                            for j = 1:3
                                if(j<i)
                                    coord(obj.operatorOrder(j))={':'};
                                elseif(j==i)
                                    coord(obj.operatorOrder(j))={ind(i2,i3)};
                                elseif(j==i+1)
                                    coord(obj.operatorOrder(j))={i3};
                                else
                                    coord(obj.operatorOrder(j))={i2};
                                end
                            end
                            mask(coord{:}) = 1;
                        end
                    end
                    
					w = w .* mask;
                end
            end
            
            if(isempty(Fmin))
                Fmin = F;
            end
            
            obj.FOM = F
            obj.w = w;
            obj.Fmin = Fmin;
        end        
    end
    
end

