classdef ModeMatch < Monitor
    % ModeMatch class used for mode-matching merit function
 
    properties
        % Mode profile
        Em; % mode electric field on internal mesh
        Hm; % mode magnetic field on internal mesh (Em OR Hm specified)
        surface_normal; %1 = x, 2 = y, 3 = z
    end
    
    methods
        %% CONSTRUCTOR
        function obj = ModeMatch(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, params, testFlag)
            obj = obj@Monitor(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, testFlag);
            
            grid = params{1};
            EmGrid = params{2};
            HmGrid = params{3};
            obj.surface_normal=params{5};
            
            if (obj.surface_normal==0)
				error('Mode-Matching Monitor is not the right dimension');
			end
            if(obj.dim==3)
                error('Mode-Matching Monitor is 3D: not allowed.');
            end
            
            
            [eInt,hInt] = interpFields(obj,EmGrid,HmGrid,grid);
            obj.Em=eInt;
            obj.Hm = hInt;
           
            
			if(obj.testFlag)
				fprintf('%%%% ADDING MODE_MATCH MONITOR %%%% \n');
				fprintf('  (x0,y0,z0) = (%g,%g,%g) \n',x0,y0,z0);
				fprintf('  (numX,numY,numZ) = (%g,%g,%g) \n',obj.numX,obj.numY,obj.numZ);
				fprintf('  (eMode) = (%g) \n',obj.eMode);
			end
        end
        
        %% CALCULATE MERIT FUNCTION AND DIPOLE AMPLITUDES
        % Merit = Merit=abs(int(E*conj(Hm)).z)/abs(int(Em*conj(Hm).z)
		% *=vector product; .=scalar product
        function [F,dip,pos] =  calcMerit(obj, eField, hField, grid, powerNorm, freqInd)
            
            numFreq = length(freqInd);
            
            F = zeros(1,numFreq); 
            Ftmp=zeros(1,numFreq);
            dip = cell(1,numFreq);
            
            for i = 1:numFreq
                [eInt, hInt] = interpFields(obj,eField(i,:), hField(i,:), grid);
                wi = obj.weight(freqInd(i));
    
            	ex = eInt{1}; ey = eInt{2}; ez = eInt{3};
            	hx = hInt{1}; hy = hInt{2}; hz = hInt{3};
                    
            	if obj.surface_normal==1;
                    Fi = (ey.*conj(obj.Hm{3})-ez.*conj(obj.Hm{2}));
                    Fi2 = (conj(obj.Em{2}).*(hz)-conj(obj.Em{3}).*(hy));
                    modePow=abs(sum(sum(real((obj.Em{2}.*conj(obj.Hm{3})-obj.Em{3}.*conj(obj.Hm{2}))))))*obj.dx^obj.dim;
            	elseif obj.surface_normal==2;
                    Fi= (ez.*conj(obj.Hm{1})-ex.*conj(obj.Hm{3}));
                    Fi2=(conj(obj.Em{3}.*(hx)-conj(obj.Em{1}).*hz));
                    modePow=abs(sum(sum(real((obj.Em{3}.*conj(obj.Hm{1})-obj.Em{1}.*conj(obj.Hm{3}))))))*obj.dx^obj.dim;
            	elseif obj.surface_normal==3;
                    Fi=(ex.*conj(obj.Hm{2})-ey.*conj(obj.Hm{1}));
                    Fi2=(conj(obj.Em{1}.*(hy)-conj(obj.Em{2}).*hx));
                    modePow=abs(sum(sum(real((obj.Em{1}.*conj(obj.Hm{2})-obj.Em{2}.*conj(obj.Hm{1}))))))*obj.dx^obj.dim;
            	else
                	error('surface normal problem in mode match monitor');
            	end

	            F(1,i) = wi * 1/8* (abs(sum(sum(Fi+Fi2)* obj.dx^obj.dim)^2)/modePow/ powerNorm(i))^obj.exponent ;
    	        Ftmp(1,i) = wi * (sum(sum(Fi+Fi2)) * obj.dx^obj.dim );
    	        eCoeff = obj.exponent * abs(F(1,i))^(obj.exponent-1);

            
    	        if(obj.surface_normal==1)
                    pyDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * obj.eps0 *2* conj(obj.Hm{3})*conj(Ftmp(1,i))/modePow/8;
                    pzDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * obj.eps0 *2* (-conj(obj.Hm{2}))*conj(Ftmp(1,i))/modePow/8; 
                    pxDip = zeros(size(pyDip));
                    myDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim *2* obj.eps0/obj.mu0*(- conj(-obj.Em{3}))*conj(Ftmp(1,i))/modePow/8;
                    mzDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim *2* obj.eps0/obj.mu0*(- conj(obj.Em{2}))*conj(Ftmp(1,i))/modePow/8;
                    mxDip = zeros(size(pxDip));
    	        elseif(obj.surface_normal==2)
                    pzDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * obj.eps0 *2* conj(obj.Hm{1})*conj(Ftmp(1,i))/modePow/8;
                    pxDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * obj.eps0 *2* (-conj(obj.Hm{3}))*conj(Ftmp(1,i))/modePow/8; 
                    pyDip = zeros(size(pzDip));
                    mzDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim *2* obj.eps0/obj.mu0*(- conj(-obj.Em{1}))*conj(Ftmp(1,i))/modePow/8;
                    mxDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim *2* obj.eps0/obj.mu0*(- conj(obj.Em{3}))*conj(Ftmp(1,i))/modePow/8;
                    myDip = zeros(size(pzDip));       
    	        else
                    pxDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * obj.eps0 *2* conj(obj.Hm{2})*conj(Ftmp(1,i))/modePow/8;
                    pyDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * obj.eps0 *2* (-conj(obj.Hm{1}))*conj(Ftmp(1,i))/modePow/8; 
                    pzDip = zeros(size(pxDip));
                    mxDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim *2* obj.eps0/obj.mu0*(- conj(-obj.Em{2}))*conj(Ftmp(1,i))/modePow/8;
                    myDip = wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim *2* obj.eps0/obj.mu0*(- conj(obj.Em{1}))*conj(Ftmp(1,i))/modePow/8;
                    mzDip = zeros(size(pxDip));
            	end
                
                pxDip = obj.permuteLum(pxDip);
                pyDip = obj.permuteLum(pyDip);
                pzDip = obj.permuteLum(pzDip);
                mxDip = obj.permuteLum(mxDip);
                myDip = obj.permuteLum(myDip);
                mzDip = obj.permuteLum(mzDip);
                
                dip(1,i) = {{pxDip,pyDip,pzDip,mxDip,myDip,mzDip}};
            end
            
            xInt = obj.x0 + obj.dx *(0:obj.numX-1);
            yInt = obj.y0 + obj.dx *(0:obj.numY-1);
            zInt = obj.z0 + obj.dx *(0:obj.numZ-1);
            
            pos = {{xInt,yInt,zInt}};
        end
        
    end
end

