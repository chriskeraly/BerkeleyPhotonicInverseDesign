classdef Transmission < Monitor
    %Summary
    % Note: Changed interpolation method, should do this everywhere
    %       interpolation is done!
    
    properties
        % Parameters
        nx; ny; nz;
    end
    
    methods
        %% CONSTRUCTOR
        % Check that n is unit vector
        function obj = Transmission(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, params, testFlag)
            obj = obj@Monitor(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, testFlag);
            
            nx = params(1);
            ny = params(2);
            nz = params(3);
            nNorm = sqrt(nx^2 + ny^2 + nz^2);
            obj.nx = nx/nNorm;
            obj.ny = ny/nNorm;
            obj.nz = nz/nNorm;
            
            if(obj.dim==0)
                error('Tranmission Monitor is 0D, which is not allowed.');
            end
            
            if(obj.dim==3)
                error('Tranmission Monitor is 3D, which is not allowed.');
            end
            			
			if(obj.testFlag)
				fprintf('%%%% ADDING TRANSMISSION MONITOR %%%% \n');
				fprintf('  (x0,y0,z0) = (%g,%g,%g) \n',x0,y0,z0);
				fprintf('  (numX,numY,numZ) = (%g,%g,%g) \n',obj.numX,obj.numY,obj.numZ);
				fprintf('  (nx,ny,nz) = (%g,%g,%g) \n',obj.nx,obj.ny,obj.nz);
			end
        end
        
        %% CALCULATE MERIT FUNCTION AND DIPOLE AMPLITUDES
        % eField, hField are cells of size {numFreqx1}{1x3}
        % eDip, hDip, pos = cell(1x3)[Nx1]
        % Fix to square dipole positions
        function [F,dip,pos] = calcMerit(obj,eField,hField,grid,powerNorm,freqInd)
            
            numFreq = length(freqInd);
            
            F = zeros(1,numFreq); 
            dip = cell(1,numFreq);
          
            for i = 1:numFreq
                eFieldI = eField(i,:);
                hFieldI = hField(i,:);
                wi = obj.weight(freqInd(i));
                [eInt,hInt] = interpFields(obj,eFieldI,hFieldI,grid);
                ex = eInt{1}; ey = eInt{2}; ez = eInt{3};
                hx = hInt{1}; hy = hInt{2}; hz = hInt{3};
                Fi = 0.5*real( obj.nx*(ey.*conj(hz)-ez.*conj(hy)) ...
                    + obj.ny*(ez.*conj(hx)-ex.*conj(hz)) ...
                    + obj.nz*(ex.*conj(hy)-ey.*conj(hx)) );
                
                %figure(10); imagesc(Fi); axis equal;
                F(1,i) = wi * ( sum(sum(Fi)) / powerNorm(i) * obj.dx^obj.dim )^obj.exponent;
                eCoeff = obj.exponent * F(1,i)^(obj.exponent-1);
                
                % Check this math, prediction is closer with weird factor
                % of 1/3
                pxDip = 1/3*wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * 0.5 * obj.eps0 * ( obj.nz*conj(hy) - obj.ny*conj(hz) );
                pyDip = 1/3*wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * 0.5 * obj.eps0 * ( obj.nx*conj(hz) - obj.nz*conj(hx) );
                pzDip = 1/3*wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * 0.5 * obj.eps0 * ( obj.ny*conj(hx) - obj.nx*conj(hy) );
                mxDip = 1/3*-wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * 0.5 * obj.eps0/obj.mu0 * ( obj.ny*conj(ez) - obj.nz*conj(ey) );
                myDip = 1/3*-wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * 0.5 * obj.eps0/obj.mu0 * ( obj.nz*conj(ex) - obj.nx*conj(ez) );
                mzDip = 1/3*-wi*eCoeff*(1/powerNorm(i)) * obj.dx^obj.dim * 0.5 * obj.eps0/obj.mu0 * ( obj.nx*conj(ey) - obj.ny*conj(ex) );
                
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
