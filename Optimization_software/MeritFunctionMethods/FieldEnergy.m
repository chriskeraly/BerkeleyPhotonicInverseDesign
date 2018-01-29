classdef FieldEnergy < Monitor
    
    properties
        % Parameters
        type; % 1 = e-field intensity, 2 = h-field intensity, 3 = absorption
        freq;  
    end
    
    methods
        % Constructor
        function obj = FieldEnergy(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, params, testFlag)
            obj = obj@Monitor(weight, exponent, x0, y0, z0, xLengthReal, yLengthReal, zLengthReal, dx, testFlag);
            
            obj.type = params(1);
            if(obj.type >= 3)
                if(length(params)~=2)
                    error('FieldEnergy: must pass frequency vector as params{2}');
                end
                obj.freq = params(2);
            end
            
            if(obj.testFlag)
                fprintf('%%%% ADDING FieldEnergy MONITOR %%%% \n');
                fprintf('  (x0,y0,z0) = (%g,%g,%g) \n',x0,y0,z0);
                fprintf('  (numX,numY,numZ) = (%g,%g,%g) \n',obj.numX,obj.numY,obj.numZ);
                fprintf('  (type) = (%g) \n',obj.type);
            end
        end
        
        %% CALCULATE MERIT FUNCTION AND DIPOLE AMPLITUDES
        % eField, hField are cells of size {numFreqx1}{1x3}
        % eDip, hDip, pos = cell(1x3)[Nx1]
        % Fix to square dipole positions
        function [F,dip,pos] = calcMerit(obj,eField,hField,grid,eps_,eps_pos,powerNorm,freqInd)
            
            numFreq = length(freqInd);
            
            F = zeros(1,numFreq); 
            dip = cell(1,numFreq);
            
            for i = 1:numFreq
                eFieldI = eField(i,:);
                hFieldI = hField(i,:);
                [eInt,hInt] = interpFields(obj,eFieldI,hFieldI,grid);
                
                if( obj.type == 1 )
                    fieldInt = eInt;
                    eps_=[];
                elseif( obj.type == 2 )
                    fieldInt = hInt;
                    eps_=[];
                else % type = 3 or 4
                    fieldInt = eInt;
                    epsI = eps_(i,:);
                    [epsInt,~] = interpFields(obj,epsI,epsI,eps_pos);
                    epsInt = 1/3*( epsInt{1} + epsInt{2} + epsInt{3} );
                    epsInt(imag(epsInt)<.1)=0;
                end
                
                wi = obj.weight(freqInd(i));
                fieldX = fieldInt{1}; fieldY = fieldInt{2}; fieldZ = fieldInt{3};
                
                Fi = abs(fieldX).^2 + abs(fieldY).^2 + abs(fieldZ).^2;
                
                %figure(41); plot(Fi);
                
                if( (obj.type == 1) || (obj.type == 2) )
                    norm = 1 / (obj.numX*obj.numY*obj.numZ);
                else
                    norm = 0.5 * 8.854e-12 * 2*pi * obj.freq(i) * imag(epsInt) * obj.dx^obj.dim / powerNorm(i);
                end
                
                if(obj.dim>2)
                    FiSum = trapz(trapz(trapz(Fi .* norm)));
                elseif(obj.dim>1)
                    FiSum = trapz(trapz(Fi .* norm));
                elseif(obj.dim>0)
                    FiSum = trapz(Fi .* norm);
                else
                    FiSum = Fi .* norm;
                end
                F(1,i) = wi * FiSum^obj.exponent;
                eCoeff = obj.exponent * FiSum^(obj.exponent-1);
                
                if( (obj.type == 1) || (obj.type == 3) )
                    pxDip = wi * eCoeff * obj.eps0 * conj(fieldX) .* norm;
                    pyDip = wi * eCoeff * obj.eps0 * conj(fieldY) .* norm;
                    pzDip = wi * eCoeff * obj.eps0 * conj(fieldZ) .* norm;
                    mxDip = 0*pxDip;
                    myDip = 0*pxDip;
                    mzDip = 0*pxDip;
                elseif (obj.type == 2)
                    mxDip = -wi * eCoeff * obj.eps0/obj.mu0 * conj(fieldX) * norm;
                    myDip = -wi * eCoeff * obj.eps0/obj.mu0 * conj(fieldY) * norm;
                    mzDip = -wi * eCoeff * obj.eps0/obj.mu0 * conj(fieldZ) * norm;
                    pxDip = 0*mxDip;
                    pyDip = 0*mxDip;
                    pzDip = 0*mxDip;
                else % type = 4
                    maskGrid = imag(epsInt)~=0;
                    
                    zLen = size(maskGrid,3);
                    maskGrid2 = maskGrid;
                    maskGrid = 0*maskGrid;
                    for z=1:zLen
                        for bnd=1:8
                            maskGrid2(:,:,z) = maskGrid2(:,:,z) - getBoundaryIn(maskGrid2(:,:,z));
                            if(bnd>2)
                                maskGrid(:,:,z) = maskGrid(:,:,z) + getBoundaryIn(maskGrid2(:,:,z));
                            end
                        end
                    end
                    
                    pxDip = wi * eCoeff * obj.eps0 * conj(fieldX) .* norm .* maskGrid;
                    pyDip = wi * eCoeff * obj.eps0 * conj(fieldY) .* norm .* maskGrid;
                    pzDip = wi * eCoeff * obj.eps0 * conj(fieldZ) .* norm .* maskGrid;
                    mxDip = 0*pxDip;
                    myDip = 0*pxDip;
                    mzDip = 0*pxDip;
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
