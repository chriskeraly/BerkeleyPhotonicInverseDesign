classdef GeometryTemplate
    
    properties
        
    end
    
    methods
        
        function obj = GeometryTemplate()
            
        end
        
        %% GET GEOMETRY DATA
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % (x0,y0) is the bottom left corner coordinate (meters)
        % z0 is the center z coordinate (meters)
        % thickness of extruded planar geometry (meters)
        % dx is the uniform spacing of epsGrid (meters)
        % epsGrid is a binary bitmap of permittivity (1 = eps_, 0 = epsOut)
        % eps_ = relative permittivity or 'Lumerical material name'
        % epsOut = relative permittivity or 'Lumerical material name'
        function [epsGrid, eps_, epsOut, x0, y0, z0, dx, thickness] = getGeometry(obj)
            
            % Output geometry given internal properties
            
        end
        
        %% UPDATE SHAPES
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % F0 = current Merit Function
        % dFdxBnd = Gradient with respect to a boundary perturbation
        % dFdxSpace = Gradient w.r.t a pertubation anywhere in space
        % dF = predicted change to the Merit Function
        function [obj, dF] = updateShapes(obj, F0, dFdxBnd, dFdxSpace)
            
            % Update internal properties given the Gradient
            
        end

        %% CHANGE STEP-SIZE
        % *** MUST EXIST, DO NOT CHANGE RETURN ARGUMENTS ***
        % Optimization wil automatically change the step-size
        % factor < 1 for step-size reduction
        % factor >1 for step-size increase
        function obj = changeStepSize(obj, factor)
            
            % Change internal properties given the step-size change
            
        end
    end
end