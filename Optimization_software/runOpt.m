% function [opt,data] = runOpt(varargin)

close all; clearvars -except varargin;

% Load runOpt parameters
runOpt_params;

% Copy setupFile and baseFile into current directory
% if(length(varargin)<2)
    [setupFile,folderPath] = uigetfile('setup*.m','Select a setupFile.m');
    [~,setupFile,ext] = fileparts(setupFile);
    if(~strcmp(folderPath(1:end-1), pwd))
        copyfile([folderPath setupFile ext],'./');
    end
    if(strfind(setupFile,'EZSiPh'))
        setupFile='setup_EZSiPh';
    end
    [baseFile,folderPath] = uigetfile('base*.fsp','Select a baseFile.fsp');
    [~,baseFile,ext] = fileparts(baseFile);
    if(~strcmp(folderPath(1:end-1), pwd))
        copyfile([folderPath baseFile '*'],'./');
    end
% else
%     setupFile = varargin{1};
%     [folderPath,setupFile,ext] = fileparts(setupFile);
%     if(~strcmp(ext,'.m'))
%         error('setupFile must have extension .m');
%     end
%     if(~isempty(folderPath) && ~strcmp(folderPath, pwd))
%         copyfile([folderPath '/' setupFile ext],'./');
%     end
%     baseFile = varargin{2};
%     [folderPath,baseFile,ext] = fileparts(baseFile);
%     if(~strcmp(ext,'.fsp'))
%         error('baseFile must have extension .fsp');
%     end
%     if(~isempty(folderPath) && ~strcmp(folderPath, pwd))
%         copyfile([folderPath '/' baseFile ext],'./');
%     end
% end

% Initialize optimization
addpath('./examples');
addpath('./InverseDesignMethods');
addpath('./GeometryMethods');
addpath('./GeometryMethods/FreeFormMethods');
addpath('./GeometryMethods/LevelSetMethods');
addpath('./LumericalMethods');
addpath('./MeritFunctionMethods');

%multiWaitbar('Close all');
fprintf(['Setup (' setupFile ')...\n']); %multiWaitbar(['Setup (' setupFile ')'],'Busy');

eval(setupFile);
save('SetupBaseLum.mat','baseFile','queueName','velocityMon','indexMon','userSim','shapeType','merits','-v7.3');
runLumericalScript(lumerical, 'SetupBaseLum.mat', 'LumericalMethods/SetupBaseLum.lsf');
load('SetupBaseLum.mat','sim3d','xMon', 'xMonLength', 'yMon', 'yMonLength', 'zMon', 'zMonLength');
Init;
numMon = opt.mf.numMonitors; monNames = opt.mf.monNames;
save('SetupBaseLum.mat','newRun','numMon','monNames','-append','-v7.3');
notFinished = 1; first=0;
%multiWaitbar(['Setup (' setupFile ')'],'Close');
%multiWaitbar('Iteration','Color','g');

while(notFinished)
    
    % Run optimization sub-routine
    opt = opt.runOptStep;
    iter = opt.iter;
    simFlag = opt.simFlag;
    notFinished = opt.notFinished;
    %multiWaitbar('Iteration','Value',iter/opt.numIter);
    
    % Save data (optional)
    data = opt.data;
    geo = opt.geo;
    save('optResults.mat','-struct','data','-v7.3');
    pause(5);
    save('optResults.mat','-append','geo','-v7.3');
    clear data geo;
    
    % Run Lumerical simulation
    if(notFinished && (simFlag==1))
        fprintf('\nIteration %g of %g\n',iter,opt.numIter); %multiWaitbar('Iteration','Value',iter/opt.numIter);
        %multiWaitbar('Export Geometry','Reset');
        %multiWaitbar(['Forward Simulation (' queueName ')'],'Reset');
        %multiWaitbar(['Adjoint Simulation (' queueName ')'],'Reset');
        %multiWaitbar('Update Geometry','Reset');
        
        fprintf('Export Geometry...\n'); %multiWaitbar('Export Geometry','Busy');
        geoData = opt.exportGeometry;
        save('AddGeometry.mat','geoData','-v7.3');
        clear geoData;
        runLumericalScript(lumerical, 'AddGeometry.mat', 'LumericalMethods/AddGeometry.lsf');
        %multiWaitbar('Export Geometry','Reset');
        %multiWaitbar('Export Geometry','Value',1);
        
        [~, numFreqSims, numUserSims] = size(opt.forwardSolves);
        
        fprintf(['Forward Simulation (' queueName ')...\n']); %multiWaitbar(['Forward Simulation (' queueName ')'],'Busy');
        i = 1; % numMonSims always = 1 for forward simulations
        for j = 1:numFreqSims
            freqVec = opt.forwardSolves{1,j,1}.freqVec;
            freqInd = opt.forwardSolves{1,j,1}.freqInd;
            if(iter==0)
                InitMaterials;
            end
            for k = 1:numUserSims
                freqVec = opt.forwardSolves{1,j,k}.freqVec;
                runUserSim = opt.forwardSolves{1,j,k}.runUserSim;
                freqInd = opt.forwardSolves{1,j,k}.freqInd;
                numFreq = length(freqVec);
                save('runForwardSolves.mat','first','iter','j','k','freqVec','runUserSim','freqInd','numFreq','-v7.3');
                runLumericalScript(lumerical, 'runForwardSolves.mat', 'LumericalMethods/runForwardSolves.lsf');
                pause(5);load('runForwardSolves.mat');
                opt = opt.updateMeritData(merit_E,merit_H,merit_pos,merit_eps,merit_eps_pos,normParam,freqInd,k);
                opt = opt.updateGeoData(EField,pos_E,epsReal,pos_eps,1,freqInd,k,0);
                clear freqVec runUserSim freqInd numFreq;
                clear eCell hCell gridCell normParam EField pos_E epsReal pos_eps;
            end
        end
        %multiWaitbar(['Forward Simulation (' queueName ')'],'Reset');
        %multiWaitbar(['Forward Simulation (' queueName ')'],'Value',1);
        clear numFreqSims numUserSims;
        first = 0;
    elseif(notFinished && (simFlag==2))
        [numMonSims, numFreqSims, numUserSims] = size(opt.adjointSolves);
        approxFreq = opt.approxFreq;
        dipolePos = opt.mf.pos;
        
        fprintf(['Adjoint Simulation (' queueName ')...\n']); %multiWaitbar(['Adjoint Simulation (' queueName ')'],'Busy');
        for i = 1:numMonSims
            for j = 1:numFreqSims
                for k = 1:numUserSims
                    freqInd = opt.adjointSolves{i,j,k}.freqInd;
                    freqVec = opt.adjointSolves{i,j,k}.freqVec;
                    monInd = opt.adjointSolves{i,j,k}.monInd;
                    runUserSim = opt.adjointSolves{i,j,k}.runUserSim;
                    currWeights = opt.getWeights(monInd,freqInd,k);
                    currDipoles = opt.getDipoles(monInd,freqInd,k);
                    runSim = sum(abs(currWeights)) > 0;
                    numFreq = length(freqVec);
                    if(runSim)
                        save('runAdjointSolves.mat','first','iter','i','j','k','approxFreq','dipolePos','freqInd','freqVec','monInd','runUserSim','numFreq','currWeights','currDipoles','-v7.3');
                        runLumericalScript(lumerical, 'runAdjointSolves.mat', 'LumericalMethods/runAdjointSolves.lsf');
                        pause(5);load('runAdjointSolves.mat');
                        if(sourceCnt>0)
                            opt = opt.updateGeoData(EField,pos_E,epsReal,pos_eps,i,freqInd,k,1);
                        end
                    end
                    clear freqVec runUserSim freqInd numFreq monInd runUserSim currWeights currDipoles;
                    clear EField pos_E epsReal pos_eps;
                end
            end
        end
        %multiWaitbar(['Adjoint Simulation (' queueName ')'],'Reset');
        %multiWaitbar(['Adjoint Simulation (' queueName ')'],'Value',1);
        fprintf(['Update Geometry (' queueName ')...\n']); %multiWaitbar('Update Geometry','Busy');
        clear numMonSims numFreqSims numUserSims approxFreq dipolePos;
        first=0;
    end
    
end

%multiWaitbar('Close all');
%multiWaitbar([num2str(iter) 'Iterations Complete'],'Color','g');
fprintf('\nOptimization Complete\n'); %multiWaitbar([num2str(iter) 'Iterations Complete'],'Value',1);

% end
% exit;
