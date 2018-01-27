%  DO NOT CHANGE: Main setup file

% Extra setup params that are not in the setup file
restartFlag = 0; % 1 = restart optimization from the following .mat file
if(restartFlag)
    load('savedOptimizationData.mat');
    restartIter = 25; % Iteration within old optimization to restart from
end
approxFreq = 0; % (Beta) 1 = Try to do all frequencies of adjoint simulation at once

if(strcmp(shapeType,'LevelSet'))
    addpath('./LevelSetMethods');
    addpath('./LSM');
end
set(0,'DefaultFigureColor',[1 1 1]);
set(0,'DefaultLineLineWidth',2);
newRun = 1;

if(restartFlag)
    opt = opt.resume(restartIter,numIter);
    iter = opt.iter;
else
    % Originally, 'merits' only contains {mf_type, 'monName', mf_weight, mf_exp, mf_dx, params}
    merits = [merits(:,1:4),num2cell(xMon.'),num2cell(yMon.'),num2cell(zMon.'),num2cell(xMonLength.'),num2cell(yMonLength.'),num2cell(zMonLength.'),merits(:,5:end)];
    validateSetupParams;
    
    for i=1:size(merits,1)
        % If monitoring absorption, then adjust weights by .5*eps0*w
        if( strcmp(merits{i,1},'fieldenergy') && (merits{i,end}(1)>=3) )
            %merits{i,end}(2) = merits{i,end}(2) * 0.5 * 8.854e-12 * 2*pi .* linspace(freq(1),freq(2),freq(3));
            merits{i,end}(2) = linspace(freq(1),freq(2),freq(3));
        end
        mf = mf.addMerit(merits(i,:));
    end
    opt = Optimizer(geo,grad,mf,freq,approxFreq,numUserSims,dataList,numIter,alpha,testFlag);

    % PlotFlags
    plot_FOMHist = 1;
    plot_FHist = 1;
    plot_Geo = 1;
    plot_dFHist = 1;
    plot_Ex = 0;
    plot_Ey = 0;
    plot_Ez = 0;
    plot_E2 = 1;
    plot_Hx = 0;
    plot_Hy = 0;
    plot_Hz = 0;
    plot_H2 = 0;
    
    opt = opt.setPlotFlags([plot_FOMHist,plot_FHist,plot_Geo,plot_dFHist]);
    clear mf geo;
end
