% 6 sections to include in an optimization setup file
% 1. Lumerical Simulation
% 2. Optimization
% 3. Frequency
% 4. Geometry Properties
% 5. Optimization Region
% 6. Figure of Merit

%% LUMERICAL SIMULATION ---------------------------------------------------------------------------
velocityMon = 'Velocity';
indexMon = 'Index';
userSim = ''; % name of .lsf to handle special user-defined co-optimization parameters
numUserSims = 0; % Number of special user-defined co-optimization parameters

%% OPTIMIZATION -----------------------------------------------------------------------------------
numIter = 10; % Number of optimization iterations to run

%% FREQEUNCY --------------------------------------------------------------------------------------
% (choose a subsequent option)

%% FREQEUNCY >> SINGLE FREQUENCY
f = 3e8/900e-9; % frequency
freq = [f f 1];

%% FREQEUNCY >> SMALL BANDWIDTH
f1 = 3e8/900e-9; % lowest of frequency range
f2 = 3e8/700e-9; % highest of frequency range (f2>=f1)
numFreq = 3; % Number of frequency data points between f1 & f2 (numFreq>=1)
freq = [f1 f2 numFreq];

%% FREQEUNCY >> LARGE BANDWIDTH (rarely used)
f1 = 3e8/900e-9; % lowest of first frequency range
f2 = 3e8/700e-9; % highest of first frequency range (f2>=f1)
f3 = 3e8/680e-9; % lowest of second frequency range
f4 = 3e8/500e-9; % highest of second frequency range (f2>=f1)
numFreq1 = 3; % Number of frequency data points between f1 & f2 (numFreq>=1)
numFreq2 = 5; % Number of frequency data points between f3 & f4 (numFreq>=1)
freq = [f1 f2 numFreq1; ...
        f3 f4 numFreq2];
% Note: In some cases, it is necessary to use the 'Large Bandwidth' setup. The software will automatically 
% run separate simulations for every sub-'freqeuncy range' for better accuracy from the simulation.

%% GEOMETRY PROPERTIES ---------------------------------------------------------------------------
% (choose a subsequent option)

%% GEOMETRY PROPERTIES >> FREEFORM
% Materials
eps_ = 'Ta2O5'; % Object permittivity or 'Material Name'
epsOut = 'SiO2 (silica) :: 300-2000nm'; % Cladding permittivity or 'Material Name'

% Constraints
dx = 5e-9; % Geometrical precision
minDimension = 300e-9; % Minimum dimension of each shape
minPadding = 0; % Minimum dimension of the gap between shapes
radiusCurv = 150e-9; % Soft radius of curvature constraint, radiusCurv >= radiusCurvHard
radiusCurvHard = 150e-9; % Hard radius of curvature constraint, radiusCurvHard <= .5*minDimension
xBC = 0; % 1 = force x-symmetry (horizontal symmetry about x = x_min+.5*x_span)
yBC = 0; % 1 = force y-symmetry (vertical symmetry about y = y_min+.5*y_span)
diagBC = 0; % 1 = force symmetry along the diagonal

% Automatically added shapes
newShapeCreation = 0; % 0 = new shapes are not allowed, 1 = allowed
newShapeRad = 150e-9; % Radius of automatically added shapes
newShapePad = newShapeRad + minPadding + dx; % Min distance b/w newly added shapes

% Optimization step size
maxMove = 40e-9; % Maximium boundary movement per iteration
maxArea = 4*(2*newShapeRad); % Maximium total area change per iteration
initialShape = 0; % (Beta Feature) massive + random-ish first iteration guess

% Import dimensions (x,y,z) and initial shape (epsGridImport) from baseFile 
[epsGridImport,x_min,x_span,y_min,y_span,z_center,z_span] = autoImportGeometry(eps_, freq, lumerical, baseFile, velocityMon, queueName, dx);

% *Optional* 
% Force 1D geometry, extruded in y (rarely used)
% z_center = y_min + .5*y_span;
% z_span = y_span;
% y_min = 0;
% y_span = dx;

% Make Geometry Object (must be called 'geo')
geo = FreeForm(x_min, y_min, z_center, x_span, y_span, dx, z_span, eps_, epsOut, xBC, yBC, diagBC, newShapeCreation, newShapeRad, newShapePad, maxMove, maxArea, minPadding, radiusCurv, radiusCurvHard, minDimension, initialShape);

% INITIAL GEOMETRY
% Options for initial geometry:
% 1. Automatic import from the baseFile.fsp
% 2. Manually define a bitmap representing the geometry
% 3. Load a geometry from a previous optimization

% *Optional*
% Set intial shape to the imported geometry from baseFile
geo = geo.setGeometry(epsGridImport, geo.xGrid, geo.yGrid);

% *Optional*
% Manually contruct a bit map to set the initial shape
% 1 = FreeForm Object (eps_), 0 = Cladding (epsOut)
% Coordinates of bitmap = geo.xGrid, geo.yGrid
% wg1 = abs(geo.yGrid-600e-9)<=150e-9 & (geo.xGrid+1900e-9)<-1100e-9;
% wg2 = abs(geo.yGrid-200e-9)<=150e-9 & (geo.xGrid+1900e-9)>1100e-9;
% wg3 = abs(geo.yGrid-400e-9)<=350e-9 & abs(geo.xGrid+1900e-9)<= 1100e-9;
% epsGridManual = wg1 + wg2 + wg3;
% geo = geo.setGeometry(epsGridManual, geo.xGrid, geo.yGrid);

% *Optional*
% Load geometry from a previous optimization
% load('folder/optResults.mat','GeoHist');
% epsGridOld = GeoHist{end}.epsGrid;
% geo = geo.setGeometry(epsGridOld, geo.xGrid, geo.yGrid);

% NON-DESIGNABLE REGION
% *Optional*
% Contruct a bit map to define non-designable regions
% 1 = designable, 0 = non-designable
wgIn = geo.xGrid<-2400e-9;
wgOut = geo.xGrid>2900e-9;
maskGrid0 = ~(wgIn|wgOut);
geo = geo.setDesignableRegion(maskGrid0, geo.xGrid, geo.yGrid);

% *Optional*
% Show preview of geometry and designable region
figure(1); imagesc(geo.xGrid(1,:),geo.yGrid(:,1),geo.epsGrid+.1*geo.maskGrid); axis equal;

%% GEOMETRY PROPERTIES >> LEVEL SET

% Materials
eps_ = 3.47^2; % Object permittivity or 'Material Name'
epsOut = 1.46^2; % Cladding permittivity or 'Material Name'

% Constraints
dx = 20e-9; % Geometrical precision
radiusCurv = 200e-9; % Soft radius of curvature constraint, radiusCurv >= radiusCurvHard
xBC = 0; % 1 = force x-symmetry (horizontal symmetry about x = x_min+.5*x_span)
yBC = 1; % 1 = force y-symmetry (vertical symmetry about y = y_min+.5*y_span)
diagBC=0; % 1 = force symetry along the diagonals

% Optimization step size
maxArea = 0.5*pi*(250e-9)^2; % Maximium total area change per iteration
maxlsIter = 300; % Maximum updates of the level set function
sideAnchoring = 1;  % Anchors the boarder points of the initial geometry so they don't move during the optimization

% Import dimensions (x,y,z) and initial shape (epsGridImport) from baseFile
[epsGridImport,x_min,x_span,y_min,y_span,z_center,z_span] = autoImportGeometry(eps_, freq, lumerical, baseFile, velocityMon, queueName, dx);

% Make Geometry Object (must be called 'geo')
geo = levelSet(x_min, y_min, z_center, x_span, y_span, dx, z_span, eps_, epsOut, xBC, yBC, maxArea, radiusCurv,diagBC,maxlsIter,sideAnchoring);

% INITIAL GEOMETRY
% Options for initial geometry:
% 1. Automatic import from the baseFile.fsp
% 2. Manually define a bitmap representing the geometry
% 3. Load a geometry from a previous optimization

% *Optional*
% Set intial shape to the imported geometry from baseFile
phi=messySD(0.5-epsGridImport,geo.dx);

% Manually contruct a bit map to set the initial shape
% 1 = LevelSet Object (eps_), 0 = Cladding (epsOut)
% Coordinates of bitmap = geo.xGrid, geo.yGrid
% epsGrid0=your_bitmap;
% phi=messySD(0.5-epsGrid0,obj.dx);

% Load geometry from a previous optimization
% load('restartphi')
% phi=restartphi;
% phi=messySD(restartphi,geo.dx);

% Set the geometry
geo = geo.setGeometry(phi);

%% GEOMETRY PROPERTIES >> INDEX MAP

% Materials
eps_ = 3^2; % Object permittivity or 'Material Name'
epsOut = 2^2; % Cladding permittivity or 'Material Name'

% Constraints
dx = 20e-9; % Geometrical precision
xBC = 0; % 1 = force x-symmetry (horizontal symmetry about x = x_min+.5*x_span)
yBC = 0; % 1 = force y-symmetry (vertical symmetry about y = y_min+.5*y_span)
diagBC = 0; % 1 = force symmetry along the diagonal

% Optimization step size
maxMove = 0.02;

% Import dimensions (x,y,z) from baseFile
[~,x_min,x_span,y_min,y_span,z_center,z_span] = autoImportGeometry(eps_, freq, lumerical, baseFile, velocityMon, queueName, dx);

% Make Geometry Object (must be called 'geo')
geo = IndexMap(x_min, y_min, z_center, x_span, y_span, dx, z_span, eps_, epsOut, xBC, yBC,diagBC,maxMove);

% Initial Geometry is automatically epsOut everywhere in space

%% OPTIMIZATION REGION

shapeType='Extruded'; % Type of geometry export to Lumerical (Do not change, default = 'Extruded');
dataList = {'GeoHist'}; % data to save every iteration {'GeoHist','EHist','VHist'} (geometry object, electric fields, gradient)
alpha = .005; % allowable decrease (will automatically change step-size if a greater decrease occurs)
testFlag = 0; % 1 = print a lot of output

% Boundary approximation parameters (Advanced, default = 1*dx)
eraseSize = 1*dx; % Distance from boundary of badly-interpolated FDTD data
velPadding = 1*dx; % Use this data to better interpolate boundary data

% Make Gradient Object (must be called 'grad')
% Recommended: Use same x_min, y_min, z_min... as the Geometry object
grad = Gradient(x_min, y_min, z_center, x_span, y_span, dx, z_span, eraseSize, velPadding);

%% FIGURE OF MERIT -----------------------------------------------------------------------------

%% SPECTRAL SPLIT of 3 frequencies

% Define a Function over the 3 main optimization parameters:
% parameters = {monitor, frequency, user-defined}
% operators = {sum, prod, minFast, min}
% left-most operator is evaluated first
mf_operatorOrder = {'user-defined','monitor','frequency'}; % order is customizable
mf_operators = {'sum','sum','minFast'};

% Create MeritFunction object
mf = MeritFunction(mf_operatorOrder, mf_operators);

% Monitors
mf_type = 'transmission';
mf_weight = [1 0 0]; % frequency dependent weights (only weight first frequency)
mf_exp = 1; % exponent
mf_dx = 10e-9; % discritization (meters) (don't make it too small)
params = [0 -1 0]; % normal vector = [nx ny nz]
mf_mon = 'merit1'; % Lumerical monitor name
merit1 = {mf_type, mf_mon, mf_weight, mf_exp, mf_dx, params};

mf_weight = [0 1 0]; % only weight second frequency
mf_mon = 'merit2'; % Lumerical monitor name
merit2 = {mf_type, mf_mon, mf_weight, mf_exp, mf_dx, params};

mf_weight = [0 0 1]; % only weight third frequency
mf_mon = 'merit3'; % Lumerical monitor name
merit3 = {mf_type, mf_mon, mf_weight, mf_exp, mf_dx, params};

% Define list of monitors: merits = [merit1; merit2; ...]
merits = [merit1; merit2; merit3];

%% TRANSMISSION of 1 frequency

% Define a Function over the 3 main optimization parameters:
% parameters = {monitor, frequency, user-defined}
% operators = {sum, prod, minFast, min}
% left-most operator is evaluated first
mf_operatorOrder = {'user-defined','frequency','monitor'}; % order is customizable
mf_operators = {'sum','sum','sum'};

% Create MeritFunction object
mf = MeritFunction(mf_operatorOrder, mf_operators);

% Monitors
mf_type = 'transmission';
mf_weight = 1; % frequency dependent weights
mf_exp = 1; % exponent
mf_dx = 10e-9; % Discritization (meters) (don't make it too small)
params = [1 0 0]; % normal vector = [nx ny nz]
mf_mon = 'merit1'; % Lumerical monitor name
merit1 = {mf_type, mf_mon, mf_weight, mf_exp, mf_dx, params};

% Define list of monitors: merits = [merit1; merit2; ...]
merits = merit1;

%% ABSORPTION in volume of 3 frequency

% Define a Function over the 3 main optimization parameters:
% parameters = {monitor, frequency, user-defined}
% operators = {sum, prod, minFast, min}
% left-most operator is evaluated first
mf_operatorOrder = {'user-defined','monitor','frequency'}; % order is customizable
mf_operators = {'sum','sum','minFast'};

% Create MeritFunction object
mf = MeritFunction(mf_operatorOrder, mf_operators);

mf_type = 'fieldenergy';
mf_weight = [1 1 1]; % frequency dependent weights
mf_exp = 1; % exponent
mf_dx = 3e-9; % discritization (meters) (don't make it too small)
params = 3; % 1 = E intensity, 2 = H intensity, 3 = absorption
mf_mon = 'merit1'; % Lumerical monitor name
merit1 = {mf_type, mf_mon, mf_weight, mf_exp, mf_dx, params};

% Define list of monitors: merits = [merit1; merit2; ...]
merits = merit1;

%% MODE MATCH at a plane

% Define a Function over the 3 main optimization parameters:
% parameters = {monitor, frequency, user-defined}
% operators = {sum, prod, minFast, min}
% left-most operator is evaluated first
mf_operatorOrder = {'user-defined','frequency','monitor'}; % order is customizable
mf_operators = {'sum','sum','sum'};

% Create MeritFunction object
mf = MeritFunction(mf_operatorOrder, mf_operators);

mf_type = 'modematch';
mf_weight = 1; % frequency dependent weights
mf_exp = 1; % exponent
mf_dx = 5e-9; % discritization (meters) (don't make it too small)
load('base2D_SiPh_waveguideCrossing_modeData.mat','modeData');
params = modeData; % modeData struct created by LumericalMethods/saveModeData.lsf
mf_mon = 'modeOut';
meritMode = {mf_type, mf_mon, mf_weight, mf_exp, mf_dx, params};

% Define list of monitors: merits = [merit1; merit2; ...]
merits = meritMode;
