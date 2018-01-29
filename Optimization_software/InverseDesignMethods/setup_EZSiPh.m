%% LUMERICAL SIMULATION ---------------------------------------------------------------------------
setup_EZSiPh_params;
velocityMon = 'Velocity';
indexMon = 'Index';
meritMon='merit1';
runfile='runfile';

save('EZSiPh.mat');
baseFile = runfile;
[status,result] = system([lumerical ' -run LumericalMethods/EZSiPh.lsf']);
load('extracted_opt_parms');

userSim = ''; % name of .lsf to handle special user-defined co-optimization parameters
numUserSims = 0; % Number of special user-defined co-optimization parameters

%% OPTIMIZATION -----------------------------------------------------------------------------------
%numIter = 100; % Number of optimization iterations to run

%% FREQEUNCY --------------------------------------------------------------------------------------
f1 = 3e8/wavelength_stop;
f2 = 3e8/wavelength_start;
numFreq = numFreq;
freq = [f1 f2 numFreq];

%% GEOMETRY PROPERTIES >> LEVEL SET

x0 = vel_x-vel_x_span/2; % left-most x
y0 = vel_y-vel_y_span/2; % bottom-most y
z0 = vel_z; % center z
xLengthReal = vel_x_span;
yLengthReal = vel_y_span;
if (dim=='2D')
thickness = 0;
else
    thickness=vel_z_span;
end
dx = meshsize; % Geometrical precision

% Materials
eps_ = eps_; % Object permittivity or 'Material Name'
epsOut = epsOut; % Cladding permittivity or 'Material Name'

% Constraints
radiusCurv = radiuscurv; % Soft radius of curvature constraint, radiusCurv >= radiusCurvHard
xBC = xBC; % 1 = force x-symmetry (vertical symmetry about x = x0+.5*xLengthReal)
yBC = yBC; % 1 = force y-symmetry (vertical symmetry about y = y0+.5*yLengthReal)
diagBC=diagBC; % 1 = force symetry along the diagonals

% Optimization step size
maxArea = maxArea; % Maximium total area change per iteration
maxlsIter=maxlsIter; % Maximum updates of the level set function
sideAnchoring=sideAnchoring;  % Anchors the boarder points of the initial geometry so they don't move during the optimization

% Make Geometry Object (must be called 'geo')
% Recommended: Use same x0, y0, z0... as the Geometry object
geo = levelSet(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eps_, epsOut, xBC, yBC, maxArea, radiusCurv,diagBC,maxlsIter,sideAnchoring);

%%%%%%%INITIAL GEOMETRY%%%%%%%%%%%
%Either automatic import using autoImportGeometry or define your own
%epsGrid0, or your own level set initial function

%%%% Automatic import of geometry
[epsGridImport,~,~,~,~,~,~] = autoImportGeometry(eps_, freq, lumerical, baseFile, velocityMon, queueName, dx);
phi=messySD(0.5-epsGridImport,geo.dx);
%%%%

%%%% Defining from a Bitmap corresponding to geo.xGrid and geo.yGrid
%epsGrid0=your_bitmap;
%phi=messySD(0.5-epsGrid0,obj.dx);
%%%%

%%%% Defining your own level set corresponding to geo.xGird and geo.yGrid,
%%%% by hand or by import

% *Optional*
% Contruct a bit map to set the initial shape
% negative =  (eps_), positive = Cladding (epsClad)
% wg=abs(geo.yGrid)-250E-9-350E-9/2E-6*(geo.xGrid+30E-6);
% wg1=(abs(geo.yGrid)<250e-9);
% wg2=(abs(geo.xGrid)<250e-9);
% wg3=((abs(geo.xGrid)+abs(geo.yGrid))<3e-6);
% wg=wg1 + wg2+wg3;
% wg=0.5-wg;

% Load a geometry from a previous optimization
%  load('restartphi')
%  phi=restartphi;
%  phi=messySD(restartphi,geo.dx);

%%%%

geo = geo.setGeometry(phi);

%% OPTIMIZATION REGION
shapeType='Extruded'; % Type of geometry export to Lumerical (Do not change, default = 'Extruded');
dataList = {'GeoHist'};
alpha = 0.0015; % allowable decrease (will automatically change step-size if a greater decrease occurs)
testFlag = 1; % 1 = print a lot of output

% Boundary approximation parameters (Advanced, default = 1*dx)
eraseSize = 1*dx; % Distance from boundary of badly-interpolated FDTD data
velPadding = 1*dx; % Use this data to better interpolate boundary data

% Make Gradient Object (must be called 'grad')
grad = Gradient(x0, y0, z0, xLengthReal, yLengthReal, dx, thickness, eraseSize, velPadding);

%% MERIT FUNCTION

% Define a Function over the 3 main optimization parameters:
% 1 = monitor, 2 = frequency, 3 = user-defined
mf_operatorOrder = [3 2 1]; % left-most operator is evaluated first
mf_operators = {'sum','sum','sum'}; % sum, min, minFast
mf_freq = freq;
mf_numUser = numUserSims;

mf = MeritFunction(mf_operatorOrder, mf_operators);

%Monitors
mf_type = 'modematch';

load('modeinfo');

mf_weight = mf_weight;
mf_exp = 1;
mf_dx = 30e-9;
params = modeinfo; 
merit1 = {mf_type, 'merit1', mf_weight, mf_exp, mf_dx, params};

% mf_type = 'transmission';
% mf_weight = [1 1 1];
% mf_exp = 1;
% mf_dx = 10e-9;
% params = [1 0 0]; % [nx ny nz] normal vector
% merit2 = {mf_type, 'merit2', mf_weight, mf_exp, mf_dx, params};

merits = [merit1];
