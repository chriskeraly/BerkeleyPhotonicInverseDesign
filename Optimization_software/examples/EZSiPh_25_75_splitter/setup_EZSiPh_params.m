% EZ SILICON PHOTONICS
eps_ = 2.8^2; % Object permittivity or 'Material Name'
epsOut = 1.46^2; % Cladding permittivity or 'Material Name'

nb_sources_merit=2; % number of sources used to make mode

meshsize=20e-9; % mesh size you want your geometry to be defined on

radiuscurv=200e-9; % minimal radius of curvature to enforce

maxArea = 3*pi*(250e-9)^2;  % Area change at each iteration

numIter = 100; % Number of iterations to run

alpha = 0.005; % allowable decrease (will automatically change step-size if a greater decrease occurs)

%%% ADVANCED PARAMETERS

xBC = 0; % 1 = force x-symmetry (horizontal symmetry about x = x0+.5*xLengthReal)
yBC = 0; % 1 = force y-symmetry (vertical symmetry about y = y0+.5*yLengthReal)
diagBC=0; % 1 = force symmetry along the diagonals

maxlsIter=3000; % Maximum updates of the level set function
sideAnchoring=1;  % Anchors the boarder points of the initial geometry so they don't move during the optimization

numFreq=1;  % number of frequencies to optimize over
mf_weight=[1 1 1 1 1 1]; % weight of the different frequencies
