%Fiber based cell migration code - Multiple cell persistent random walk
%code written by Ben Yeoman
%created 10/10/2016
%code updated and edited by Tyler Collins
%last edited 2/28/2019
clear all%.........................Clears the entire workspace of variables
rng shuffle%..........................Shuffles the RNG to seed off the time
runID = num2str(randi(10000,1)-1,'%04d');%....Generates random 4 digit code
..............................................for uniquely IDing save files
close all%...............................................Closes all figures

% c = parcluster;%............................Creates local parcluster object
% c.NumWorkers = 5;%........Sets the number of workers for parallel computing
% saveProfile(c);%..........................................Saves the profile
% parpool(5);%.....................Initializes the parallel computing cluster

tic%..............................................Begins tracking real time

%Model Definitions
runTimeLimit = 6*3600;%......................Total runtime for model (sec)
dt = 2;%.......................................Size of each time step (sec)
elmtSize = 0.01;%............Size of binding site element (um) - max is 0.3
polarize = 0.2;%..............................................Polarize cell
nCells = 3;%...............................Number of cells running
oriPos = 7;%...................Distance (+/-) cells initialize from origin
% dataPnts = 1;%......Number of points in the parametric sweep of parameters
samples = 2;%........Number of repeated simulations for a given data point

filename = sprintf('MCMCluster2_%dcell_HiAI_%s_%s%s_%s',nCells,date,datestr(now,'HH'),...
    datestr(now,'MM'),runID);
% rng(datenum(date));%......................................For repeatability

fprintf('Running %d cell/s over %2.1f day/s. Start Time:  %s \n',nCells,...
    runTimeLimit/(24*3600),datestr(now))%.String announcing simulation init

%Min and max values for each parameter
% pr = [3.5 6;%.....................................Gel concentration (mg/ml)
%     0.0005 0.002;%..............................Fiber density (fibers/um^3)
%     1/16 1/2;%.........................Pseudopod extension frequency (s^-1)
%     0 0.8;%............................................Fiber alignment (AI)
%     5 8;%.........................Binding Site (RGD peptides/tropocollagen)
%     1e-9 1e-7;%..................................Adhesive binding force (N)
%     3e-9 3e-7;%...............................Compressive binding force (N)
%     0.2 1];%...........................................Friction Coeffecient

%Code Sectional Toggles 
% animation = 0;..............(1) Creates 3D Animation or (0) Creates 3D plot
gpuOn = 0;%................................Uses the GPU for matrix division
runAvgVelocity = 1;%...............Enable(1) or Disable(0) Average Velocity
runFitCurves = 1;%.......................Enable(1) or Disable(0) Fit Curves
runStatAnalysis = 1;%.....On(1)/Off(0) Stat Analysis[Avg Vel, MSD, Pers. L]

%Constant Cell Attributes
Vpseudo = 0.45;%...................Velocity of extending pseudopod (um/sec)
Lsearch = 0.3;%...........Area of searching pseudopod in contact with fiber
Fmax = 10e-9;%........................Max force provided by single cell (N)
Dcell = 15;%.............................................Cell Diameter (um)
Vcell = (4/3)*pi*(Dcell/2)^3;%.........................Cell's volume (um^3)
dRecp = 0.5*10^6;%..............................Integrin receptors per cell
tC = Lsearch/Vpseudo;%............................Binding site contact time
Bhalf = 200;%....................Number of bound ligands at half max. force
Bmin = 100;%...............Min bound ligands needed for pseudopod outgrowth
Bmax = 530;%.............Min bound ligands needed for pseudopod contraction
kcell = 1e-8;%..................................Cell spring constant (N/um)
contactRange = Dcell+5.5;%..............Contact range from cell center (um)

%Given Constants Palsson Papers
% etaSurf = 8e-8;%..............................Surface viscosity; u0(N*s/um)
etaCC = 2.5e-7;%.........Cell-cell viscosity; roughly 3*etaSurf; uC(N*s/um)
lambda = 7;%.................................................Given constant
xNot = sqrt(1/(2*lambda));%.................................Constant in Eqn
vNot = xNot*exp(-1*lambda*xNot^2);%.........................Constant in Eqn
sf = 001;%.........................................................Constant

%Constant ECM Attributes
Ltropo = 0.3;%.................................Length of tropocollagen (um)
Dtropo = 0.0015;%............................Diameter of tropocollagen (um)
Lfiber = 75;%.............................................Fiber length (um)
otAngle = 0;%.......Average angle between ECM fibers and x-axis AKA otAngle
vRef = [1 0 0];%.........................Fiber orientation referance vector
eta = 1e-11;%............................Extracellular viscosity (N*s/um^2)
kI = 0.25e-9;%........................................Bond stiffness (N/um)
koff = 5;%.............................Bond dissociation rate .1-100 (s^-1)
kon = 1.4e-3;%............................Bond association rate (um^2/s^-1)
kB = 1.38e-17;%.........................Boltzmann's contant (um^2*kg/s^2*K)
T = 311;%...................................................Temperature (K)

% % % Parameters %$%&*%)^(*&*&)(^*%*&^$&*^%(*&#)&^$)*#^)(*)^
% % CGel = 5;%........................................Gel concentration (mg/ml)
% % fiberDens = 0.0015;%.............................Fiber density [fibers/um^3] 
% % avgSearchTime = 16;%...............................Pseudopod extension time
% % FiberAI = 0;%..........................................Fiber alignment (AI)
% % RGD = 5.8;%................................RGD Peptides/monomer/tropocollagen

% Array prealocation && Output Parameter Preallocations
% MSD(:,:,samp) = msdTmp;
D = zeros(samples,nCells,1);%..........Random motility coefficient
alp = zeros(samples,nCells,1);%..............................Alpha
R1 = zeros(samples,nCells,1);%................LpR Goodness of fit
R2 = zeros(samples,nCells,1);%................MSD Goodness of fit
LpR = zeros(samples,nCells,1);%.......Persistence length by <R^2>
R3 = zeros(samples,nCells,1);%.................LpC Goodness of fit
TTC = zeros(samples,nCells,1);%..............Total Time Collective
TTS = zeros(samples,nCells,1);%..................Total Time Single
avgVel = zeros(samples,nCells,1);%............Average cell speed
celldX =  zeros(samples,nCells);%...Rebranded cellDist in stats
cMag = zeros(samples,nCells);
% celldX = zeros(samples,nCells);
steps = zeros(samples,nCells);
retr = zeros(samples,nCells);
outg = zeros(samples,nCells);
cont = zeros(samples,nCells);
folw = zeros(samples,nCells);
gpRT = zeros(runTimeLimit,nCells,samples);
avgRat = zeros(samples);
cluTime = zeros(samples);
tocTime = zeros(samples);
avgCluVel = zeros(samples);
R1C = zeros(samples);
cMagClu = zeros(samples);
celldXClu = zeros(samples);
stepsClu = zeros(samples);
LpRClu = zeros(samples);
R3C = zeros(samples);
Parameters = zeros(7,1);

% Multi-cell Preallocation - Most are unused/will be purged
% sumSurfA = zeros(nCells,1);
% gFiberDir = zeros(1,3,nCells);

% rotAx = zeros(nCells,3);
% rotAng = zeros(nCells,1);
% uvGFiberDir = zeros(nCells,3,2);
% uvPV = zeros(nCells,3);
% gLpseudo = zeros(1,nCells);
% gBonds = zeros(1,nCells);
% vLpseudo = zeros(nCells,3);

% fStep = zeros(nCells,1);
% fnPos = ones(nCells,1);

% tsurf = zeros(nCells,1);
% nSurf = zeros(nCells,1);
Fx = zeros(nCells,1);%...................X-component of a cell's net force
Fy = zeros(nCells,1);%...................Y-component of a cell's net force
Fz = zeros(nCells,1);%...................Z-component of a cell's net force
% dMMax = zeros(nCells,1);
% Lps2 = zeros(nCells,1);

% pVTest = zeros(nCells,3);
shpRat = zeros(runTimeLimit/dt,samples);%Shape Ratio aka Aspect Ratio
rWalkUn = zeros(1+runTimeLimit/dt,3,nCells,samples);%Random walk end storage
cWalkUn = zeros(1+runTimeLimit/dt,3,1,samples);%Cluster random walk end storage
MSDClu = zeros(floor((1+runTimeLimit/dt)/8),samples);%MSD for cluster
DClu = zeros(samples);
alpClu = zeros(samples);
R2C = zeros(samples);

%=========================================================================%
%Parametric Sweep
%=========================================================================%
    %Adjustable Parameters
    C_gel = 5;%........................................Gel concentration (mg/ml)
    fiberDens = 0.0015;%.............................Fiber density [fibers/um^3] 
    avgSearchTime = 16;%...............................Pseudopod extension time
    FiberAI = 0;%..........................................Fiber alignment (AI)
    RGD = 5.8;%................................RGD Peptides/monomer/tropocollagen
    Fadhv = 1600e-10;%...............................Adhesive strength constant (N) 1e-8
    Fcomp = 150e-10;%........................Compressive force during contact (N) 30e-9
    fbCoeff = 0.3;%.......................Coefficient for binding friction
    
    p = [C_gel,fiberDens,avgSearchTime,FiberAI,...
        RGD,Fadhv,Fcomp];%......Vector that contains the 
    ......................................parameters tested during sampling

    Fiber_D = sqrt((C_gel*6.022e20*(Ltropo+0.067)*Dtropo^2)/...
    (1e12*fiberDens*805000*0.7*Lfiber*0.9));%........Average fiber diameter
    dL = (RGD*.7*.9)/((Ltropo+0.067)*Dtropo);%RGD ligand density [ligand/um^2]
    sigma = sqrt(-1642*reallog(FiberAI));%....Angle deviation from orientation angle 
    if isinf(sigma)%When AI is 0 then sigma --> inf; Makes sigma usable
        sigma = 1e10;
    end
    crslnkPerL = ((4.439*sigma)/(sigma+4.55))/8;%......Crosslinks per fiber
    crslnkPerF = Lfiber*crslnkPerL;
    crslnkDens = fiberDens*crslnkPerF;%...................Crosslink Density

    %Gel stiffness as function of crosslink density - eqn from Lin et. al. 2015
    if crslnkPerF < 3.5
        gel_stiffness = 1039.9*crslnkPerF-1992.9;
    else
        gel_stiffness = 5247.9*crslnkPerF-16274;
    end%.................................................Gel stiffness (Pa)

    %Converts gel stiffness to Pa/um^2
    kecm = ones(nCells,1)*((1/crslnkDens^(1/3))*gel_stiffness/1e12);%......
    .............................................ECM spring constant (N/um)
    % kecm = ones(nCells,1)*0.5095e-05;

    %=========================================================================%
    %Sample Sweep - repeats simulation for number of given samples
    %=========================================================================%
    for samp=1:samples %USE FOR SERIAL PROCESSING
%     parfor samp=1:samples %USE FOR PARALLEL PROCESSING..... (Can be much
%     faster)
    
        rng shuffle%................Reseeds because parfor loop resets the seed
        tic%..........................................Begins tracking real time
        disp(samp)

        %=========================================================================%
        %Position and Phase Matrices
        %=========================================================================%
        cellPos = zeros(1,3,nCells);%..................Cell position matrix
        phaseRT = zeros(runTimeLimit,nCells);%....................Pseudopod phase runtime matrix
        pseudoVect = zeros(nCells,3);%...........Pseudopod direction vector
        contPos = zeros(1,3,nCells);%..........Cell contact position matrix
        cellDist = zeros(nCells,3,nCells);%.......Cell-Cell distance vector
        nCellDist = zeros(nCells);%...............Normal cell-cell distance     

        %=========================================================================%
        %Tracking variables
        %=========================================================================%
        Lpseudo = zeros(nCells,1);%..............Tracks length of pseudopod
        newAngle = ones(nCells,1);%........Tracks when cell finds new fiber
        retracting = ones(nCells,1);%...............Tracks retracting phase
        outgrowth = zeros(nCells,1);%................Tracks outgrowth phase
        contracting = zeros(nCells,1);%............Tracks contracting phase
        following = zeros(nCells,1);%................Tracks following phase
        time = zeros(nCells,1);%.......Tracks pseudopod extension frequency
        pos = 1;%....................Tracks cell position at each time step
        cellBehav = zeros(nCells,4);%.......Tracks the behavior of the cell 
        cPos = 1;%................................Step for contact position
        timeCol = zeros(nCells,1);%.........................Time Collective
%         mFollow = zeros(nCells);
        cluPos = zeros(1,3,nCells);%.......................Cluster centroid
        cluSize = zeros(1,nCells);%..............Number of cells in cluster
        noCluSize = zeros(1,nCells);%..............Number of a cluster size
%         avgRatio = zeros(1);
%         aRatioG = zeros(1,nCells);
        aRatio = zeros(1);%...............................Aspect Ratio temp
        aRat = zeros(1);%.................................Aspect Ratio temp

        %=========================================================================%
        %Temporary Variables
        %=========================================================================%
        fiberDir = zeros(nCells,3);%...................Unit vector of guiding fiber
        fiberAngle = zeros(nCells,1);%.........Angle of fiber from reference vector
        searchTime = zeros(nCells,1);%.....................Pseudopod extension time
        numElmts = zeros(nCells,1);%.....................Discrete elements of fiber
        angle = 0;%...................................................Angle tracker
        Bsites = zeros(nCells,30);%.......................Distributed binding sites
        searchLength = zeros(nCells,1);%.......Number of elements searched per step
        bonds = zeros(nCells,1);%...................................Filopodia bonds
        Bpseudo = zeros(nCells,1);%.............Number bonds along entire pseudopod
        dMove = zeros(nCells,1);%................Distance cell moves in contraction
        localFibers = zeros(nCells,1);%..............Number of fibers touching cell
        cellAx = ones(nCells,3)*Dcell/2;%..............................Axes of cell
        polAngle = 0;%...............................Cell polarity angle from fiber
        dl = zeros(nCells,1);%................................Pseudopod search area
        acAng = 0; 
        bondsTotal = zeros(nCells,1);
%         dist_byV = zeros(nCells,3);
        vPrev = zeros(nCells,3);
        otX = zeros(nCells);%............................Orientation factor
        nuCellDist = zeros(nCells,3,nCells);%.........Cell-cell unit vector
        fb = zeros(nCells,1);%................................Bond friction
        fv = zeros(nCells,1);%.............................Viscous friction
        dR = zeros(nCells,1);%.............................Receptor density
        F0 = 0;
        sl = 0;
        ve = zeros(nCells,3);%Elongation vector
        veAngle = zeros(nCells,1);
%         gVe = zeros(1,3,nCells);
        % gVeAngle = zeros(1,nCells);
        SAcell = zeros(nCells,1);%Cell surface area (um^2)
        a = ones(nCells,1)*Dcell/2;%Major axis
        b = ones(nCells,1)*Dcell/2;%Minor Axis
        c = ones(nCells,1)*Dcell/2;%Mirrored minor Axis
        rCell = zeros(nCells,1);%Cell Radius
        dmin = zeros(nCells,1);%Minimum distance threshold
        % dR(1:nCells,1) = dRecp/saellipsoid(15/2);%...Integrin Receptor density (receptors/um^2)
        Kprime = zeros(nCells,1);%Drag adjustment factor


        %New Variables
        cellContact = zeros(nCells,nCells);%Tracks if cell is within contact range
%         normCDist = zeros(1,3);
        pseudoCellAng = zeros(nCells);%.Angle between pseudopod vector and cell-cell vector
        angAvoid = zeros(nCells,nCells);%Angle range for anterior detection of cells
        dirAng = zeros(nCells);%Angle between movement and 
        tFNet = zeros(nCells,3);%Tracks FNet
        di = zeros(nCells);%Distance from cell center to cell edge towards another cell
        d = zeros(nCells);%Distance between cell edges
        surfA = zeros(nCells);%Common surface area
        ecks = zeros(nCells);%Tracks the x variable
        SE = zeros(nCells,1);%Friction of self (diagonal)
        SN = zeros(nCells);%Friction between cells
        SNSE = zeros(nCells);%SE and SN combined
        nFp = 0;%Passive force tracker
        Fp = zeros(nCells,3,nCells);%Passive force between cells
        F = zeros(nCells,1);%Active force magnitude
        nFNet = zeros(nCells,1);%FNet magnitude
        FNet = zeros(nCells,3);%FNet vector
%         FNET = zeros(1,3,nCells);%.....For graphing
%         gpd = zeros(1,3,nCells);
%         gCellAx = zeros(1,3,nCells);
%         pV2 = zeros(nCells,3);
%         gPV = zeros(1,3,nCells);
%         rWalker = zeros(1,3,nCells);
        Vinst = zeros(nCells,1);%Instantaneous velocity
        cVel = zeros(nCells,3);%Temp cell velocity
        vDMove = zeros(nCells,3);%Movement vector
        Vx = zeros(nCells,1);%Cell velocity in X
        Vy = zeros(nCells,1);%Cell velocity in Y
        Vz = zeros(nCells,1);%Cell velocity in Z
        nuVel = zeros(nCells,3);%Unit vector velocity
        bondsMax = (Bmax-1)*ones(nCells,1);%Bonds max limit
        bondsMin = zeros(nCells,1);%Bonds min tracker
        groupNo = zeros(nCells,1);%Group number for a cell
        
        %Statistical Variables
        avgVelTmp = zeros(nCells,1);
        R1Tmp = zeros(nCells,1);
        DTmp = zeros(nCells,1);
        alpTmp = zeros(nCells,1);
        R2Tmp = zeros(nCells,1);
        LpRTmp = zeros(nCells,1);
        R3Tmp = zeros(nCells,1);
        cMagTmp = zeros(nCells,1);
        celldXTmp = zeros(nCells,1);
        stepsTmp = zeros(nCells,1);
        retrTmp = zeros(nCells,1);
        outgTmp = zeros(nCells,1);
        contTmp = zeros(nCells,1);
        folwTmp = zeros(nCells,1);
        searchEnd = round(Lsearch/elmtSize)-1;
        avgCluVelTmp = zeros;
        R1CTmp = zeros;
        DCluTmp = zeros;
        alpCluTmp = zeros;
        R2CTmp = zeros;
        cMagCluTmp = zeros;
        celldXCluTmp = zeros;
        stepsCluTmp = zeros;
        LpRCluTmp = zeros;
        R3CTmp = zeros;
%         contChk = ones(nCells);
%         ratSkip = 0;
        msdTmp = zeros(floor((1+runTimeLimit/dt)/8),nCells);
        msdCluTmp = zeros(floor((1+runTimeLimit/dt)/8),1);

        %Predetermined cell init locations (Cluster shape)
        cellPosClu = [0     0     0;
                      12.31 7.5   7.5;
                      12.31 -7.5  0;
                      12.31 7.5  -7.5;
                      25   -7.5   -7.5;
                      25   7.5  0;
                      25   -7.5  7.5;
                      38 7.5   7.5;
                      38 -7.5  0;
                      38 7.5  -7.5];

        %Generates cell's initial position/polarity
        for cell=1:nCells
%             cellPos(1,:,cell) = oriPos*(2*rand(1,3)-1);%Random init
%             cellPos(1,:,cell) = [cell*15 rand/2 rand/2];%Linear init
             cellPos(1,:,cell) = cellPosClu(cell,:);%Cluster init

            %Calculates Fiber's Angle from Reference Vector #1
            if sigma == 1e10
                fiberAngle(cell) = (2*rand-1)*360;%........Randomly aligned fibers
            else
                fiberAngle(cell) = (180*(randi(2)-1))+normrnd(otAngle,sigma);%.......Aligned fibers
            end

            %Sets Initial Distance to Fiber Intersection
            Ltemp  = exprnd(1/crslnkPerL);
            if Ltemp <= Lsearch
                Ltemp = Lsearch;%...........Limited to length of pseudopod tip
            end

            %Sets Initial Polarity
            vPrev(cell,:) = 2*rand(1,3)-1;
            vPrev(cell,:) = vPrev(cell,:)/norm(vPrev(cell,:));

            %Temp Initialization
            cellAxTemp = cellAx(cell,:);
            SAcell(cell) = saellipsoid(cellAxTemp);
%             SAcell(cell) = saellipsoid(cellAx(cell,:));
%             contChk(cell,cell) = 0;
            
            %Cell-Cell Distance Vector, Magnitude, Unit Vector, and Contact Range
            for cellOne = 1:nCells
                for cellTwo = 1:nCells
                    if cellOne == cellTwo
                        continue
                    end
                    %Rij Distance Vector 
                    cellDist(cellOne,:,cellTwo) = cellPos(pos,:,cellTwo)- ...
                        cellPos(pos,:,cellOne);
                    %Rij Length
                    nCellDist(cellOne,cellTwo) = norm(cellDist(cellOne,:,cellTwo));
                    %Rij Unit Vector
                    nuCellDist(cellOne,:,cellTwo) = cellDist(cellOne,:,cellTwo)/...
                        nCellDist(cellOne,cellTwo);

                    %Establishes Cell Contact based on C-C distance
                    if nCellDist(cellOne,cellTwo) <= contactRange && nCellDist(cellOne,cellTwo) > 0
                        cellContact(cellOne,cellTwo) = 1;
        %                 timeCol(cellOne) = timeCol(cellOne) + 1;
                    else
                        cellContact(cellOne,cellTwo) = 0;
                    end
                end
            end
        end
        cluPos(1,:,1) = sum(cellPos(1,:,:),3)/nCells;%Cluster position is the centroid

        %=========================================================================%
        % DeltaT Loop - Calculates path of cell
        %=========================================================================%
        for i = 1:dt:runTimeLimit
        try
            if length(unique(groupNo)) > 1%When more than 2 groups of cells exist, break
                break
            end
            F = zeros(nCells,1);

            %Cell-Cell Distance Vector, Magnitude, Unit Vector, and Contact Range
            for cellOne = 1:nCells
                for cellTwo = 1:nCells
                    if cellOne == cellTwo
                        continue
                    end
                    %Rij Distance Vector 
                    cellDist(cellOne,:,cellTwo) = cellPos(pos,:,cellTwo)- ...
                        cellPos(pos,:,cellOne);
                    %Rij Length
                    nCellDist(cellOne,cellTwo) = norm(cellDist(cellOne,:,cellTwo));
                    %Rij Unit Vector
                    nuCellDist(cellOne,:,cellTwo) = cellDist(cellOne,:,cellTwo)/...
                        nCellDist(cellOne,cellTwo);

                    %Establishes Cell Contact based on C-C distance
                    if nCellDist(cellOne,cellTwo) <= contactRange && nCellDist(cellOne,cellTwo) > 0
                        cellContact(cellOne,cellTwo) = 1;
                    else
                        cellContact(cellOne,cellTwo) = 0;
                    end
                end
                if nnz(cellContact(cellOne,:)) > 0 
                    timeCol(cellOne) = timeCol(cellOne)+1;
                end
            end
            %Init group assignment and cluster size 
            if pos == 1
                groupNo = [1:nCells]';
                for ii = 1:nCells-1
                    for jj = ii+1:nCells
                        if ii == jj
                            continue;
                        end
                        if cellContact(ii,jj)
                            groupNo(jj) = groupNo(ii);
                        end
                    end
                end

                for jj = 1:length(unique(groupNo))
                    for ii = 1:nCells
                        if groupNo(ii) == jj
                            cluSize(pos,jj) = cluSize(pos,jj) + 1;
                        end
                    end
                end
            end
        %=========================================================================%
        % Cell Loop  
        %=========================================================================%
            for cell = 1:nCells
                % Generates new fiber direction
                if newAngle(cell) %&& not(following(cell))

                    %Number of fibers in contact with cell - must be atleast 1
                    localFibers(cell) = poissrnd(Vcell*fiberDens);
                    if localFibers(cell) < 1 || outgrowth(cell)
                        localFibers(cell) = 1;
                    end

                    %Reset fiber, angle, and length matrices
                    fiber = zeros(1,3);%.................Fiber direction matrix
                    fbf = zeros(1,3);
                    angle = zeros(1,1);%.....................Fiber angle matrix
                    acAng = zeros(1,1);%...............Fiber acute angle matrix
                    searchLength(cell) = round(Lsearch/elmtSize);%Length of fiber the cell searches for binding sites
                    L = zeros(localFibers(cell),nCells);%.........Distance to intersection
                    cs = zeros;
                    sn = zeros;

                    %Generates local fibers
                    for j=1:localFibers(cell)

                        %Calculates Distance to Next Crosslink
                        L(j,cell)  = exprnd(1/(crslnkDens)^(1/3));
                        if L(j,cell) <= Lsearch
                            L(j,cell) = Lsearch;
                        end

                        %Calculates random fiber angle
                        if sigma == 1e10
                            newDir = 2*rand(1,3)-1;%........New fiber direction
                            fiber(j,:) = newDir/norm(newDir);%..New fibers

                        %Calculates aligned fiber direction
                        else
                            fiberAngle(cell) = (180*(randi(2)-1))+normrnd(otAngle,sigma);%Aligned fibers
                            v1 = (L(j,cell)*cosd(fiberAngle(cell))*vRef)/norm(vRef);%'x' component of new vector
                            vRand = 2*rand(1,3)-1;%Random vector to get new vector rotation
                            xProd = (cross(vRef,vRand));%Gets new vector's rotation about ref vector
                            v2 =(L(j,cell)*sind(fiberAngle(cell))*xProd)/norm(xProd);%'y' component of new vector
                            newDir = (v1+v2);%............New fiber direction
                            fiber(j,:) = newDir(1,:)/norm(newDir);%Unit vector for new fibers
                        end

                        %Get fiber vectors in direction of previous direction
                        fbf(j,:) = fiber(j,:);
                        an = acosd(dot(vPrev(cell,:),fbf(j,:))/(norm(vPrev(cell,:))*norm(fbf(j,:))));
                        if an > 90
                            fbf(j,:) = -fbf(j,:);
                        end
                    end

                    %Elongation Vector
                    ve(cell,:) = sum(fbf)+vPrev(cell,:)/norm(vPrev(cell,:));

                    %Cell follows acute angle between ve and prev direction
                    veAngle(cell) = acosd(dot(vPrev(cell,:),ve(cell,:))/(norm(vPrev(cell,:))*norm(ve(cell,:))));
                    if veAngle(cell) > 90
                        ve(cell,:) = -ve(cell,:);
                    end

                    %Define axes of the cell
        %             b(cell) = ((3*Vcell*(kcell+keq))/(4*pi*((a0/b0)*keq+kcell)))^(1/3);
        %             c(cell) = b(cell);
        %             a(cell) = (3*Vcell)/(4*pi*b(cell)^2);%Major axis of cell
                    if not(following(cell)) && any(following)
                        a(cell) = 10.8791;%Major axis of cell
                        b(cell) = 6.2272;
                        c(cell) = b(cell);
                    else
                        a(cell) = 7.51;%Major axis of cell
                        b(cell) = 7.48;
                        c(cell) = b(cell);
                    end
                    cellAx(cell,:) = [a(cell),b(cell),c(cell)];%..........Cell axes
                    rCell(cell) = 7.5;%.........Radius in Palsson Eqns
                    dmin(cell) = rCell(1)*0.1*(-1);%.................Min distance from 
                    ...........cell center to non-deformed edge of neighboring cell

                    %Calculate cell sphericity
                    SAcell(cell) = saellipsoid(cellAx(cell,:));%..........Cell surface area 
                    sl = (pi^(1/3)*(6*Vcell)^(2/3))/SAcell(cell);%..Cell sphericity
                    dR(cell) = dRecp/SAcell(cell);%.......Membrane integrin density
                    polAngle = 180*sl^2;%............................Polarity angle

                    %Fiber angles with respect to previous direction
                    for j=1:localFibers(cell)
                        angle(j,1) = acosd(dot(fiber(j,:),ve(cell,:))./...
                            (norm(fiber(j,:))*norm(ve(cell,:))));
                    end

                    %Get acute angle for all fibers
                    for j=1:size(angle,1)
                        if angle(j) > 90
                            acAng(j,1) = 180-angle(j);
                        else
                            acAng(j,1) = angle(j);
                        end
                    end

                    %Determines which fiber to follow based on cell's polarity
                    cellPol = polAngle*(2*rand-1);
                    ang = abs(acAng-cellPol);%Fiber angles with respect to cell polarity
                    min_idx = ang == min(ang);
                    fiberDir(cell,:) = fiber(min_idx,:);

                    %Limits reversal of polarity if new pseudopod
                    if retracting(cell) && angle(min_idx) > 90
                        if rand >= polarize
                            fiberDir(cell,:) = -1*fiberDir(cell,:);
                        end

                    %Determines which way cell moves at crosslink
                    elseif outgrowth(cell) && angle(min_idx) > 90 
                        if rand >= 0.01
                            fiberDir(cell,:) = -1*fiberDir(cell,:);
                        end
                    end

                    %Gets distance to next crosslink
                    crossL = L(min_idx,cell);%................Distance to next crosslink
                    numElmts(cell) = floor(crossL/elmtSize);%Elements to search along fiber

                    %Distributes # of adhesion sites in each element - assumes cell touches half of fiber surface
                    n = numElmts(cell);
                    Fdm = normrnd(Fiber_D,0.02);
                    if (Fdm <= 0)
                        Fdm = 0.01;
                    end
                    dmax = max([dL,dR(cell)]);%Greater of ligand or receptor density
                    Pb = ((kon*dmax)/(kon*dmax+koff))*(1-exp(-(kon*dmax+koff)*tC));%Binding probability
                    Bsites(cell,1:n) = poissrnd(Pb*RGD*0.7*0.9*...
                        (elmtSize/(Ltropo+0.067))*(Fdm/Dtropo)*...
                        (pi/2),numElmts(cell),1);%........Distributed binding sites
                    Bpseudo(cell) = Bpseudo(cell)+...
                        sum(Bsites(cell,1:searchLength(cell)),2);%Binding sites on pseudopod

                    %Set time until new pseudopod extends
                    if retracting(cell)
                        searchTime(cell) = exprnd(avgSearchTime);%...Extension time
                        if searchTime(cell) > 2*avgSearchTime
                            searchTime(cell) = 2*avgSearchTime;%....Extension limit
                        end
                    end

                    %Reset values
                    newAngle(cell) = 0;
                    retracting(cell) = 0;
                    outgrowth(cell) = 0;
                end
            end

%                 gVe(pos,:,cell) = ve(cell,:);
            for cell = 1:nCells
                %Calculates bonds formed for each time step
                if not(contracting(cell))
                    dl(cell) = 0;
                    while(searchLength(cell) < numElmts(cell) && dl(cell) < dt*Vpseudo)%Pseudopod searches each element within dl
                        searchLength(cell) = searchLength(cell)+1;
                        dl(cell) = dl(cell)+elmtSize;%.....Increments search area by elmtSize
                        bonds(cell) = sum(Bsites(cell,(searchLength(cell)-...
                            searchEnd):searchLength(cell)),2);%....................
                        ...............................Bond density for 30 elements
                        Bpseudo(cell) = Bpseudo(cell)+Bsites(cell,searchLength(cell));

                        %Break if contracting or retracting phase
                        if bonds(cell) >= Bmax || bonds(cell) < Bmin
                            break
                        end
                    end
                end
            end
            
            %Make contracting cell remain leader while new becomes
            %follower. 
            bIdx = bonds >= Bmax;
            if any(contracting) == 0 && sum(bIdx) > 1
                if nnz(bIdx) > 1
                    fprintf('');
                    rngBonds = 0;
                    ii = 1;
                    for cellTwo = 1:nCells
                       if bonds(cellTwo)>= Bmax
                          rngBonds(ii) = cellTwo;
                          ii = ii + 1;
                       end
                    end
                    following(:) = 1;
                    following(rngBonds(randperm(length(rngBonds),1))) = 0;
                end
            elseif any(contracting) == 0 && sum(bIdx) == 1
                following(:) = 1;
                following(bIdx == 1) = 0;
            elseif any(contracting(bIdx == 1))
                following(:) = 1;
                following(contracting == 1) = 0;
            else
                following(:) = 0;
            end

            bondsMin(:) = min([bonds(:) bondsMax],[],2);
            bonds(following == 1) = bondsMin(following == 1); 

        %--------------------------------------------------------------------------
            % Angle Between Major Axis and Neighboring Cell Center & C-C Distance
            % And forcing retraction if a cell is infront of another cell
            revFiber = zeros(nCells);
            for cellOne = 1:nCells
                for cellTwo = 1:nCells
                    if cellOne == cellTwo
                        continue
                    end
                    if norm(pseudoVect(cellOne,:)) == 0
                        pseudoCellAng(cellOne,cellTwo) = atan2d(norm...
                            (cross(cellDist(cellOne,:,cellTwo),...
                            fiberDir(cellOne,:))),dot(cellDist(cellOne,...
                            :,cellTwo),fiberDir(cellOne,:)));
                    elseif bonds(cellOne) >= Bmax
                        pseudoCellAng(cellOne,cellTwo) = atan2d(norm...
                            (cross(cellDist(cellOne,:,cellTwo),...
                            pseudoVect(cellOne,:))),dot(cellDist(cellOne,...
                            :,cellTwo),pseudoVect(cellOne,:)));
                    else
                        pseudoCellAng(cellOne,cellTwo) = atan2d(norm...
                            (cross(cellDist(cellOne,:,cellTwo),...
                            pseudoVect(cellOne,:))),dot(cellDist(cellOne,...
                            :,cellTwo),pseudoVect(cellOne,:)));
                    end
                    angAvoid(cellOne,cellTwo) = atand(3.18*rCell(cellOne)/nCellDist(cellOne,cellTwo));
                    if pseudoCellAng(cellOne,cellTwo) < angAvoid(cellOne,cellTwo) &&...
                        cellContact(cellOne,cellTwo) == 1
                        if rand > 0.5
                            %Reverse polarity
                            revFiber(cellOne) = 1;
                        else
                            %Reset cell
                            revFiber(cellOne) = 2;
                        end
                    end
                end
            end
        %-------------------------------------------------------------------------
            for cell = 1:nCells
                if revFiber(cell) == 1
                    fiberDir(cell,:) = -1*fiberDir(cell,:);
                    if pseudoVect(cell,:) == -dl(cell)*fiberDir(cell,:)
                        Lpseudo(cell) = 0; 
                        time(cell) = 0;
                        bonds(cell) = 0;
                        Bpseudo(cell) = 0;
                        pseudoVect(cell,:) = zeros(1,3);
                        contracting(cell) = 0;     
                    end 
                elseif revFiber(cell) == 2
                    %Resets values
                    Lpseudo(cell) = 0;
                    time(cell) = 0;
                    bonds(cell) = 0;
                    Bpseudo(cell) = 0;
                    pseudoVect(cell,:) = zeros(1,3);
                    newAngle(cell) = 1;
                    contracting(cell) = 0;
                end

%                 if cell > 1
                if following(cell) == 1
                    bonds(cell) = min([bonds(cell) Bmax-1]); 
                end
                %=================================================================%
                % RETRACTING PHASE - pseudopod will retract if weak bonds are 
                % formed, or if new pseudopod starts extending
                %=================================================================%
        %         if bonds(cell) < Bmin || (bonds(cell) < Bmax && time(cell) > searchTime(cell))
                if bonds(cell) < Bmin || time(cell) > searchTime(cell)
                    retracting(cell) = 1;
                    phaseRT(pos,cell) = 1;

                    %Resets values
                    Lpseudo(cell) = 0;%Reset pseudopod length
                    time(cell) = 0;%Reset time searching
                    bonds(cell) = 0;%Reset bonds at tip of pseudopod
                    Bpseudo(cell) = 0;%Reset bonds along pseudopod
                    pseudoVect(cell,:) = zeros(1,3);%Reset pseudopod direction
                    newAngle(cell) = 1;%Select new fiber
                    contracting(cell) = 0;%Resets the contracting phase

                %=================================================================%
                % OUTGROWTH PHASE - pseudopod will continue growing as long as 
                % enough stable bonds are formed to allow for actin polymerization
                %=================================================================%
                elseif bonds(cell) >= Bmin && bonds(cell) < Bmax
                    outgrowth(cell) = 1;
                    phaseRT(pos,cell) = 2;
                    contracting(cell) = 0;

                    %Tracks distance and time pseudopod extends along current fiber
                    pseudoVect(cell,:) = pseudoVect(cell,:)+(dl(cell)*fiberDir(cell,:));
                    Lpseudo(cell) = norm(pseudoVect(cell,:));%..Length of pseudopod
                    time(cell) = time(cell)+dt;%...Increment time spent searching 

                    %Get new fiber if cell reaches cross-link
                    if searchLength(cell) >= numElmts(cell)
                        if rand < 0.5
                            newAngle(cell) = 1; 
                        else
                            Ltemp  = exprnd(1/(crslnkDens)^(1/3));%................
                            ....................Distance to next fiber intersection
                            if Ltemp <= Lsearch
                                Ltemp = Lsearch;%..................................
                                .................Limited to length of pseudopod tip
                            end

                            %Gets distance to next crosslink
                            crossL = L(min_idx,cell);%...Distance to next crosslink
                            numElmts(cell) = floor(crossL/elmtSize);%Number of elements to search along fiber

                            %Distributes # of adhesion sites in each element - assumes 
                            %cell touches half of fiber surface
                            n = numElmts(cell);
                            Fdm = normrnd(Fiber_D,0.02);...Randomize fiber diameter
                            if (Fdm <= 0)
                                Fdm = 0.01;
                            end
                            dmax = max([dL,dR(cell)]);%............................
                            ..................Greater of ligand or receptor density
                            Pb = ((kon*dmax)/(kon*dmax+koff))*...
                                (1-exp(-(kon*dmax+koff)*tC));%..Binding probability 
                            Bsites(cell,1:n) = poissrnd(Pb*RGD*0.7*0.9*...
                                (elmtSize/(Ltropo+0.067))*(Fdm/Dtropo)*(pi/2)...
                                ,numElmts(cell),1);%......Distributed binding sites
                            Bpseudo(cell) = Bpseudo(cell)+sum(Bsites(cell,...
                                1:searchLength(cell)),2);
                        end
                    end

                %=====================================================================%
                % CONTRACTING PHASE - pseudopod will contract when enough bonds are
                % formed to withstand the acto-myosin contractile force
                %=====================================================================%
                elseif bonds(cell) >= Bmax
                    phaseRT(pos,cell) = 3;

                    %Calculates final pseudopod direction, and magnitude
                    if not(contracting(cell))
                        contracting(cell) = 1;
                        pseudoVect(cell,:) = pseudoVect(cell,:)+(dl(cell)*...
                            fiberDir(cell,:));
                        Lpseudo(cell) = norm(pseudoVect(cell,:));
                    end
                end
                %=================================================================%
                % END Retraction, Outgrowth, and Contraction Loop
                %=================================================================%

                %Angle between major axis and cell-vector
                for cellOne = 1:nCells 
                    for cellTwo = 1:nCells
                        if cellOne == cellTwo
                            continue
                        end
                        if pos == 1
                            dirAng(cellOne,cellTwo) = atan2d(norm...
                                (cross(cellDist(cellOne,:,cellTwo),...
                                fiberDir(cellOne,:))),dot(cellDist(cellOne,...
                                :,cellTwo),fiberDir(cellOne,:)));
%                             gpd(pos,:,cellOne) = fiberDir(cellOne,:);
                        elseif bonds(cellOne) >= Bmax
                            dirAng(cellOne,cellTwo) = atan2d(norm...
                                (cross(cellDist(cellOne,:,cellTwo),...
                                pseudoVect(cellOne,:))),dot(cellDist(cellOne,...
                                :,cellTwo),pseudoVect(cellOne,:)));
%                             gpd(pos,:,cellOne) = pseudoVect(cellOne,:);                    
                        else
                            if norm(tFNet(cellOne,:)) > 1e-9
                                dirAng(cellOne,cellTwo) = atan2d(norm...
                                    (cross(cellDist(cellOne,:,cellTwo),...
                                    tFNet(cellOne,:))),dot(cellDist(cellOne,...
                                    :,cellTwo),tFNet(cellOne,:)));
%                                 gpd(pos,:,cellOne) = tFNet(cellOne,:);
                            else
%                                 gpd(pos,:,cellOne) = gpd(pos-1,:,cellOne);
                            end
                        end
                    end
                end

                %Prolate ellipsoid
                if a(cell) < b(cell)
                    beta = b(cell)/a(cell);%.....................Inverse sphericity
                    Kprime(cell) = ((4/3)*(beta^2-1))/(((2*beta^2-1)/...
                        sqrt(beta^2-1))*reallog(beta+sqrt(beta^2-1))-...
                        beta);%..............................Drag adjustment factor

                %Oblate ellipsoid
                else
                    beta = a(cell)/b(cell);
                    Kprime(cell) = ((4/3)*(beta^2-1))/((beta*(beta^2-2)/...
                        sqrt(beta^2-1))*atan(sqrt(beta^2-1))+...
                        beta);%..............................Drag adjustment factor
                end

                Ltouching = mean(cellAx(cell,:))/4;%...............................
                ...............................Estmd length of fibers touching cell
                bondsTotal(cell) = localFibers(cell)*RGD*(0.7*Ltouching/Ltropo)*...
                    (0.9*Fiber_D/Dtropo)*pi/2;%...............Bonds at trailing end
                F0 = (Fmax*Bpseudo(cell))/(Bhalf+Bpseudo(cell));%....Adhesion force

                %Effective equivalent of kBT
                if Lpseudo(cell) == 0
                    kT_eff = 0;
                else
                    kT_eff = (F0*kecm(cell)^2*Lpseudo(cell)^3)/(F0+kecm(cell)*Lpseudo(cell))^2;
                end

                if isnan(kT_eff) || isinf(kT_eff)
                    error('NaN Inf')
                end
                fb(cell) = bondsTotal(cell)*((kecm(cell)*kI)/((kecm(cell)+kI)*...
                koff));%*(exp(-(kT_eff)/(bondsTotal(cell)*kB*T)));%................
                ..................................Friction due to bond dissociation 
                fv(cell) = 6*pi*eta*a(cell)*Kprime(cell);%.........Viscous friction

                %Internal radii from cellOne to cellTwo
                for cellOne = 1:nCells 
                    for cellTwo = 1:nCells
                        if cellOne == cellTwo
                            continue
                        end
                        if contracting(cellOne) || pos == 1
                                if dirAng(cellOne,cellTwo) < 90
                                di(cellOne,cellTwo) = a(cellOne)*b(cellOne)/(sqrt((b(cellOne)*cosd(...
                                    dirAng(cellOne,cellTwo)))^2 ...
                                    +(a(cellOne)*sind(dirAng(cellOne,cellTwo)))^2));
                                else
                                    di(cellOne,cellTwo) = 6.2272;% rCell(cellOne);
                                end
                        else
        %                     di(cellOne,cellTwo) = a(cellOne)*b(cellOne)/(sqrt((b(cellOne)*cosd(...
        %                             dirAng(cellOne,cellTwo)))^2 ...
        %                             +(a(cellOne)*sind(dirAng(cellOne,cellTwo)))^2));
                            di(cellOne,cellTwo) = 7.5;
                        end
                    end
                end

                %Distance between cell membranes
                for cellOne = 1:nCells 
                    for cellTwo = 1:nCells
                        if cellOne == cellTwo
                            continue
                        end
                        d(cellOne,cellTwo) = nCellDist(cellOne,cellTwo)-...
                            di(cellOne,cellTwo)-di(cellTwo,cellOne);
                        ecks(cellOne,cellTwo) = (d(cellOne,cellTwo)-dmin(cellOne))/rCell(cellOne);
                    end
                end

                %Common Surface Area between cellOne and cellTwo
                for cellOne = 1:nCells 
                    for cellTwo = 1:nCells
                        if not(cellContact(cellOne,cellTwo)) || cellOne == cellTwo
                            surfA(cellOne,cellTwo) = 0;
                        elseif cellContact(cellOne,cellTwo)
                            surfA(cellOne,cellTwo) = 0.25*exp(-5*...
                                (ecks(cellOne,cellTwo)-dmin(cellOne))^2)*...
                                ((sf+rCell(cellOne)*((1/di(cellOne,cellTwo))+...
                                (1/di(cellTwo,cellOne))))/(2+sf));
                            surfA(cellOne,cellTwo) = surfA(cellOne,cellTwo);%.......Normalizing surface-surface
                        end
                    end
                end

        %         %SE diagonal component of SNSE and norm
        %         for cellTwo = 1:nCells        
        %             sumSurfA(cell) = sum(surfA(cell,:));
        %             if tsurf(cell) == 0
        %                 tsurf(cell) = 0;
        %             else
        %                 tsurf(cell) = sum(surfA(cell,:),2)/SAcell(cell);%Not
        %                     ...needed because we use fb and fv for SE.
        %             end
        %             if tsurf(cell) < 1
        %                 eay = 1 - tsurf(cell);
        %             else
        %                 eay = 0;
        %             end
        %         end
        %         SE(cell) = eay*etaSurf + etaCC*tsurf(cell);%Palsson version of SE

                %SE friction calculated based on contact with other cells
                if following(cell) == 1
                    SE(cell) = (fb(cell)+fv(cell))*fbCoeff;%*0.3;
                else
                    SE(cell) = (fb(cell)+fv(cell));
                end

                if SE(cell) < 0
                    error('SE is fucking up')
                end

                %SN triangular component of matrix SNSE
                for cellTwo = 1:nCells
                    if not(cellContact(cell,cellTwo)) || cell == cellTwo
                        SN(cell,cellTwo) = 0;
                    else
                        SN(cell,cellTwo) = -1*surfA(cell,cellTwo)*etaCC;
                    end
                    if SN(cell,cellTwo) > 0
                        error('SN is fucking up')
                    end
                end

                %SNSE matrix generated
                for cellOne = 1:nCells
                    for cellTwo = 1:nCells
                        SNSE(cellOne,cellTwo) = SN(cellOne,cellTwo);
                    end
                end
                for cellOne = 1:nCells
                    SNSE(cellOne,cellOne) = SE(cellOne);
                end

                %Passive Force generation between cells (from Palsson Paper)
                for cellTwo = 1:nCells
                    if cell == cellTwo
                        continue
                    end
                    otX(cell,cellTwo) = (rCell(cell)/2)*(1/di(cell,cellTwo)+...
                        1/di(cellTwo,cell));

                    if cellContact(cell,cellTwo) == 0
                        Fp(cell,:,cellTwo) = [0 0 0];
                        
                    %Compressive force
                    elseif ecks(cell,cellTwo) < 0 
                        nFp = -1*Fcomp*otX(cell,cellTwo)*...
                            (-1*ecks(cell,cellTwo))^(5/2);
        %                 nFp = -1*Fcomp*otX(cell,cellTwo)*...
        %                     (-1*ecks(cell,cellTwo))^(3/2);%....Palson version
                        nFp = min([nFp Fmax]);
                        Fp(cell,:,cellTwo) = nFp*nuCellDist(cell,:,cellTwo);
                        if nFp > 0
                            fprintf('nFp is POS and shouldnt be for compression')
                        end

                    %Adhesive force
                    else
                        nFp = -1*Fadhv*otX(cell,cellTwo)*(((ecks(cell,cellTwo)+...
                            xNot)*exp(-1*lambda*(ecks(cell,cellTwo)+xNot)^2))-...
                            (vNot*exp(-1*lambda*ecks(cell,cellTwo)^2)));
                        Fp(cell,:,cellTwo) = nFp*nuCellDist(cell,:,cellTwo);
                        if nFp < 0
                            fprintf('nFp is NEG and shouldnt be for adhesive')
                        end
                    end
                    if isnan(Fp(1,1,cell)) == 1 || isinf(Fp(1,1,cell)) == 1
                        error('Fp is NaN or Inf')
                    end
                end
                FNet(cell,:) = zeros(1,3);
                if contracting(cell)
                    %Active force generation for contracting cell
                    F(cell) = (F0*kecm(cell)*Lpseudo(cell))/...
                        (F0+kecm(cell)*Lpseudo(cell));

                    %Sum of active and passive forces
                    FNet(cell,:) = F(cell)*(pseudoVect(cell,:)/...
                        norm(pseudoVect(cell,:)))+sum(Fp(cell,:,:),3);%............
                    ...........................Force applied in pseudopod direction
                    nFNet(cell) = norm(FNet(cell,:));
                else
                    FNet(cell,:) = sum(Fp(cell,:,:),3);
                    nFNet(cell) = norm(FNet(cell,:));
                end
                if norm(FNet(cell,:)) ~= 0
                    tFNet(cell,:) = FNet(cell,:);
                end 
                

%                 FNET(pos,:,cell) = (FNet(cell,:).*1e9);%.............For graphing purposes
            end
            %=====================================================================%
            % END cell loop but still in deltaT
            %=====================================================================%
            if nnz(contracting) == 0 && nnz(cellContact) == 0%..No rxns/No movement 
                cellPos(pos+1,:,:) = cellPos(pos,:,:);
            else
                %Init GPU array (not really ever used)
                if gpuOn == 1
                    Fx = FNet(:,1);
                    Fy = FNet(:,2);
                    Fz = FNet(:,3);
                    SNSE2 = gpuArray(SNSE);
                else
                    Fx = FNet(:,1);%Splits the FNet force into X, Y, and Z
                    Fy = FNet(:,2);
                    Fz = FNet(:,3);
                end
                for x = 1:nCells
                    for y = 1:nCells
                        if isnan(SNSE(x,y)) == 1
                            error('SNSE NaN')
                        end
                    end
                end
                for x = 1:nCells
                    for y = 1:nCells
                        if isinf(SNSE(x,y)) == 1
                            error('SNSE Inf %d %d',x,y)
                        end
                    end
                end
                if nnz(FNet) == 0%If nothing has a force then V = 0
                    Vx(:) = 0;
                    Vy(:) = 0;
                    Vz(:) = 0;
                else
                    if gpuOn == 1 
                        Vx2 = SNSE2\Fx;
                        Vy2 = SNSE2\Fy;
                        Vz2 = SNSE2\Fz;
                        Vx = gather(Vx2);
                        Vy = gather(Vy2);
                        Vz = gather(Vz2);
                    else
                        %Velocity calculation with matrix division
                        Vx = mldivide(SNSE,Fx);
                        Vy = mldivide(SNSE,Fy);
                        Vz = mldivide(SNSE,Fz);
                    end
                end
                for x = 1:nCells
                    if isnan(Vx(x)) || isinf(Vx(x))
                        sprintf('isnan(Vx(x)) || isinf(Vx(x)) samp %d',samp)
                    end
                    if isnan(Vy(x)) || isinf(Vy(x))
                        sprintf('isnan(Vy(x)) || isinf(Vy(x)) samp %d',samp)
                    end
                    if isnan(Vz(x)) || isinf(Vz(x))
                        sprintf('isnan(Vz(x)) || isinf(Vz(x)) samp %d',samp)
                    end
                end

                cVel = [Vx,Vy,Vz];%..........Cell velocity vector
                vDMove = cVel*dt;%Cell movement vector
                for u = 1:nCells
                    nuVel(u,:) = cVel(u,:)/norm(cVel(u,:));%Cell velocity unit vector
%                     Vinst(u) = norm(cVel(u,:));%.........For single cell comparison
%                     vDMove(u,:) = cVel(u,:)*dt;
                end

                %Cell Movement
                for cell = 1:nCells
                    if contracting(cell)
                        %Movement relative to pseudopod
                        pMove = dot(pseudoVect(cell,:),vDMove(cell,:))/norm(pseudoVect(cell,:));%pseudoVect component of vDMove 
                        pMoveAng = acosd(dot(pseudoVect(cell,:),vDMove(cell,:))/(norm(pseudoVect(cell,:))*norm(vDMove(cell,:))));
                        if pMove == 0
                            sprintf('pMove == 0')
                        end
                        if pMoveAng > 90
                            pMoveAng = 180-pMoveAng;
                        end
                        if pMove > 0
                            pMove = min([pMove Lpseudo(cell)]);%Limits contracting cell movement if it exceeds the pseudopod
                        end
                        if abs(pMove) > 10
                           fprintf('something something darkside') 
                           if pos > 10
                               fprintf('something darker');
                           end
                        end
                        Lpseudo(cell) = Lpseudo(cell)-pMove;%Change in Lpseudo length
                        vDMove(cell,:) = (abs(pMove)/cosd(pMoveAng))*(vDMove(cell,:)/norm(vDMove(cell,:)));%vDMove adjusted for contracting cell
                        cellPos(pos+1,:,cell) = vDMove(cell,:)+...
                            cellPos(pos,:,cell);%Change in position
                    else
%                         if norm(cVel(cell,:)) > 3
%                             cVel(cell,:) = 3*nuVel(cell,:);
%                         end
                        cellPos(pos+1,:,cell) = dt*cVel(cell,:)+...
                            cellPos(pos,:,cell);%Change in position
                    end
                end
            end
            for cell = 1:nCells
                if cellPos(pos+1,:,cell) == 0
                    sprintf('cellPos = [0 0 0] for cell %d',cell) 
                end 
            end
            
            %Cluster Stats
            groupNo = [1:nCells]';
            for ii = 1:nCells-1
                for jj = ii+1:nCells
                    if ii == jj
                        continue;
                    end
                    if cellContact(ii,jj)
                        groupNo(jj) = groupNo(ii);
                    end
                end
            end
            [~,~,rnk] = unique(groupNo);
            groupNo = rnk;
            cluSize(pos+1,1:nCells) = 0;
            cluPos(pos+1,1:3,1:nCells) = 0;
%             COMMENTED OUT 11/26/2018 cluPos for multiple clusters
            for jj = 1:length(unique(groupNo))
                for ii = 1:nCells
                    if groupNo(ii) == jj
                        cluSize(pos+1,jj) = cluSize(pos+1,jj) + 1;
%                         cluPos(pos+1,:,jj) = (cluPos(pos+1,:,jj) + cellPos(pos+1,:,ii));
                    end
                end
% % %                 if cluSize(pos+1,jj) > 0 
% % %                     cluPos(pos+1,:,jj) = cluPos(pos+1,:,jj)/cluSize(pos+1,jj);
% % %                 end
            end
            cluPos(pos+1,:,1) = sum(cellPos(pos+1,:,:),3)/nCells;


            for ii = 1:nCells
                for jj = 1:nCells
                    if ii == jj
                        continue
                    end
                    %Rij Distance Vector 
                    cellDist(ii,:,jj) = cellPos(pos+1,:,jj)- ...
                        cellPos(pos+1,:,ii);
                    %Rij Length
                    nCellDist(ii,jj) = norm(cellDist(ii,:,jj));
                    %Rij Unit Vector
                    nuCellDist(ii,:,jj) = cellDist(ii,:,jj)/...
                        nCellDist(ii,jj);

                    %Establishes Cell Contact based on C-C distance
                    if nCellDist(ii,jj) <= contactRange && nCellDist(ii,jj) > 0
                        cellContact(ii,jj) = 1;
%                         timeCol(ii) = timeCol(jj) + 1;
                    else
                        cellContact(ii,jj) = 0;
                    end
                end
            end

            %Shape Calculation
            if length(unique(groupNo)) ~= 1
                aRat(pos+1,1) = [0];
            else
                shNo = ones(nCells,1);
                maxLG = zeros(nCells,2,nCells);
                if nCells == 2
                    aRat(pos+1,1) = (nCellDist(1,2)+Dcell)/Dcell;
                else
                    for ii = 1:nCells
                        for jj = 1:nCells
        %                    if cluSize(pos+1,ii) < 3
        %                        ratSkip = 1;
        %                    end
                           if ii == jj || cluSize(pos+1,ii) < 2 %Excludes comparison against self and clusters smaller than 2
                               continue
                           end
                           if groupNo(ii) == groupNo(jj) %Comparing cells within same group
                              maxLG(shNo(ii),:,ii) = [ii jj]; %Creating pairs of cells to determine the max length
                              shNo(ii) = shNo(ii) + 1;
                           end
                       end
                    end
        %             if nCells > 2
%                     aCount = 0;
                    for kk = 1:nCells
                        if nnz(maxLG(:,:,kk)) == 0
                            continue
                        end
                        maxLGTemp = unique(vertcat(maxLG(:,1,kk),maxLG(:,2,kk)));
                        maxLGTemp(1) = [];
                %     end
                        maxL = zeros(length(maxLGTemp),length(maxLGTemp));
                        for ii = 1:nCells-1
                            for jj = ii+1:nCells
                                if any(maxLGTemp(:) == ii) && any(maxLGTemp(:) == jj)
                                    maxL(ii,jj) = nCellDist(ii,jj);
                                end
                            end
                        end

                        maxLm = max(max(maxL));
                        [x,y] = find(maxL == maxLm);
                        maxVect = cellDist(x,:,y);
                        perpD = zeros(length(maxLGTemp)-2,1);
                        aRatio = 0;
                        aRat(pos+1,1) = [0];
                        count = 1;
                        for ii = 1:length(maxLGTemp)
                            if maxLGTemp(ii) == x || maxLGTemp(ii) == y || nCells == 2
                                continue
                            end

                            perpV = cellPos(pos+1,:,x) - cellPos(pos+1,:,maxLGTemp(ii));
                            perpD(count) = norm(cross(maxVect,perpV))/norm(maxVect);
                            count = count + 1;
                        end
    %                     aCount = aCount+1;
        %                 aRatio(aCount) = norm(maxVect)/(2*mean(perpD));
        %                 aRatioG(pos+1,kk) = norm(maxVect)/(2*mean(perpD));
                        if nCells == 3 && length(unique(groupNo)) > 1% length(groupNo) < 3 || 
    %                         aRatio(aCount) = (norm(maxVect)+Dcell)/Dcell;
                            aRat(pos+1,1) = (norm(maxVect)+Dcell)/Dcell;
                        else
    %                         aRatio(aCount) = (norm(maxVect)+Dcell)/(2*mean(perpD)+Dcell);
                            aRat(pos+1,1) = (norm(maxVect)+Dcell)/(2*mean(perpD)+Dcell);
                        end
                    end
                end
                if isinf(aRat(pos+1,1)) || isnan(aRat(pos+1,1)) || aRat(pos+1,1) == 0
                    fprintf('aRat is fucked')
                end
            end
%             avgRatio(pos+1,1) = sum(aRatio)/aCount;
%             end
            
%             mFollow = zeros(nCells);
        %     following = zeros(nCells,1);

            for cell = 1:nCells 
                %Should this be just pseudoVect instead of fiberDir sometimes
                if cellPos(pos+1,:,cell) ~= cellPos(pos,:,cell)
                    if norm(pseudoVect(cell,:)) == 0
                        vPrev(cell,:) = fiberDir(cell,:);
                    else 
                        vPrev(cell,:) = pseudoVect(cell,:);
                    end
                end

                if Lpseudo(cell) <= 0.05 || Lpseudo(cell) > Dcell

                    %Set time until new pseudopod extends
                    searchTime(cell) = exprnd(avgSearchTime);
                    if searchTime(cell) > 2*avgSearchTime
                        searchTime(cell) = 2*avgSearchTime;
                    end

                    Lpseudo(cell) = 0;
                    time(cell) = 0;
                    bonds(cell) = 0;
                    Bpseudo(cell) = 0;
                    pseudoVect(cell,:) = zeros(1,3);
                    contracting(cell) = 0;
                end

                if nnz(cellContact(cell,:)) > 0
                    contPos(cPos,:,cell) = cellPos(pos,:,cell);
                end
            end
            
            leader = zeros(nCells);
            following = zeros(nCells,1);

            for kk = 1:length(unique(groupNo))
                if cluSize(pos+1,kk) < 2
                    continue
                end
                for ii = 1:nCells
%                    if groupNo(ii) == kk && bonds(ii) >= Bmax
                   if groupNo(ii) == kk && contracting(ii) >= Bmax
                       leader(ii,kk) = 1;
                       break
                   end
                end
            end
            
            %if no cell in the group is contractiong then no cell in the group
            %should be following
            for kk = 1:length(unique(groupNo))
                if cluSize(pos+1,kk) < 2
                    continue
                end
        %         grpNoTmp = groupNo(:) == kk;
        %         followTmp = grpNoTmp == leader;
                for ii = 1:nCells
                   if groupNo(ii) == kk && any(leader(:,kk) == 1)
                       if leader(ii,kk) == 0
                           following(ii) = 1;
                       else
                           following(ii) = 0;
                       end
                   end
                end
            end
            
            %Commented out 8/29/18 used for tracking, needs to be made
            %available outside of parfor if you need data
% %             gLpseudo(pos,:) = Lpseudo(:);
% %             gBonds(pos,:) = bonds(:);
        
            %Commented out 8/29/18 - Used to track and graph, needs to be
            %made available outside of parfor loop if needed
            for cell = 1:nCells
% % % %                 gFiberDir(pos,:,cell) = fiberDir(cell,:);
%                 gCellAx(pos,:,cell) = cellAx(cell,:);
%                 gPV(pos,:,cell) = pseudoVect(cell,:);
            end

            if nnz(cellContact) ~= 0
                cPos = cPos+1;
            end
            pos = pos+1;

            for x = 1:nCells
                if isnan(cellPos(pos,:,x))
                    fprintf('NaN cell position')
                end
            end
        catch
           fprintf('hidy ho') 
        end
        end
        %=========================================================================%
        % END deltaT loop
        %=========================================================================%
        try
        runTime = (pos-1)*2;
        % Initiate matrices for MSD calculation
        t = ((0:1+runTime/dt-1)');
        tMax = (floor((1+runTime/dt)/8));
        tMsd = dt*((1:tMax)')/3600;
        totTMax = (floor((1+runTimeLimit/dt)/8));
        totTMSD = dt*((1:totTMax)')/3600;
%         
%         rWalk = zeros(1+runTime/dt,3,nCells,samples);
%         cWalk = zeros(1+runTime/dt,3,nCells,samples);
        rWalk = zeros(1+runTime/dt,3,nCells);
%         cWalk = zeros(1+runTime/dt,3,1);
        
        for j = 1:size(rWalk,1)
            noCluSize(j,:) = 0;
            for k = 1:nCells
                for i = 1:nCells
                  if cluSize(j,k) == i
                      noCluSize(j,i) = noCluSize(j,i) + 1;
                  elseif nnz(noCluSize(j,i)) == 0 
                      noCluSize(j,i) = 0;
                  end
                end
            end
        end
        cluSizeTmp = cluSize;
        noCluSizeTmp = noCluSize;
        sumNoClu = sum(noCluSize,2);
        avgNoCluTmp = sum(sumNoClu)/size(cellPos,1);
        sumCluSize = sum(noCluSize);
        pctCluSize = sumCluSize/sum(sumCluSize);
        
        rWalk = cellPos;
        cWalk = cluPos(:,:,1);
%         p1(:,samp,:) = cellPos(1,:,:);
%         p1(:,samp,cell) = p1Tmp(:,cell);
%         p2(:,samp,:) = cellPos(end,:,:);%Used for graphing, might not use
        %=================================================================%
        % Collective and Individual Time
        %=================================================================%
        TTC(samp,:) = (timeCol*2)/runTime;
        TTS(samp,:) = 1-(timeCol*2)/runTime;
        catch 
            fprintf('inb4 you guessed')
        end
        %=========================================================================%
        % Statistics Section
        %=========================================================================%
        if runStatAnalysis == 1
            
            if pos > 10
                try
                for i = 1:nCells
                    rw = cellPos(:,:,i);

                    %=================================================================%
                    % Cell Velocity
                    %=================================================================%
                    if runAvgVelocity == 1
                        pos = 1;
                        w = 1;
                        tm = (1:runTime/10:runTime)./3600;%................................10 time points to get distance
                        dst = zeros;
                        for q=1:(size(rw,1))/10:size(rw,1)
                            Mag = zeros(size(rw,1),1);
                            for j=1:q
                                if rw(j+1,:) ~= rw(j,:)
                                    Mag(pos,1) = norm(rw(j+1,:)-rw(j,:));%.......Gets distances between each direction change
                                    pos = pos+1;
                                end
                            end
                            dst(w,1) = sum(Mag(:,1));%...............................Gets distance traveled at each time point
                            w = w+1;
                        end

                        if runFitCurves == 1
                            try %..................................................Curve fit to get slope for average velocity
                                [f,g] = fit(tm',dst(:,1),'poly1');
                                c = coeffvalues(f);
                                avgVelTmp(i) = c(1);%......................................Cell Speed
                                R1Tmp(i) = g.rsquare;%.....................................Goodness of fit
                            catch
                                disp('Velocity Fitting Error')
                                avgVelTmp(i) = 0;
                                R1Tmp(i) = 0;
                            end
                        end
                    end

                    %=====================================================================%
                    % Non-Overlapping MSD
                    %=====================================================================%
                    msd = zeros;
                    m = 1;
                    for tau = 1:tMax

                        %Gets nonoverlapping points for tau
                        idx = mod(t,tau+1) == 0;
                        rwIdx = rw(idx,:);

                        dspl = rwIdx(2:end,1:3)-rwIdx(1:end-1,1:3); %xyz displacement matrix
                        D_squared = sum(dspl.^2,2); %squared displacement
                        msd(m,1) = mean(D_squared); %MSD
                        m = m+1;
                    end
                    
    % %                 MSD(:,i,samp) = msd;
                    msdTmp(1:tMax,i) = msd;

                    if runFitCurves == 1
                        try
                            [f,g] = fit(tMsd,msd,'power1');
                            c = coeffvalues(f);
                            DTmp(i) = c(1);
                            alpTmp(i) = c(2);
                            R2Tmp(i) = g.rsquare;
                        catch
                            disp('MSD Fitting Error')
                            DTmp(i) = 0;
                            alpTmp(i) = 0;
                            R2Tmp(i) = 0;
                        end
                    end

                    %===================================================================================================================================%
                    %Persistence Length 
                    %===================================================================================================================================%
                    %Creates indices and variables
                    nPos = zeros;
                    nVect = zeros;
                    nMag = zeros;
%                     avgCos = zeros;
                    avgRsquared = zeros;
                    contourLength = zeros;
                    pos = 1;

                    %Gets direction vector and magnitude for each change in position
                    for j=1:size(rw,1)-1
                        if rw(j+1,:) ~= rw(j,:)
                            nPos(pos,1:3) = rw(j,:);%Position matrix for each new position
                            nVect(pos,1:3) = rw(j+1,:)-rw(j,:);%Vector matrix for each direction change
                            nMag(pos,1) = norm(rw(j+1,:)-rw(j,:));%Magnitude matrix for each vector
                            pos = pos+1;
                        end
                    end
                    nPos(pos,1:3) = rw(size(cellPos,1),1:3);
                    cMagTmp(i) = mean(nMag);
                    celldXTmp(i) = sum(nMag);
                    stepsTmp(i) = runTime/dt;

    % %                 xl = 60; %floor(sum(n_mag)/8);%.....Contour length limit for Lp plots
                    xl = 10;
                    n = 0;%.........Tracks steps for L0

                    %Loop for increasing L0
                    for L0=0:(xl/300):xl

                        n = n+1;%............................Increments each L0 step
                        endPos = 1;%..............First direction vector
                        m = 0;%..............Tracks cos(theta)/R for each L0 step
                        stop = 0;%..............Breaks while loop if pos equals size of vector matrix
                        totalCos = 0;
                        totalR = 0;

                        %Calculates average cos(theta)/R^2 for given L0
                        while endPos < size(nVect,1)

                            Length = 0;%.........Resets L after moving step length of L0
                            m = m+1;%.........Increments steps for cos(theta)/R
                            strtPos = endPos;%...Gets starting position of L0 contour

                            %Gets the contour position at end of length L0
                            while Length <= L0
                                if endPos == size(nVect,1)%Breaks loop at end of contour length
                                    stop = 1;
                                    break;
                                end

                                Length = Length + nMag(endPos);%Calculates length of contour
                                endPos = endPos + 1;%Gets end position of L0 contour
                            end

                            if stop == 0
            %                     %cos(theta) between vectors at start and end of contour
            %                     cos_theta = dot(n_vect(end_pos-1,:),n_vect(strt_pos,:))/(n_mag(end_pos-1)*n_mag(strt_pos));
            %                     total_cos = total_cos + cos_theta;

                                %Distance R between start and end points of contour
                                Rsquared = (norm(nPos(endPos,:)-nPos(strtPos,:)))^2;
                                totalR = totalR + Rsquared;
                            end
                        end

    %                     if totalR == 0
    %                         Rsquared = (norm(nPos(endPos,:)-nPos(strtPos,:)))^2;
    %                         totalR = totalR + Rsquared;
    %                     end
            %             avg_cos_step = total_cos/m;
                        avgRstep = totalR/m;

            %             avg_cos(n,1) = avg_cos_step;                                  %Average cos(theta) for given contour length
                        avgRsquared(n,1) = avgRstep;                             %Average R squared for given contour

                        contourLength(n,1) = L0;


                    end
                    if runFitCurves == 1
                        ft = fittype('2*(b1^2)*(exp(-(x/b1))-1+x/b1)');
                        try %...............................................Curve fit to get persistence length
                            [f,g] = fit(contourLength,avgRsquared,ft,'StartPoint',1);
                            c = coeffvalues(f);
                            LpRTmp(i) = c(1);%.........................................Persistence length
                            R3Tmp(i) = g.rsquare;%.....................................Goodness of fit
                        catch
                            LpRTmp(i) = 0;
                            R3Tmp(i) = 0;
                        end
                    end

                    retrTmp(i) = sum(phaseRT(:,i) == 1);
                    outgTmp(i) = sum(phaseRT(:,i) == 2);
                    contTmp(i) = sum(phaseRT(:,i) == 3);
                    folwTmp(i) = sum(phaseRT(:,i) == 4);

            %         con_L(:,i) = contour_length;
            %         aR_sq(:,i) = avg_R_squared;
            %         av_cs(:,i) = avg_cos;
            % 
            %         dist1(i) = sum(n_mag);
            %         avgVel(i) = sum(n_mag)/runTime;
            % 
            %             con_L(:,i,sample) = contour_length;
            %             aR_sq(:,i,sample) = avg_R_squared;
            %             av_cs(:,i,sample) = avg_cos;
            %         
            %             dist1(i,sample) = sum(nMag);
            %             avgVel(i,sample) = sum(nMag)/run_time;
            %         

                end
        catch
            fprintf('single stats suck')
%             error('')
        end
    %=========================================================================
    %Cluster Stats
    %=========================================================================
        try
                cw = cluPos(:,:,1);
    %=========================================================================            
                %Cluster Velocity
    %=========================================================================           
                pos = 1;
                w = 1;
                tm = (1:runTime/10:runTime)./3600;%10 time points to get distance
                dst = zeros;
                for q=1:(size(cw,1))/10:size(cw,1)
                    Mag = zeros(size(cw,1),1);
                    for j=1:q
                            if cw(j+1,:) ~= cw(j,:)
                                Mag(pos,1) = norm(cw(j+1,:)-cw(j,:));%.......Gets distances between each direction change
                            end
                        pos = pos+1;
                    end
                    dst(w,1) = sum(Mag(:,1));%...............................Gets distance traveled at each time point
                    w = w+1;
                end

                if runFitCurves == 1
                    try %..................................................Curve fit to get slope for average velocity
                        [f,g] = fit(tm',dst(:,1),'poly1');
                        c = coeffvalues(f);
                        avgCluVelTmp = c(1);%......................................Cell Speed
                        R1CTmp = g.rsquare;%.....................................Goodness of fit
                    catch
                        disp('Velocity Fitting Error')
                        avgCluVelTmp = 0;
                        R1CTmp = 0;
                    end
                end

    %=========================================================================            
                %Cluster Non-Overlapping MSD
    %=========================================================================            
                msdClu = zeros;
                m = 1;
                for tau = 1:tMax

                    %Gets nonoverlapping points for tau
                    idx = mod(t,tau+1) == 0;
                    cwIdx = cw(idx,:);

                    dspl = cwIdx(2:end,1:3)-cwIdx(1:end-1,1:3); %xyz displacement matrix
                    D_squared = sum(dspl.^2,2); %squared displacement
                    msdClu(m,1) = mean(D_squared); %MSD
                    if m > 1
                        if mean(D_squared) - msdClu(m-1,1) > 100
                            sprintf('MSD is being a bitch')
                        end
                    end
                    m = m+1;
                end
                msdCluTmp(1:tMax,1) = msdClu;
                if runFitCurves == 1
                    try
                        [f,g] = fit(tMsd,msdClu,'power1');
                        c = coeffvalues(f);
                        DCluTmp = c(1);
                        alpCluTmp = c(2);
                        R2CTmp = g.rsquare;
                    catch
                        disp('MSD Fitting Error')
                        DCluTmp = 0;
                        alpCluTmp = 0;
                        R2CTmp = 0;
                    end
                end

    %=========================================================================
                %Cluster Persistence Length
    %=========================================================================            

                %Creates indices and variables
                nPos = zeros;
                nVect = zeros;
                nMag = zeros;
%                 avgCos = zeros;
                avgRsquared = zeros;
                contourLength = zeros;
                pos = 1;

                %Gets direction vector and magnitude for each change in position
                for j=1:size(cw,1)-1
                    if cw(j+1,:) ~= cw(j,:)
                        nPos(pos,1:3) = cw(j,:);%Position matrix for each new position
                        nVect(pos,1:3) = cw(j+1,:)-cw(j,:);%Vector matrix for each direction change
                        nMag(pos,1) = norm(cw(j+1,:)-cw(j,:));%Magnitude matrix for each vector
                        pos = pos+1;
                    end
                end
                nPos(pos,1:3) = cluPos(size(cluPos,1),1:3);
    %             nPos(pos,1:3) = cw(size(cluPos,1),1:3);
                cMagCluTmp = mean(nMag);
                celldXCluTmp = sum(nMag);
                stepsCluTmp = runTime/dt;

                xl = 10; %floor(sum(n_mag)/8);%.....Contour length limit for Lp plots
                n = 0;%.........Tracks steps for L0

                %Loop for increasing L0
                for L0=0:(xl/300):xl

                    n = n+1;%............................Increments each L0 step
                    endPos = 1;%..............First direction vector
                    m = 0;%..............Tracks cos(theta)/R for each L0 step
                    stop = 0;%..............Breaks while loop if pos equals size of vector matrix
                    totalCos = 0;
                    totalR = 0;

                    %Calculates average cos(theta)/R^2 for given L0
                    while endPos < size(nVect,1)

                        Length = 0;%.........Resets L after moving step length of L0
                        m = m+1;%.........Increments steps for cos(theta)/R
                        strtPos = endPos;%...Gets starting position of L0 contour

                        %Gets the contour position at end of length L0
                        while Length <= L0
                            if endPos == size(nVect,1)%Breaks loop at end of contour length
                                stop = 1;
                                break;
                            end

                            Length = Length + nMag(endPos);%Calculates length of contour
                            endPos = endPos + 1;%Gets end position of L0 contour
                        end

                        if stop == 0
        %                     %cos(theta) between vectors at start and end of contour
        %                     cos_theta = dot(n_vect(end_pos-1,:),n_vect(strt_pos,:))/(n_mag(end_pos-1)*n_mag(strt_pos));
        %                     total_cos = total_cos + cos_theta;

                            %Distance R between start and end points of contour
                            Rsquared = (norm(nPos(endPos,:)-nPos(strtPos,:)))^2;
                            totalR = totalR + Rsquared;
                        end
                    end

    % % %                 if totalR == 0%CHEESE
    % % %                     Rsquared = (norm(nPos(endPos,:)-nPos(strtPos,:)))^2;
    % % %                     totalR = totalR + Rsquared;
    % % %                 end

        %             avg_cos_step = total_cos/m;
                    avgRstep = totalR/m;

        %             avg_cos(n,1) = avg_cos_step;                                  %Average cos(theta) for given contour length
                    avgRsquared(n,1) = avgRstep;                             %Average R squared for given contour

                    contourLength(n,1) = L0;
                end
                if runFitCurves == 1
                    ft = fittype('2*(b1^2)*(exp(-(x/b1))-1+x/b1)');
                    try %...............................................Curve fit to get persistence length
                        [f,g] = fit(contourLength,avgRsquared,ft,'StartPoint',1);
                        c = coeffvalues(f);
                        LpRCluTmp = c(1);%.........................................Persistence length
                        R3CTmp = g.rsquare;%.....................................Goodness of fit
                    catch
                        LpRCluTmp = 0;
                        R3CTmp = 0;
                    end
                end
        catch
            fprintf('cluster stats sucks')
%             error('')
        end
            else
                %Single
                avgVelTmp(:) = 0;
                R1Tmp(:) = 0;
                DTmp(:) = 0;
                alpTmp(:) = 0;
                R2Tmp(:) = 0;
                LpRTmp(:) = 0;
                R3Tmp(:) = 0;
                            
                %Cluster
                avgCluVelTmp = 0;
                R1CTmp = 0;
                DCluTmp = 0;
                alpCluTmp = 0;
                R2CTmp = 0;
                LpRCluTmp = 0;
                R3CTmp = 0;
            end
try 
%         error('Error Chk')

%                 cMagClu(samp,:) = cMagCluTmp;
%                 celldXClu(samp,:) = celldXCluTmp;
%                 stepsClu(samp,:) = stepsCluTmp;
%                 LpRClu(samp,:) = LpRCluTmp;
%                 R3C(samp,:) = R3CTmp;
%             ndZero = 2+((runTimeLimit-runTime)/dt);
%             ndZero = (runTimeLimit/dt)-(runTime/dt);
            aRat(1) = [];
            ndZero = (runTimeLimit-runTime)/dt;
            if ndZero > 0 
%                 continue
%             else
%                 for cell = 1:nCells
    %                 rWalk(end:1+runTimeLimit/dt,:,cell) = [0,0,0];
    %                 FNET(end:1+runTimeLimit/dt,:,cell) = [0,0,0];
    %                 gpd(end:1+runTimeLimit/dt,:,cell) = [0,0,0];
    %                 gCellAx(end:1+runTimeLimit/dt,:,cell) = [0,0,0];
    %                 gPV(end:1+runTimeLimit/dt,:,cell) = [0,0,0];
                    rWalk(end+1:1+runTimeLimit/dt,:,:) = zeros(ndZero,3,nCells);
                    cWalk(end+1:1+runTimeLimit/dt,:,1) = zeros(ndZero,3,1);
%                     FNET(end+1:runTimeLimit/dt,:,:) = zeros(ndZero,3,nCells);
%                     gpd(end+1:runTimeLimit/dt,:,:) = zeros(ndZero,3,nCells);
%                     gCellAx(end+1:runTimeLimit/dt,:,:) = zeros(ndZero,3,nCells);
%                     gPV(end+1:runTimeLimit/dt,:,:) = zeros(ndZero,3,nCells);
%                     phaseRT(end+1:runTimeLimit/dt,:) = zeros(ndZero,nCells);
%                 end
                %Usable for only Cluster
%                 if ratSkip == 1
%                     aRat = zeros(1+runTimeLimit/dt,1);
%                 elseif nCells == 2
%                     aRat = 2*ones(1+runTimeLimit/dt,1);
%                 else
                    aRat(end+1:runTimeLimit/dt,1) = zeros(ndZero,1);    
%                 end
            end
            
%             error('Error Chk')
            try
            rWalkUn(:,:,:,samp) = rWalk;
            cWalkUn(:,:,1,samp) = cWalk;
            catch
                fprintf('1')
            end
            try
%             aRat(end:1+runTimeLimit/dt) = 0;
            shpRat(:,samp) = aRat;
            avgRat(samp) = sum(aRat,1)/sum(aRat~=0,1);
            cluTime(samp) = runTime;
            tocTime(samp) = toc;
            catch
                fprintf('2')
            end
            try
            %Temporary Var unload to Par and Samp Var
            avgVel(samp,:) = avgVelTmp;
            R1(samp,:) = R1Tmp;
            avgCluVel(samp) = avgCluVelTmp;
            R1C(samp) = R1CTmp;
            catch
                fprintf('3')
            end
            try 
%             MSD(:,:,samp) = msdTmp;
            D(samp,:) = DTmp;
            alp(samp,:) = alpTmp;
            R2(samp,:) = R2Tmp;
%             MSDClu(:,samp) = msdCluTmp;
            DClu(samp) = DCluTmp;
            alpClu(samp) = alpCluTmp;
            R2C(samp) = R2CTmp;
            catch
                fprintf('4')
            end
            try
            cMag(samp,:) = cMagTmp;
            celldX(samp,:) = celldXTmp;
            steps(samp,:) = stepsTmp;
            LpR(samp,:) = LpRTmp;
            R3(samp,:) = R3Tmp;
            cMagClu(samp) = cMagCluTmp;
            celldXClu(samp) = celldXCluTmp;
            stepsClu(samp) = stepsCluTmp;
            LpRClu(samp) = LpRCluTmp;
            R3C(samp) = R3CTmp;
            catch
                fprintf('5')
            end
            try
            retr(samp,:) = retrTmp;
            outg(samp,:) = outgTmp;
            cont(samp,:) = contTmp;
            folw(samp,:) = folwTmp;
            catch
                fprintf('6')
            end
            
            %For graphing
%             FNETUn(:,:,:,samp) = FNET;
%             gpdUn(:,:,:,samp) = gpd;
%             gCellAxUn(:,:,:,samp) = gCellAx;
%             gPVUn(:,:,:,samp) = gPV;
            try
            gpRT(:,:,samp) = phaseRT;
            catch
                fprintf('7')
            end
    %         pV2 = zeros(nCells,3);
            
%                     error('Error Chk')
%             cellPosUn(:,:,:,samp) = 
catch
    fprintf('end of the line')
end
        end
    end
    Parameters(:) = [C_gel,fiberDens,avgSearchTime,FiberAI,RGD,Fadhv,Fcomp];
    tclock = [floor(toc/3600), mod(floor(toc/60),60), mod(toc,60)];
    disp(tclock)
% end
% Stores all cell migration data
% ===================================================================================================================================%
try
totSamp = samples*nCells;
avgVelRS = reshape(avgVel,[totSamp 1]);
LpRRS = reshape(LpR,[totSamp 1]);
DRS = reshape(D,[totSamp 1]);
alpRS = reshape(alp,[totSamp 1]);
retr1 = reshape(retr,[totSamp 1]);
outg1 = reshape(outg,[totSamp 1]);
cont1 = reshape(cont,[totSamp 1]);
folw1 = reshape(folw,[totSamp 1]);
R1r = reshape(R1,[totSamp 1]);
R2r = reshape(R2,[totSamp 1]);
R3r = reshape(R3,[totSamp 1]);
            
v1 = avgVelRS(1:totSamp);
% v2 = avgVelRS(:,totSamp/2+1:totSamp);

p1 = LpRRS(1:totSamp);
% p2 = LpRRS(:,totSamp/2+1:totSamp);

d1 = DRS(1:totSamp);
% d2 = DRS(:,totSamp/2+1:totSamp);

a1 = alpRS(1:totSamp);
% a2 = alpRS(:,totSamp/2+1:totSamp);
endTime = toc;
% save(filename,'avgVel','celldX','LpR','D','alp','R1','R2','R3','R1r','R2r','R3r','outg','cont','retr','folw','outg1','cont1','retr1','folw1','p1','TTC','TTS','v1','p1','d1','a1','totSamp','endTime','shpRat')

save(filename);
catch
   fprintf('save game failed') 
end
% tclock = [floor(toc/3600), mod(floor(toc/60),60), mod(toc,60)];
% Commented out 8/28/18 to deal with less shizzle
% % head = gobjects(nCells,1);
% % vHead = gobjects(nCells,1);
% % curve = gobjects(nCells,1);
% % 
% % for cell = 1:nCells
% %     curve(cell) = animatedline('Color',abs(rand(3,1)),'LineWidth',2);
% % end
% % 
% % axis([-50 50 -50 50 -50 50])
% % % axis([-300 300 -300 300 -300 300])
% % 
% % grid on
% % daspect([1 1 1])
% % xlabel('x[\mum]')
% % ylabel('y[\mum]')
% % zlabel('z[\mum]')
% % hold on
% % orng  = [255 175 0]/255;
% % rw = cellPos;
% % if animation == 1
% %     for j=1:pos-1 %3D Animation plot
% %         for cell = 1:nCells
% %             if phaseRT(j,cell) == 1
% %                 C = [1 0 0];%...........................................Red
% %             elseif phaseRT(j,cell) == 2
% %                 C = [1 1 0];%........................................Yellow
% %             elseif phaseRT(j,cell) == 3
% %                 C = [0 1 0];%.........................................Green
% %             else
% %                 C = [.61 .51 .74];%..................................Purple
% % %                 C = [0 0 0];%.........................................Black
% %             end
% %             
% %             [x,y,z] = ellipsoid(0,0,0,gCellAx(j,3,cell),gCellAx(j,2,cell),gCellAx(j,1,cell),10);
% %             addpoints(curve(cell),rw(j,1,cell),rw(j,2,cell),rw(j,3,cell))
% %             head(cell) = surf(x+rw(j,1,cell),y+rw(j,2,cell),z+rw(j,3,cell),'FaceColor',C);
% % 
% % %             if phaseRT(j,cell) == 3 || j == 1
% % %             if j == 1
% % %             if norm(gPV(j,:,cell)) ~= 0
% %                 pV2(cell,:) = gpd(j,:,cell);
% % %                 pV2(cell,:) = gPV(j,:,cell);gpd
% % %             else
% % %                 pV2(cell,:) = gPV(j,:,cell);
% % %                 pV2(cell,:) = gFiberDir(j,:,cell);
% % %                 pV2(cell,:) = gVe(j,:,cell);
% % %             end
% % %             end
% %             rotate(head(cell),[0 0 1],atand(pV2(cell,2)/pV2(cell,1)),rw(j,:,cell))
% %             rotate(head(cell),cross([0 0 1],pV2(cell,:)),acosd(pV2(cell,3)/norm(pV2(cell,:))),rw(j,:,cell))
% %             alpha(head(cell),0.25);
% %             
% %             if FNET(j,:,cell) ~= rw(j,:,cell)
% %                 vHead(cell) = quiver3(rw(j,1,cell),rw(j,2,cell),rw(j,3,cell),FNET(j,1,cell),FNET(j,2,cell),FNET(j,3,cell),'LineWidth',3);
% %             end
% %         end
% %         drawnow
% %         delete(vHead(:));
% %         delete(head(:))
% %         
% %         if j>1
% %             for r=0:log10(j-1)
% %                 fprintf('\b'); % delete previous counter display
% %             end
% %         end
% %         
% % %         fprintf('\b')
% %         fprintf('%d',j);
% %     end
% % else
% %     for i=1:nCells % 3D scatter plot
% %         scatter3(rw(1,1,i),rw(1,2,i),rw(1,3,i),25,'filled','MarkerEdgeColor','g','MarkerFaceColor','g','LineWidth',2)
% %         plot3(rw(:,1,1),rw(:,2,1),rw(:,3,1),'Color','m')
% %         if nCells > 1 
% %             plot3(rw(:,1,2),rw(:,2,2),rw(:,3,2),'Color','b')
% %         end
% %         if i > 2
% %             plot3(rw(:,1,i),rw(:,2,i),rw(:,3,i),'Color','k')
% %         end
% %         scatter3(rw(end,1,i),rw(end,2,i),rw(end,3,i),25,'MarkerEdgeColor','r','MarkerFaceColor','r')
% %         if contPos(:,:,i) == 0 
% %             continue
% %         else
% %             scatter3(contPos(:,1,i),contPos(:,2,i),contPos(:,3,i),5,'*','MarkerFaceColor',orng)
% %         end
% %     end 
% % end

% Saves all data for statistical analysis
% % timeTotal = dt*((1:tMax)')/3600;
% endTime = toc;
% % fprintf('Complete!....  End Time:  %s \n',datestr(now))
% % fprintf('Elapsed time is %d seconds. \n',round(endTime))
% fname = sprintf('v%d_n%d_%2.1f_%3.4f_%d_%2.1f_%2.1f.mat',13,nCells,CGel,fiberDens,avgSearchTime,AI,RGD);
% save(fname);
% save(filename,'random_walk')
% save(filename1)

delete(gcp('nocreate'))
% end
% end
tclock = [floor(toc/3600), mod(floor(toc/60),60), mod(toc,60)];
    disp(tclock)
fprintf('Complete!  End Time:  %s \n',datestr(now))
% References
%Influence of Crosslink Density and Stiffness on Mechanical Properties of
%Type I Collagen Gel, Lin et al, 2015
