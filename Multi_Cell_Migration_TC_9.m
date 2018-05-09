%Fiber based cell migration code - Multiple cell persistent random walk
%code written by Ben Yeoman
%created 10/10/2016
%last edited 10/17/2016
%code updated and edited by Tyler Collins
%last edited 5/7/2018

clear ALL
close all
rng shuffle

% File settings or data storage
filename = 'Multi_Cell_No_Interaction';
filename1 = 'Multi_Cell_No_Interaction_Data';
% rng(datenum(date));%......................................For repeatability
tic

% Model Definitions
runTime = 2*24*3600;%.......................Total runtime for model (sec)
dt = 2;%.......................................Size of each time step (sec)
elmtSize = 0.01;%.........................Size of binding site element (um)
polarize = 0.2;%..............................................Polarize cell
nCells = 3;%........................................Number of cells running
oriPos = 10;%.......................Distance cells initialize from origin

fprintf('Running %d cell/s over %d day/s. Start Time:  %s \n',nCells,...
    round(runTime/(24*3600)),datestr(now))

%Code Sectional Toggles 
gpuOn = 0;%................................Uses the GPU for matrix division
runAvgVelocity = 0;%...............Enable(1) or Disable(0) Average Velocity
runFitCurves = 0;%.......................Enable(1) or Disable(0) Fit Curves
runStatAnalysis = 0;%.................Enable(1) or Disable(0) Stat Analysis 
.....[Avg Velocity, MSD, Persistence L] Disabled it speeds up significantly

% Initiate matrices for MSD calculation
% random_walk = zeros(1+run_time/delta_t,3,num_cells);
t = ((0:1+runTime/dt-1)');
t_max = (floor((1+runTime/dt)/4));
t_msd = dt*((1:t_max)')/3600;

% Constant Cell Attributes
Vpseudo = 0.45;%...................Velocity of extending pseudopod (um/sec)
Lsearch = 0.3;%...........Area of searching pseudopod in contact with fiber
Fmax = 10e-9;%........................Max force provided by single cell (N)
Dcell = 15;%.............................................Cell Diameter (um)
Vcell = (4/3)*pi*(Dcell/2)^3;%.........................Cell's volume (um^3)
dRecp = 0.5*10^6;%..............................Integrin receptors per cell
% dR(1:nCells,1) = dRecp/saellipsoid(15/2);%...Integrin Receptor density (receptors/um^2)
dR = zeros(nCells,1);
tC = Lsearch/Vpseudo;%............................Binding site contact time
Bhalf = 200;%....................Number of bound ligands at half max. force
Bmin = 100;%...............Min bound ligands needed for pseudopod outgrowth
Bmax = 530;%.............Min bound ligands needed for pseudopod contraction
kcell = 1e-8;%..................................Cell spring constant (N/um)
contactRange = Dcell/2 + Dcell;%........Contact range from cell center (um)
% contactRange = 0;
angAvoid = 10;%..........Radial angle of avoidance from extending pseudopod Need paper value

% Given Constants Palsson Papers
etaSurf = 8e-8;%..............................Surface viscosity; u0(N*s/um)
etaCC = 2.5e-7;%.........Cell-cell viscosity; roughly 3*etaSurf; uC(N*s/um)
lambda = 7;%.................................................Given constant
xNot = sqrt(1/(2*lambda));%.................................Constant in Eqn
vNot = xNot*exp(-1*lambda*xNot^2);%.........................Constant in Eqn
sf = 25;%..................................................Mystery Constant
alpha = 1e-8;%...............................Adhesive strength constant (N) 1e-8
Fcomp = 30e-9;%........................Compressive force during contact (N)

% Constant ECM Attributes
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

% Array prealocation && Output Parameter Preallocations
MSD = zeros(t_max,nCells,1);
con_L = zeros(301,nCells,1);
aR_sq = zeros(301,nCells,1);
av_cs = zeros(301,nCells,1);
D = zeros;%.....................................Random motility coefficient
alp = zeros;%.........................................................Alpha

dist1 = zeros(nCells,1);
avgVel = zeros(nCells,1); %.............................Average cell speed

Lp_cos = zeros(nCells,1);
R1 = zeros(nCells,1); %.................................LpR Goodness of fit
R1_O = zeros(nCells,1);
R2 = zeros(nCells,1); %.................................MSD Goodness of fit
Lp_R = zeros(nCells,1); %.......................Persistence length by <R^2>
R3 = zeros(nCells,1);%..................................LpC Goodness of fit
R3_O = zeros(nCells,1);
snseCount  = 0;

% Multi-cell Preallocation
SE = zeros(nCells,1);
SN = zeros(nCells);
SNSE = zeros(nCells);
pseudoCellAng = zeros(nCells);
di = zeros(nCells);
d = zeros(nCells);
surf = zeros(nCells);
ecks = zeros(nCells);
nCellDist = zeros(nCells);
cellDist = zeros(nCells,3,nCells);

tsurf = zeros(nCells,1);
nSurf = zeros(nCells,1);
Vinst = zeros(nCells,1);
fb = zeros(nCells,1);
nFp = 0;
Fp = zeros(nCells,3,nCells);
F = zeros(nCells,1);
cVel = zeros(nCells,3);
a = ones(nCells,1)*Dcell/2;
b = ones(nCells,1)*Dcell/2;
c = ones(nCells,1)*Dcell/2;
rCell = zeros(nCells,1);
dmin = zeros(nCells,1);
SAcell = zeros(nCells,1);
Kprime = zeros(nCells,1);
FNet = zeros(nCells,3);
Fx = zeros(nCells,1);
Fy = zeros(nCells,1);
Fz = zeros(nCells,1);
Vx = zeros(nCells,1);
Vy = zeros(nCells,1);
Vz = zeros(nCells,1);
nuVel = zeros(nCells,3);

% Parameters
CGel = 5;%........................................Gel concentration (mg/ml)
fiberDens = 0.001;%.............................Fiber density [fibers/um^3] 
avgSearchTime = 8;%...............................Pseudopod extension time
FiberAI = 0;%..........................................Fiber alignment (AI)
RGD = 5.7;%................................RGD Peptides/monomer/tropocollagen
Fiber_D = sqrt((CGel*6.022e20*(Ltropo+0.067)*Dtropo^2)/...
    (1e12*fiberDens*805000*0.7*Lfiber*0.9));%........Average fiber diameter
dL = (RGD*.7*.9)/((Ltropo+0.067)*Dtropo);%RGD ligand density [lingand/um^2]
AI = 0;%....................................................Fiber Alignment
sigma = sqrt(-1642*reallog(AI));%....Angle deviation from orientation angle 
if isinf(sigma)
    sigma = 1e10;
end

crslnkPerL = ((4.439*sigma)/(sigma+4.55))/8;
crslnkPerF = Lfiber*crslnkPerL;
crslnkDens = fiberDens*crslnkPerF;%2.09-4.439 crosslinks per fiber.........
%............Lin et. al. 2015, Influence of Crosslink Density and Stiffness

%Gel stiffness as function of crosslink density - eqn from Lin et. al. 2015
if crslnkPerF < 3.5
    gel_stiffness = 1039.9*crslnkPerF-1992.9;
else
    gel_stiffness = 5247.9*crslnkPerF-16274;
end

%Converts gel stiffness to Pa/um^2
kecm = ones(nCells,1)*((1/crslnkDens^(1/3))*gel_stiffness/1e12);

%=========================================================================%
%Position and Phase Matrices
%=========================================================================%
cellPos = zeros(1,3,nCells);%..........................Cell position matrix
phaseRT = zeros;%............................Pseudopod phase runtime matrix
pseudoVect = zeros(nCells,3);%.............................Pseudopod vector
contPos = zeros(1,3,nCells);%..................Cell contact position matrix

%=========================================================================%
%Tracking variables
%=========================================================================%
Lpseudo = zeros(nCells,1);%......................Tracks length of pseudopod
newAngle = ones(nCells,1);%................Tracks when cell finds new fiber
retracting = ones(nCells,1);%.......................Tracks retracting phase
outgrowth = zeros(nCells,1);%........................Tracks outgrowth phase
contracting = zeros(nCells,1);%....................Tracks contracting phase
time = zeros(nCells,1);%...............Tracks pseudopod extension frequency
pos = 1;%............................Tracks cell position at each time step
cellBehav = zeros(nCells,4);%...............Tracks the behavior of the cell 
cPos = 1;%........................................Step for contact position

%=========================================================================%
%Temporary Variables
%=========================================================================%
fiberDir = zeros(nCells,3);%...................Unit vector of guiding fiber
fiberAngle = zeros(nCells,1);
searchTime = zeros(nCells,1);%.....................Pseudopod extension time
numElmts = zeros(nCells,1);%.....................Discrete elements of fiber
angle = 0;
Bsites = zeros(nCells,30);%.......................Distributed binding sites
searchLength = zeros(nCells,1);%.......Number of elements searched per step
bonds = zeros(nCells,1);%...................................Filopodia bonds
Bpseudo = zeros(nCells,1);%.............Number bonds along entire pseudopod
dMove = zeros(nCells,1);%................Distance cell moves in contraction
localFibers = zeros(nCells,1);%..............Number of fibers touching cell
cellAx = zeros(nCells,3);%.....................................Axes of cell
polAngle = 0;%...............................Cell polarity angle from fiber
dl = zeros(nCells,1);%................................Pseudopod search area
acAng = 0; 
bondsTotal = zeros(nCells,1);
dist_byV = zeros(nCells,3);
vPrev = zeros(nCells,3);
otX = zeros(nCells);%....................................Orientation factor
nuCellDist = zeros(nCells,3,nCells);
fv = zeros(nCells,1);
F0 = 0;
sl = 0;
ve = zeros(nCells,3);

%New Variables
cellContact = zeros(nCells,nCells);%Tracks if cell is within contact range
normCDist = zeros(1,3);

searchEnd = round(Lsearch/elmtSize)-1;

%Generates cell's initial position/polarity (+/- 50um from origin)
for cell=1:nCells

    cellPos(1,:,cell) = oriPos*(2*rand(1,3)-1);

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
end

%=========================================================================%
% DeltaT Loop - Calculates path of cell
%=========================================================================%
for i = 1:dt:runTime
   
    %Cell-Cell Distance Vector, Magnitude, Unit Vector, and Contact Range
    for cellOne = 1:nCells
        for cellTwo = 1:nCells
            if cellOne == cellTwo
                continue
            end
            cellDist(cellOne,:,cellTwo) = cellPos(pos,:,cellOne)- ...
                cellPos(pos,:,cellTwo);
            nCellDist(cellOne,cellTwo) = norm(cellDist(cellOne,:,cellTwo));
            nuCellDist(cellOne,:,cellTwo) = cellDist(cellOne,:,cellTwo)/...
                nCellDist(cellOne,cellTwo);
            
            %Establishes Cell Contact based on C-C distance
            if nCellDist(cellOne,cellTwo) <= contactRange && nCellDist(cellOne,cellTwo) > 0
                cellContact(cellOne,cellTwo) = 1;
            else
                cellContact(cellOne,cellTwo) = 0;
            end
        end
    end
%=========================================================================%
% Cell Loop  
%=========================================================================%
    for cell = 1:nCells
        
        % Generates new fiber direction
        if newAngle(cell)

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
                    newDir = 2*rand(1,3)-1;
                    fiber(j,:) = newDir(1,:)/norm(newDir);
                    
                %Calculates new fiber direction
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
            ve = sum(fbf)+vPrev(cell,:)/norm(vPrev(cell,:));

            %Cell follows acute angle between ve and prev direction
            veAngle = acosd(dot(vPrev(cell,:),ve)/(norm(vPrev(cell,:))*norm(ve)));
            if veAngle > 90
                ve = -ve;
            end

            %Get ratio of a0 to b0
            fbp = vertcat(vPrev(cell,:),fbf);%...........All fibers cell touches
            for j=1:size(fbp,1)
                cs(j,1) = abs(dot(ve,fbp(j,:))/(norm(ve)*norm(fbp(j,:))));              
                sn(j,1) = abs(sind(acosd(dot(ve,fbp(j,:))/(norm(ve)*...
                    norm(fbp(j,:)))))); 
            end
            a0 = sum(cs);   
            b0 = sum(sn);
            
            bondsTotal(cell) = localFibers(cell)*RGD*(0.7*(Dcell/2)/Ltropo)*(0.9*Fiber_D/Dtropo)*pi/2;
            keq = (bondsTotal(cell)*kI*kecm(cell))/(bondsTotal(cell)*kI+kecm(cell));

            %Define axes of the cell
            b(cell) = ((3*Vcell*(kcell+keq))/(4*pi*((a0/b0)*keq+kcell)))^(1/3);%Major axis of cell
            c(cell) = b(cell);
            a(cell) = (3*Vcell)/(4*pi*b(cell)^2);
            cellAx(cell,:) = [a(cell),b(cell),c(cell)];%.................Cell axes
            rCell(:) = 7.5;%.........Assuming radius needed in Palsson Eqns
            dmin(:) = rCell(1)*0.1*(-1);%.................Min distance from 
            ...........cell center to non-deformed edge of neighboring cell

            %Calculate cell sphericity
            SAcell(cell) = saellipsoid(cellAx(cell,:));%..........Cell surface area 
            sl = (pi^(1/3)*(6*Vcell)^(2/3))/SAcell(cell);%..Cell sphericity
            dR(cell) = dRecp/SAcell(cell);%.......Membrane integrin density
            polAngle = 180*sl^2;%............................Polarity angle
            
            %Fiber angles with respect to previous direction
            for j=1:localFibers(cell)
                angle(j,1) = acosd(dot(fiber(j,:),ve)/...
                    (norm(fiber(j,:))*norm(ve)));
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
            dmax = max([dL,dR(cell)]);
            Pb = ((kon*dmax)/(kon*dmax+koff))*(1-exp(-(kon*dmax+koff)*tC));
            Bsites(cell,1:n) = poissrnd(Pb*RGD*0.7*0.9*...
                (elmtSize/(Ltropo+0.067))*(Fdm/Dtropo)*...
                (pi/2),numElmts(cell),1);%........Distributed binding sites
            Bpseudo(cell) = Bpseudo(cell)+...
                sum(Bsites(cell,1:searchLength(cell)),2);
            
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
        
        % Angle Between Major Axis and Neighboring Cell Center & C-C Distance
        % And forcing retraction if a cell is infront of another cell
        for cellOne = 1:nCells
            for cellTwo = 1:nCells
                if cellOne == cellTwo
                    continue
                end
                pseudoCellAng(cellOne,cellTwo) = atan2d(norm...
                    (cross(cellDist(cellOne,:,cellTwo),...
                    fiberDir(cellOne,:))),dot(cellDist(cellOne,...
                    :,cellTwo),fiberDir(cellOne,:)));
                if pseudoCellAng(cellOne,cellTwo) < angAvoid &&...
                    cellContact(cellOne,cellTwo) == 1
                    bonds(cellOne) = 0; 
                end
            end
        end

        %=================================================================%
        % RETRACTING PHASE - pseudopod will retract if weak bonds are 
        % formed, or if new pseudopod starts extending
        %=================================================================%
        if bonds(cell) < Bmin || time(cell) > searchTime(cell)
            retracting(cell) = 1;
            phaseRT(pos,cell) = 1;

            %Resets values
            Lpseudo(cell) = 0;
            time(cell) = 0;
            bonds(cell) = 0;
            Bpseudo(cell) = 0;
            pseudoVect(cell,:) = zeros(1,3);
            newAngle(cell) = 1;
            contracting(cell) = 0;

        %=================================================================%
        % OUTGROWTH PHASE - pseudopod will continue growing as long as 
        % enough stable bonds are formed to allow for actin polymerization
        %=================================================================%
        elseif bonds(cell) >= Bmin && bonds(cell) < Bmax
            outgrowth(cell) = 1;
            phaseRT(pos,cell) = 2;
            contracting(cell) = 0;

            %Tracks distance and time pseudopod extends along current fiber
            pseudoVect(cell,:) = pseudoVect(cell,:)+dl(cell)*(fiberDir(cell,:));
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
                        (elmtSize/Ltropo+0.067)*(Fdm/Dtropo)*(pi/2)...
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
                pseudoVect(cell,:) = pseudoVect(cell,:)+dl(cell)*...
                    fiberDir(cell,:);
                Lpseudo(cell) = norm(pseudoVect(cell,:));
            end
        end
        %=================================================================%
        % END Retraction, Outgrowth, and Contraction Loop
        %=================================================================%
        
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

        Ltouching = mean(cellAx(cell,:))/4;%...Estmd length of fibers touching cell
        bondsTotal(cell) = localFibers(cell)*RGD*(0.7*Ltouching/Ltropo)*...
            (0.9*Fiber_D/Dtropo)*pi/2;
        F0 = (Fmax*Bpseudo(cell))/(Bhalf+Bpseudo(cell));%....Adhesion force
        if Lpseudo(cell) == 0
            kT_eff = 0;
        else
            kT_eff = (F0*kecm(cell)^2*Lpseudo(cell)^3)/(F0+kecm(cell)*Lpseudo(cell))^2;
        end
        if isnan(kT_eff) || isinf(kT_eff)
            error('NaN Inf')
        end
        fb(cell) = bondsTotal(cell)*((kecm(cell)*kI)/((kecm(cell)+kI)*...
            koff))*exp(-(kT_eff*kI)/(bondsTotal(cell)*kB*T));%................
        ..................................Friction due to bond dissociation
        fv(cell) = 6*pi*eta*a(cell)*Kprime(cell);%.........Viscous friction
        if ~isreal(fv(cell))
            disp(fv(cell))
        end

        %Internal radii from cellOne to cellTwo
        for cellOne = 1:nCells 
            for cellTwo = 1:nCells
                if cellOne == cellTwo
                    continue
                end
                di(cellOne,cellTwo) = sqrt((b(cellOne)*cos(pseudoCellAng(cellOne,cellTwo)))^2 ...
                    +(a(cellOne)*sin(pseudoCellAng(cellOne,cellTwo)))^2);
            end
        end

        %Distance between cell membranes
        for cellOne = 1:nCells 
            for cellTwo = 1:nCells
                if cellOne == cellTwo
                    continue
                end
                d(cellOne,cellTwo) = nCellDist(cellOne,cellTwo)...
                    -di(cellOne,cellTwo)-di(cellTwo,cellOne);
                ecks(cellOne,cellTwo) = (d(cellOne,cellTwo)-dmin(cellOne))/rCell(cellOne); %/rCell(cellOne)-a(cellOne);
            end
        end

        %Common Surface Area between cellOne and cellTwo
        for cellOne = 1:nCells 
            for cellTwo = 1:nCells
                if not(cellContact(cellOne,cellTwo)) || cellOne == cellTwo
                    surf(cellOne,cellTwo) = 0;
                elseif cellContact(cellOne,cellTwo)
                    surf(cellOne,cellTwo) = 0.25*exp(-5*...
                        (ecks(cellOne,cellTwo)-dmin(cellOne))^2)*...
                        (((sf+rCell(cellOne)*((1/di(cellOne,cellTwo))+...
                        (1/di(cellTwo,cellOne))))/(2+sf)));
                    surf(cellOne,cellTwo) = surf(cellOne,cellTwo)/...
                        SAcell(cellOne);%.......Normalizing surface-surface
                end
            end
        end
        
        %SE diagonal component of SNSE and norm
%         for cellTwo = 1:nCells
%                 if surf(cellTwo,:) == 0
%                     tsurf(cell) = 0;
%                 else
%                     tsurf(cell) = sum(surf(cell,:),2)/SAcell(cell);%Not
%                     ...needed because we use fb and fv for SE.
%                 end
%         end
%             SE(cell) = dmin(cell)*etaSurf + etaCC*tsurf(cell);%Palsson version of SE

        %SE calculated based on contact with other cells
        if nnz(cellContact(cell,:)) == 0
            SE(cell) = fb(cell)+fv(cell);
        else
            SE(cell) = fb(cell)*1;
        end

        %SN triangular component of matrix SNSE
        for cellTwo = 1:nCells
            if not(cellContact(cell,cellTwo)) || cell == cellTwo
                SN(cell,cellTwo) = 0;
            else
                SN(cell,cellTwo) = -1*surf(cell,cellTwo)*etaCC;
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
            otX(cell,cellTwo) = (rCell(cell)/2)*(1/di(cell,cellTwo) +...
                1/di(cellTwo,cell));
%             if  d(cell,cellTwo) > contactRange%ecks(cell,cellTwo)
            if ~cellContact(cell,cellTwo)
                Fp(cell,:,cellTwo) = [0 0 0];
            elseif ecks(cell,cellTwo) < 0 
                nFp = Fcomp*otX(cell,cellTwo)*(-1*ecks(cell,cellTwo))^(3/2);
                Fp(cell,:,cellTwo) = nFp*nuCellDist(cell,:,cellTwo);
            else
                nFp = -1*alpha*otX(cell,cellTwo)*((ecks(cell,cellTwo)+...
                    xNot)*exp(-1*lambda*(ecks(cell,cellTwo)+xNot)^2)-...
                    vNot*exp(-1*lambda*ecks(cell,cellTwo)^2));
                Fp(cell,:,cellTwo) = nFp*nuCellDist(cell,:,cellTwo);
            end
            if isnan(Fp(1,1,cell)) == 1 || isinf(Fp(1,1,cell)) == 1
                error('Fp is NaN or Inf')
            end
        end

        if contracting(cell)
            
            %Active force generation for contracting cell
            if F0 == 0 || Lpseudo(cell) == 0 
                F(cell) = 0;
            else
                F(cell) = (F0*kecm(cell)*Lpseudo(cell))/...
                    (F0+kecm(cell)*Lpseudo(cell));
            end
            
            %Sum of active and passive forces
            FNet(cell,:) = F(cell)*(pseudoVect(cell,:)/...
                norm(pseudoVect(cell,:)))+sum(Fp(cell,:,:),3);%.............................
            ...........................Force applied in pseudopod direction
        else
            FNet(cell,:) = sum(Fp(cell,:,:),3);
        end
    end
    %=====================================================================%
    % END cell loop but still in deltaT
    %=====================================================================%
    if nnz(contracting) == 0 && nnz(cellContact) == 0%..No rxns/No movement 
        cellPos(pos+1,:,:) = cellPos(pos,:,:);
    elseif nnz(cellContact) == 0%..........No cell contact/No need for SNSE
        for cell = 1:nCells
            if contracting(cell)
                Vinst(cell) = norm(FNet(cell,:))/SE(cell);
                dMove(cell) = Vinst(cell)*dt;%.............................
                .............Distance traveled during contractile time step
                dMove(cell) = min([Lpseudo(cell),dMove(cell)]);%...........
                .......................................Eliminates overshoot
                Lpseudo(cell) = Lpseudo(cell)-dMove(cell);  %specifically for contraction

                %Stores new cell position and direction
                cellPos(pos+1,:,cell) = dMove(cell,:)*(pseudoVect(cell,:)/...
                norm(pseudoVect(cell,:)))+cellPos(pos,:,cell);
            else
                cellPos(pos+1,:,cell) = cellPos(pos,:,cell);
            end
        end
    else
        if gpuOn == 1
            Fx = FNet(:,1);
            Fy = FNet(:,2);
            Fz = FNet(:,3);
%             Fx2 = gpuArray(Fx);
%             Fy2 = gpuArray(Fy);
%             Fz2 = gpuArray(Fz);
            SNSE2 = gpuArray(SNSE);
        else
            Fx = FNet(:,1);
            Fy = FNet(:,2);
            Fz = FNet(:,3);
        end
        snseCount = snseCount+1;

%         Fx = sparse(FNet(:,1));
%         Fy = sparse(FNet(:,2));
%         Fz = sparse(FNet(:,3));
%         SNSE = sparse(SNSE);
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

        if nnz(FNet) == 0
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
%                 Vx = linsolve(SNSE,Fx);
%                 Vy = linsolve(SNSE,Fy);
%                 Vz = linsolve(SNSE,Fz);
                
                Vx = mldivide(SNSE,Fx);
                Vy = mldivide(SNSE,Fy);
                Vz = mldivide(SNSE,Fz);
            end
        end

        for x = 1:nCells
            if isnan(Vx(x)) || isinf(Vx(x))
                error('isnan(Vx(x)) || isinf(Vx(x))')
%                     Vx(x) = 0;
            end
            if isnan(Vy(x)) || isinf(Vy(x))
                error('isnan(Vy(x)) || isinf(Vy(x))')
%                     Vy(x) = 0;
            end
            if isnan(Vz(x)) || isinf(Vz(x))
                error('isnan(Vz(x)) || isinf(Vz(x))')
%                     Vz(x) = 0;
            end
        end

        for u = 1:nCells
            Vinst(u) = norm([Vx(u),Vy(u),Vz(u)]);
            cVel(u,:) = [Vx(u),Vy(u),Vz(u)];%..........Cell velocity vector 
        end

        for v = 1:nCells
            if nnz(cVel(v,:)) == 0
                nuVel(v,:) = [0,0,0];
            else
                nuVel(v,:) = cVel(v,:)/norm(cVel(v,:));
            end
        end

        for cell = 1:nCells
%             if Vinst(cell) > 1
%                 pause;
% %                 error('high velocity')
%             end
            if Vinst(cell) > 1
%                 pause;
                error('high velocity')
            end
            dMove(cell) = Vinst(cell)*dt;%.................................
            .................Distance traveled during contractile time step
%             dMove(cell) = min([Lpseudo(cell),dMove(cell)]);%...............
            ...........................................Eliminates overshoot
            if contracting(cell)
%                 dMove(cell) = min([Lpseudo(cell),dMove(cell)]);%...............
            ...........................................Eliminates overshoot
                Lpseudo(cell) = Lpseudo(cell)-dMove(cell);  %specifically for contraction
            end
            %Stores new cell position and direction
            cellPos(pos+1,:,cell) = dMove(cell,:)*nuVel(cell,:)+cellPos(pos,:,cell);
        end
    end
    
    for cell = 1:nCells
        if cellPos(pos+1,:,cell) == 0
            error('cellPos = [0 0 0] for cell %d',cell) 
        end 
    end
    
    for cell = 1:nCells
        if cellPos(pos+1,:,cell) ~= cellPos(pos,:,cell)
            vPrev(cell,:) = cellPos(pos+1,:,cell)-cellPos(pos,:,cell);
        end

        if contracting(cell)
            if Lpseudo(cell) <= 0.05

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
        end
        if nnz(cellContact(cell,:)) ~= 0
            contPos(cPos,:,cell) = cellPos(pos,:,cell);
        end
    end
    
    %Reset Variables for next time step
    Vx = zeros(nCells,1);
    Vy = zeros(nCells,1);
    Vz = zeros(nCells,1);
    FNet = zeros(nCells,3);
    Vinst = zeros(nCells,1);
    SNSE = zeros(nCells,nCells);
    dMove = zeros(nCells,1);
    F = zeros(nCells,1);
    
    if nnz(cellContact) ~= 0
        cPos = cPos+1;
    end
    pos = pos+1;
    
    if pos == runTime/2*.25
        fprintf('25%%...');
    elseif pos == runTime/2*.5
        fprintf('50%%...');
    elseif pos == runTime/2*.75
        fprintf('75%%...');
    end
    
    for x = 1:nCells
        if isnan(cellPos(pos,:,x))
            error('NaN cell position')
        end
    end
    
    %Forcing retraction for following cells while leading cell is
    %contracting
    for cellOne = 1:nCells
        if contracting(cellOne)
            for cellTwo = 1:nCells
                if cellContact(cellOne,cellTwo)
                    
                    %Set time until new pseudopod extends
                    searchTime(cellTwo) = exprnd(avgSearchTime);
                    if searchTime(cellTwo) > 2*avgSearchTime
                        searchTime(cellTwo) = 2*avgSearchTime;
                    end
                    
                    Lpseudo(cellTwo) = 0;
                    time(cellTwo) = 0;
                    bonds(cellTwo) = 0;
                    Bpseudo(cellTwo) = 0;
                    pseudoVect(cellTwo,:) = zeros(1,3);
%                     newAngle(cellTwo) = 1;
                    contracting(cellTwo) = 0;
%                     retracting(cellTwo) = 1;
                end
            end
        end
    end
end
%=========================================================================%
% END deltaT loop
%=========================================================================%
randomWalk = cellPos;

%=========================================================================%
% Statistics Section
%========================================================================%
if runStatAnalysis == 1
    for i = 1 : nCells
        rw = randomWalk(:,:,i);
    
        %=================================================================%
        % Cell Velocity
        %=================================================================%
        if runAvgVelocity == 1
            pos = 1;
            w = 1;
            tm = (1:runTime/10:runTime)./3600;%................................10 time points to get distance
            dst = zeros;
            for q=1:(size(rw,1))/10:size(rw,1)
                Mag = zeros(size(rw,1),1,nCells);
                for j=1:q
                    if cellPos(j + 1,:,i) ~= cellPos(j,:,i)
                        Mag(pos,1,i) = norm(cellPos(j + 1,:,i) - cellPos(j,:,i));%.......Gets distances between each direction change
                        pos = pos + 1;
                    end
                end
                dst(w,1,i) = sum(Mag(:,1,i));%...............................Gets distance traveled at each time point
                w = w + 1;
            end
            
            if runFitCurves == 1
                try %..................................................Curve fit to get slope for average velocity
                    [f,g] = fit(tm',dst(:,1,i),'poly1');
                    c = coeffvalues(f);
                    avgVel(i) = c(1);%......................................Cell Speed
                    R1(i) = g.rsquare;%.....................................Goodness of fit
                catch
                    disp('Velocity Fitting Error')
                    avgVel(i) = 0;
                    R1(i) = 0;
                end
            end
        end

        %=====================================================================%
        % Non-Overlapping MSD
        %=====================================================================%
        msd = zeros;
        m = 1;
        for tau = 1:t_max

            %Gets nonoverlapping points for tau
            idx = mod(t,tau + 1) == 0;
            rw_idx = rw(idx,:);

            displacement = rw_idx(2:end,1:3) - rw_idx(1:end-1,1:3); %xyz displacement matrix
            D_squared = sum(displacement.^2,2); %squared displacement
            msd(m,1) = mean(D_squared); %MSD
            m = m + 1;
        end

        MSD(:,i) = msd;
        
        if runFitCurves == 1
            [f,g] = fit(t_msd,msd,'power1');
            c = coeffvalues(f);
            D(i) = c(1);
            alp(i) = c(2);
            R2(i) = g.rsquare;
        end

        %===================================================================================================================================%
        %Persistence Length code - Calculates Lp using <cos> and <R^2> for a cell's random walk
        %code written by Ben Yeoman and Parag Katira
        %created 10/8/2015
        %edited 11/9/2015 - BY
        %===================================================================================================================================%
        %Creates indices and variables
        n_pos = zeros;
        n_vect = zeros;
        n_mag = zeros;
        avg_cos = zeros;
        avg_R_squared = zeros;
        contour_length = zeros;

        %clear n_pos n_vect n_mag
        pos = 1;

        %Gets direction vector and magnitude for each change in position
        for j=1:size(rw,1)-1
            if rw(j+1,:) ~= rw(j,:)
                n_pos(pos,1:3) = rw(j,:);                                    %Position matrix for each new position
                n_vect(pos,1:3) = rw(j+1,:) - rw(j,:);                   %Vector matrix for each direction change
                n_mag(pos,1) = norm(rw(j+1,:) - rw(j,:));                %Magnitude matrix for each vector

                pos = pos+1;
            end
        end
        n_pos(pos,1:3) = rw(size(rw,1),1:3);

        xl = 60; %floor(sum(n_mag)/8);                                                      %x limit for Lp plots
        n = 0;                                                                     %Tracks steps for L0

        %Loop for incresing L0
        for L0=0:(xl/300):xl

            n = n + 1;                                                             %Increments each L0 step
            end_pos = 1;                                                           %First direction vector
            m = 0;                                                                 %Tracks cos(theta)/R for each L0 step
            stop = 0;                                                              %Breaks while loop if pos equals size of vector matrix
            total_cos = 0;
            total_R = 0;

            %Calculates average cos(theta)/R^2 for given L0
            while end_pos < size(n_vect,1)

                Length = 0;                                                        %Resets L after moving step length of L0
                m = m + 1;                                                         %Increments steps for cos(theta)/R
                strt_pos = end_pos;                                                %Gets starting position of L0 contour

                %Gets the contour position at end of length L0
                while Length <= L0
                    if end_pos == size(n_vect,1)                                   %Breaks loop at end of contour length
                        stop = 1;
                        break;
                    end

                    Length = Length + n_mag(end_pos);                              %Calculates length of contour
                    end_pos = end_pos + 1;                                         %Gets end position of L0 contour
                end

                if stop == 0
                    %cos(theta) between vectors at start and end of contour
                    cos_theta = dot(n_vect(end_pos-1,:),n_vect(strt_pos,:))/(n_mag(end_pos-1)*n_mag(strt_pos));
                    total_cos = total_cos + cos_theta;

                    %Distance R between start and end points of contour
                    R_squared = (norm(n_pos(end_pos,:) - n_pos(strt_pos,:)))^2;
                    total_R = total_R + R_squared;
                end
            end

            avg_cos_step = total_cos/m;
            avg_R_step = total_R/m;

            avg_cos(n,1) = avg_cos_step;                                  %Average cos(theta) for given contour length
            avg_R_squared(n,1) = avg_R_step;                             %Average R squared for given contour

            contour_length(n,1) = L0;
        end
        if runFitCurves == 1
            ft = fittype('2*(b1^2)*(exp(-(x/b1))-1+x/b1)');
            try %...............................................Curve fit to get persistence length
                [f,g] = fit(contour_length,avg_R_squared,ft,'StartPoint',1);
                c = coeffvalues(f);
                Lp_R(i) = c(1);%.........................................Persistence length
                R3(i) = g.rsquare;%.....................................Goodness of fit
            catch
                Lp_R(i) = 0;
                R3(i) = 0;
            end
        end

        con_L(:,i) = contour_length;
        aR_sq(:,i) = avg_R_squared;
        av_cs(:,i) = avg_cos;

        dist1(i) = sum(n_mag);
        avgVel(i) = sum(n_mag)/runTime;

    %         con_L(:,i,sample) = contour_length;
    %         aR_sq(:,i,sample) = avg_R_squared;
    %         av_cs(:,i,sample) = avg_cos;
    %     
    %         dist1(i,sample) = sum(n_mag);
    %         avg_vel(i,sample) = sum(n_mag)/run_time;
        %         try
        %             %Lp by cosine
        %             ft = fittype('exp(-x/Lp1)');
        %             [f,g] = fit(contour_length,avg_cos,ft,'StartPoint',0.1);
        %             c = coeffvalues(f);
        %             Lp_cos(i,sample) = c(1);
        %             R1(i,sample) = g.rsquare;
        %
        %             %Lp by average R^2
        %             ft = fittype('2*(b1^2)*(exp(-(x/b1))-1+x/b1)');
        %             [f1,g1] = fit(contour_length,avg_R_squared,ft,'StartPoint',0.1);
        %             c = coeffvalues(f1);
        %             Lp_R(i,sample) = c(1);
        %             R2(i,sample) = g1.rsquare;
        %         catch
        %         end   

    %     try
    %         %Generalized Least-Square Regression
    %         fun = @(b,x) exp(-x./b);
    %         beta0 = 1;
    %         mdl = NonLinearModel.fit(contour_length,avg_cos,fun,beta0);
    % 
    %         Lp_cos(i) = mdl.Coefficients{1,1};
    %         R1(i) = mdl.Rsquared.Adjusted;
    %         R1_O(i) = mdl.Rsquared.Ordinary;
    % 
    %         fun1 = @(b1,x) (2*(b1.^2)*(exp(-(x./b1))-1+x./b1));
    %         beta1 = 0.1;
    %         mdl = NonLinearModel.fit(contour_length,avg_R_squared,fun1,beta1);
    % 
    %         Lp_R(i) = mdl.Coefficients{1,1};
    %         R3(i) = mdl.Rsquared.Adjusted;
    %         R3_O(i) = mdl.Rsquared.Ordinary;
    % 
    % %             Lp_cos(i,sample) = mdl.Coefficients{1,1};
    % %             R1(i,sample) = mdl.Rsquared.Adjusted;
    % %             R1_O(i,sample) = mdl.Rsquared.Ordinary;
    % %     
    % %             fun1 = @(b1,x) (2*(b1.^2)*(exp(-(x./b1))-1+x./b1));
    % %             beta1 = 0.1;
    % %             mdl = NonLinearModel.fit(contour_length,avg_R_squared,fun1,beta1);
    % %     
    % %             Lp_R(i,sample) = mdl.Coefficients{1,1};
    % %             R2(i,sample) = mdl.Rsquared.Adjusted;
    % %             R2_O(i,sample) = mdl.Rsquared.Ordinary;
    %     catch
    %     end
    end
end

% Stores all cell migration data
% ===================================================================================================================================%

tclock = [floor(toc/3600), mod(floor(toc/60),60), mod(toc,60)];
%disp(tclock) %Commented out now

axis([-300 300 -300 300 -300 300])
hold on
grid on
daspect([1 1 1])
xlabel('x[\mum]')
ylabel('y[\mum]')
zlabel('z[\mum]')

orng  = [255 175 0]/255;

for i=1:nCells % 3D scatter plot
    scatter3(randomWalk(1,1,i),randomWalk(1,2,i),randomWalk(1,3,i),25,'filled','MarkerEdgeColor','g','MarkerFaceColor','g','LineWidth',2)
%     plot3(random_walk(1:50,1,i),random_walk(1:50,2,i),random_walk(1:50,3,i),'Color','c')
%     plot3(random_walk(:,1,i),random_walk(:,2,i),random_walk(:,3,i),'Color','b')
    plot3(randomWalk(:,1,1),randomWalk(:,2,1),randomWalk(:,3,1),'Color','m')
    if nCells > 1 
        plot3(randomWalk(:,1,2),randomWalk(:,2,2),randomWalk(:,3,2),'Color','b')
    end
    if i > 2
        plot3(randomWalk(:,1,i),randomWalk(:,2,i),randomWalk(:,3,i),'Color','k')
    end
    scatter3(randomWalk(end,1,i),randomWalk(end,2,i),randomWalk(end,3,i),25,'MarkerEdgeColor','r','MarkerFaceColor','r')
    if contPos(:,:,i) == 0 
        continue
    else
        scatter3(contPos(:,1,i),contPos(:,2,i),contPos(:,3,i),5,'*','MarkerFaceColor',orng)
    end
%     scatter3(cont_pos(:,1,i),cont_pos(:,2,i),cont_pos(:,3,i),'MarkerEdgeColor','r','MarkerFaceColor','r')
%     plot3(cont_pos(:,1,i),cont_pos(:,2,i),cont_pos(:,3,i),'Color','r')
end 

% Saves all data for statistical analysis
timeTotal = dt*((1:t_max)')/3600;
endTime = toc;
fprintf('Complete!....  End Time:  %s \n',datestr(now))
fprintf('Elapsed time is %d seconds. \n',round(endTime))
% save(filename,'random_walk')
% save(filename1)

% References
%Influence of Crosslink Density and Stiffness on Mechanical Properties of
%Type I Collagen Gel, Lin et al, 2015