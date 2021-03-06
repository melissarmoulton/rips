function Ustruct = rip_speed(waveH, waveDir, wavePeriod, waveDepth,...
                               depthVar, zBar, tidalElev, gammaBr,...
                               vAlong, vAdjust, ampCoeff)
% RIP_SPEED computes parameterized offshore-directed flow speed.
%  This function computes the parameterized the maximum offshore-directed
%  flow speed (Moulton et al., 2017) resulting from waves incident on non-
%  uniform bathymetry. It also returns the wave-breaking regime, sea-level
%  difference, alongshore flow speed, and suppression factor.
%
%  Syntax: Ustruct = odflowspeed(waveH, waveDir, wavePeriod, waveDepth, ...
%                            depthVar, zBar, tidalElev, gammaBr, ...
%                            vAlong, vAdjust, ampCoeff)
%
%  Inputs:
%     waveH - Wave height (m) at depth waveDepth
%     waveDir - Wave direction (deg. from shore-normal) at depth waveDepth
%     wavePeriod - Wave period (s)
%     depthVar - Vertical scale of alongshore bathymetric variability (m)
%     zBar - Mean elevation on bar relative to arbitrary datum (m)
%     tidalElev - Tidal elevation relative to same arbitrary datum (m)
%     gammaBr - Wave breaking gamma value, typically 0.3<gammaBr,0.9
%     vAlong - Breaking-wave-driven alongshore current (if 1 element) OR
%              use vAlong = [dragC bSlope] (2 elements) to compute with
%              quadratic drag law:
%                 dragC - Quadratic drag coefficient, ~3E-3
%                 bSlope - Approximate beach slope
%     vAdjust - Option to adjust sensitivity to alongshore flow, 0 to 1
%               (Beta in Moulton et al., 2016)
%     ampCoeff - Coefficient controlling amplitude,
%                Uoff=ampCoeff*sqrt(2*g*DeltaEta)*FV
%
%  Outputs:
%     Ustruct - structure with fields:
%     Uoff - Parameterized maximum offshore-directed flow speed (m/s)
%     Regime - Wave breaking regime
%       Regime = 0: shore-break
%       Regime = 1: bar-break
%       Regime = 2: saturated
%     DeltaEta - Alongshore sea-level difference
%     VAlong - Alongshore flow speed
%     Factor - Factor accounting for suppression of U by V, between 0 and 1
%     Inputs - All input fields
%         (FV in Moulton et al., 2017)
%
%  Subroutines: WAVESHOAL
%
% Melissa Moulton

%% Constants

g=9.8; % m/s^2

%% Compute depthBar and depthChan

depthBar = tidalElev-zBar;
depthBar(depthBar<0)=0; % If bar exposed, set depth to 0
depthChan = depthBar+depthVar;

%% Exclude times with wave angles between 90 and 270 from analysis
exind=intersect(find(abs(waveDir)>=90),find(abs(waveDir)<=270));
waveDir(exind)=NaN;

%% Compute wave properties at breaking (add: and at bar crest)

[waveH_Br, waveDepth_Br, waveDir_Br] = deal(NaN*ones(size(waveH)));

for ii = 1:length(waveH)
    if ~isnan(wavePeriod(ii)+waveDepth(ii)+waveH(ii)+waveDir(ii))
        wavestruct = waveshoal(wavePeriod(ii),waveDepth(ii),waveH(ii),...
                               waveDir(ii),gammaBr);
        waveH_Br(ii) = wavestruct.breaking_height;
        waveDepth_Br(ii) = wavestruct.breaking_depth;
        waveDir_Br(ii) = wavestruct.breaking_angle;
    end
end

%% Determine wave-breaking regime

R = ones(size(waveH)); % Initialize with bar-break value, 1
R(waveH_Br<gammaBr*depthBar) = 0; % Set to 0 if shore-break
R(waveH_Br>gammaBr*depthChan) = 2; % Set to 2 if saturated

%% Compute sea-level tilt

% Initialize with bar-break value,
% If waves break on the bar but not in the channel:
D = 1/16*gammaBr^2.*(cosd(abs(waveDir_Br)).^2+1/2).*...
    (waveH_Br/gammaBr-depthBar); 

% Set to zero if shore-break,
% No sea-level difference if waves break onshore of bathy variability:
D(R==0) = 0;

% Set to saturated value if saturated,
% If waves break in channel and on bar (offshore of bathy variability),
% the sea level difference are limited by the water depths
D(R==2) = 1/16*gammaBr^2.*(cosd(abs(waveDir_Br(R==2))).^2+1/2).*...
          (depthChan(R==2)-depthBar(R==2)); 

%% Compute alongshore flow

% Assign speed or drag coefficient and beach slope from vAlong

if numel(vAlong)==1
    V = vAlong;
else
   
dragC = vAlong(:,1);
bSlope = vAlong(:,2);
    
% Initialize with maximum alongshore flow for waves breaking on the bar,
% Alongshore flow speed near the bar crest is limited by depth on bar:
V = dragC.^(-1/2)*sqrt(5)/4*(gammaBr)^(3/4)*bSlope.^(1/2)*...
    (g*gammaBr*depthBar/sqrt(2)).^(1/2).*...
    sind(abs(waveDir_Br)).^(1/2).*sign(sind(waveDir_Br)); 

% For shore-break cases, fill in value based on wave height
V(R==0) = dragC.^(-1/2)*sqrt(5)/4*(gammaBr)^(3/4)*bSlope.^(1/2)*...
          (g*waveH_Br(R==0)/sqrt(2)).^(1/2).*...
          sind(abs(waveDir_Br(R==0))).^(1/2).*sign(sind(waveDir_Br(R==0)));

% Option to use value based on wave height for all regimes (uncomment if desired)
% V = dragC.^(-1/2)*sqrt(5)/4*(gammaBr)^(3/4)*bSlope.^(1/2)*...
%           (g*waveH_Br/sqrt(2)).^(1/2).*...
%           sind(abs(waveDir_Br)).^(1/2).*sign(sind(waveDir_Br));
      
end

%% Compute alongshore flow suppression factor F

F = ((1+vAdjust^2.*abs(V).^2./(2*g*D)).^(1/2)-vAdjust*abs(V)./sqrt(2*g*D));
F(D==0) = NaN; % set to NaN if D = 0

%% Compute U

U = ampCoeff*sqrt(2*g*D).*F;
U(isnan(F)) = ampCoeff*sqrt(2*g*D(isnan(F))); % Cases with D = 0

%% Assign to output structure

Ustruct.Uoff = U;
Ustruct.Regime = R;
Ustruct.DeltaEta = D;
Ustruct.VAlong = V;
Ustruct.Factor = F;

Ustruct.Inputs.waveH=waveH;
Ustruct.Inputs.waveDir=waveDir;
Ustruct.Inputs.wavePeriod=wavePeriod;
Ustruct.Inputs.waveDepth=waveDepth;
Ustruct.Inputs.depthVar=depthVar;
Ustruct.Inputs.zBar=zBar;
Ustruct.Inputs.tidalElev=gammaBr;
Ustruct.Inputs.vAlong=vAlong;
Ustruct.Inputs.vAdjust=vAdjust;
Ustruct.Inputs.ampCoeff=ampCoeff;

end

%--------------------------------------------------------------------------

function wave = waveshoal(T, h0, H0, theta0, gamma)

% Constants
g=9.8; % m/s^2

% Find a water depth near but offshore of breaking,
% shoal waves on finer depth grid onshore of this depth
hBr_approx = round(10*H0/gamma*2)/10;
h = [h0:-.2:hBr_approx (hBr_approx-0.02):-0.02:0]';

% Calculate wavelengths in deep water at depths h

% Deep water wavelength:
Ldeep = g*T.^2/(2*pi); 

% Wavelength, Ole Madsen approx:
L = Ldeep.*(1-exp(-(2.*pi.*h./Ldeep).^1.25)).^0.4;

% Calculate group and phase speeds at depths h
c = L./T; % Phase speed
k = 2.*pi./L; % Wavenumber
cg = (L./(2.*T)).*(1+2.*(k).*h./sinh(2.*(k).*h)); % Group velocity

% Calculate group and phase speeds at depth h0
c0 = c(1); % Phase speed at depth h0
cg0 = cg(1); % Phase speed at depth h0

% Compute wave height and angle at depths h
theta = asind(c./c0.*sind(theta0));
H = H0*sqrt(cg0./cg).*sqrt(cosd(abs(theta0))./cosd(abs(theta)));

% Calculate breaking variables
breaking_index = abs(H./h-gamma)==min(abs(H./h-gamma));
breaking_depth = h(breaking_index);
breaking_height = H(breaking_index);
breaking_angle = theta(breaking_index);

% Store variables
wave.h = h;
wave.L = L;
wave.Ldeep = Ldeep;
wave.H = H;
wave.c = c;
wave.cg = cg;
wave.theta = theta;
wave.breaking_depth = breaking_depth;
wave.breaking_height = breaking_height;
wave.breaking_angle = breaking_angle;

end

