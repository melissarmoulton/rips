function Ustruct = odflowspeed(waveH, waveDir, wavePeriod, waveDepth, ...
                           depthVar, zBar, tidalElev, gammaBr, ...
                           vAlong, vAdjust)
% ODFLOWSPEED computes parameterized offshore-directed flow speed.
%  This function computes the parameterized the maximum offshore-directed
%  flow speed (Moulton et al., 2016) resulting from waves incident on non-
%  uniform bathymetry. It also returns the wave-breaking regime, sea-level
%  difference, alongshore flow speed, and suppression factor.
%
%  Syntax: Ustruct = odflowspeed(waveH, waveDir, wavePeriod, waveDepth, ...
%                            depthVar, zBar, tidalElev, gammaBr, ...
%                            vAlong, vAdjust)
%
%  Inputs:
%     waveH - Wave height (m) at depth waveDepth
%     waveDir - Wave direction (deg. from shore-normal) at depth waveDepth
%     wavePeriod - Wave period (s)
%     depthVar - Vertical scale of alongshore bathymetric variability (m)
%     zBar - Mean elevation on bar relative to arbitrary datum (m)
%     tidalElev - Tidal elevation relative to same arbitrary datum (m)
%     gammaBr - Wave breaking gamma value, typically 0.3<gammaBr,0.9
%     vAlong - Breaking-wave-driven alongshore current (if 1 value) OR
%              use vAlong = [dragC bSlope] (2 values) to compute with
%              quadratic drag law:
%                 dragC - Quadratic drag coefficient, ~3E-3
%                 bSlope - Approximate beach slope
%     vAdjust - Option to adjust sensitivity to alongshore flow, 0 to 1
%               (Beta in Moulton et al., 2016)
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
%         (FV in Moulton et al., 2016)
%
%  Subroutines: WAVESHOAL, FIND_C, FINDC0, FIND_L0, FIND_L, FIND_N
%  MAT-files required: none

% Melissa Moulton, Last updated: July 6, 2016

%% Constants

g=9.8; % m/s^2

%% Compute depthBar and depthChan

depthBar = tidalElev-zBar;
depthChan = depthBar+depthVar;

%% Compute wave properties at breaking

% Shoal the waves to breaking
[waveH_Br, waveDepth_Br, waveDir_Br] = deal(NaN*waveH);

for ii = 1:length(waveH)
    if ~isnan(wavePeriod(ii)+waveDepth(ii)+waveH(ii)+waveDir(ii))
        wavestruct = waveshoal(wavePeriod(ii),waveDepth(ii),waveH(ii),waveDir(ii),gammaBr);
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

% Initialize with bar-break value
D = 1/16*gammaBr^2.*(cosd(abs(waveDir_Br)).^2+1/2).*(waveH_Br/gammaBr-depthBar); 

% Set to zero if shore-break
D(R==0) = 0;

% Set to saturated value if saturated
D(R==2) = 1/16*gammaBr^2.*(cosd(abs(waveDir_Br(R==2))).^2+1/2).*(depthChan(R==2)-depthBar(R==2)); 

%% Compute alongshore flow

% Assign speed or drag coefficient and beach slope from vAlong

if numel(vAlong)==1
    V = vAlong;
else
   
dragC = vAlong(:,1);
bSlope = vAlong(:,2);
    
% Initialize with maximum alongshore flow for waves breaking on the bar
V = dragC.^(-1/2)*sqrt(5)/4*(gammaBr)^(3/4)*bSlope.^(1/2)*(gammaBr*depthBar/sqrt(2)).^(1/2).*sind(abs(waveDir_Br)).^(1/2).*sign(sind(waveDir_Br));

% For shore-break cases, fill in value based on wave height
V(R==0) = dragC.^(-1/2)*sqrt(5)/4*(gammaBr)^(3/4)*bSlope.^(1/2)*(waveH_Br/sqrt(2)).^(1/2).*sind(abs(waveDir_Br)).^(1/2).*sign(sind(waveDir_Br));

end

%% Compute alongshore flow suppression factor F

F = ((1+vAdjust^2.*abs(V).^2./(2*g*D)).^(1/2)-vAdjust*abs(V)./sqrt(2*g*D));
F(D==0) = NaN; % set to NaN if D=0

%% Compute U

U = sqrt(2*g*D).*F;
U(isnan(F)) = sqrt(2*g*D(isnan(F)));

%% Assign to output structure

Ustruct.Uoff = U;
Ustruct.Regime = R;
Ustruct.DeltaEta = D;
Ustruct.VAlong = V;
Ustruct.Factor = F;

end

%--------------------------------------------------------------------------

function wave = waveshoal(T, h0, H0, theta0, gamma)

% Find offshore group speed, phase speed, wavelength
cg0 = find_cg(T,h0);
c0 = find_c(T,h0);
L0 = find_L0(T);

% Find approx water depth at breaking,
% shoal waves on finer depth grid near this value
hBr_approx = round(10*gamma*H0*1.5)/10;
hlist = [h0:-.2:hBr_approx (hBr_approx+0.02):-0.02:0]';

% Calculate L, H, c, cg, n:

L = find_L(T,hlist);
c = find_c(T,hlist);
cg = find_cg(T,hlist);
n = find_n(T,hlist);
theta = asind(c./c0.*sind(theta0));
H = find_H(H0,cg0,cg,theta0,theta);

% Calculate breaking variables
breaking_index = abs(H./hlist-gamma)==min(abs(H./hlist-gamma));
breaking_depth = hlist(breaking_index);
breaking_height = H(breaking_index);
breaking_angle = theta(breaking_index);

% Store variables
wave.h = hlist;
wave.L = L;
wave.L0 = L0;
wave.H = H;
wave.c = c;
wave.cg = cg;
wave.n = n;
wave.theta = theta;
wave.breaking_depth = breaking_depth;
wave.breaking_height = breaking_height;
wave.breaking_angle = breaking_angle;

end

%--------------------------------------------------------------------------

function c = find_c(T, h)
% FIND_C computs c, the phase speed
% Syntax: c=find_c(T,h)
% Inputs: T - period (s), h - depth (m) 

L = find_L(T, h);

c = L./T;

end

%--------------------------------------------------------------------------

function cg = find_cg(T, h)
% FIND_CG computs cg, the group velocity
% Syntax: cg=find_cg(T,h)
% Inputs: T - period (s), h - depth (m) 

L = find_L(T, h);
k = 2.*pi./L;

cg = (L./(2.*T)).*(1+2.*(k).*h./sinh(2.*(k).*h));

end

%--------------------------------------------------------------------------

function L0 = find_L0(T)
% FIND_L0 computes L0, the offshore wavelength
% Syntax: L0 = find_L0(T)
% Inputs: T - period (s)

g = 9.81;
L0 = g*T.^2/(2*pi);

end

%--------------------------------------------------------------------------

function L = find_L(T, h)
% FIND_L computes L, the wavelength
% Syntax: L0 = find_L0(T)
% Inputs: T - period (s), h - depth (m)

L0 = find_L0(T);
% Approximate formula, Ole Madsen, MIT Coastal Engineering Course
L = L0.*(1-exp(-(2.*pi.*h./L0).^1.25)).^0.4;

% Alternative is to solve iteratively, but this is accurate enough
% for most purposes.

end

%--------------------------------------------------------------------------

function n = find_n(T, h)
% FIND_N computes N, the ratio of the group and phase speeds
% Syntax: n = find_n(T,h)
% Inputs: T - period (s), h - depth (m)

L = find_L(T,h);
k = 2*pi*L.^-1;

n = 1/2*(1+2*k.*h./(sinh(2*k.*h)));

end

%--------------------------------------------------------------------------

