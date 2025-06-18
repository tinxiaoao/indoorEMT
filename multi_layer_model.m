function result = multi_layer_model(frequency, wall_type, t_mid, epsTab)
% MULTI_LAYER_MODEL  Strict Balanis (Sec. 5.4.4) implementation for a
% 3-layer wall under oblique incidence (ASCII-only MATLAB).
% -------------------------------------------------------------------------
% frequency : column-vector [Hz]
% wall_type : 1-external | 2-internal | 3-ceiling/floor | 4-window
% t_mid     : mid-layer thickness [m] (for window = total glass thickness)
% epsTab    : Cole-Cole parameter table  (N x 5)
% -------------------------------------------------------------------------
% Returns "result" structure with fields (size F x 1):
%   SUM_Solid_Angle_mean_R   mean reflection solid-angle integral
%   SUM_Solid_Angle_mean_T   mean transmission solid-angle integral
%   SUM_Solid_Angle_mean_R_T SUM_R - SUM_T
% -------------------------------------------------------------------------

mu0  = 4*pi*1e-7;           % vacuum permeability
eps0 = 1/(36*pi)*1e-9;      % vacuum permittivity
f    = frequency(:);        % ensure column vector
F    = numel(f);            % # frequency points

%% --- layer definition ---------------------------------------------------
switch wall_type
    case 1      % external wall  (15-19-15)
        id = [15 19 15];  t = [0.005 t_mid 0.005];
    case 2      % internal wall  (15-16-15)
        id = [15 16 15];  t = [0.005 t_mid 0.005];
    case 3      % ceiling / floor (15-19-15, mid=0.1 m)
        id = [15 19 15];  t = [0.005 0.1 0.005];
    case 4      % window glass-air-glass (13-21-13)
        id = [13 21 13];
        airGap = 0.012;
        g      = (t_mid - airGap)/2;
        t      = [g airGap g];
    otherwise
        error('wall_type must be 1-4');
end
L = 3;                             % # dielectric slabs

%% --- epsilon_r(f) for each layer ---------------------------------------
eps_r = zeros(F,L);
for k = 1:L
    p = epsTab(id(k),:);           % [e_inf e_s sigma_s tau alpha]
    eps_r(:,k) = p(1) + (p(2)-p(1))./(1 + (1i*2*pi*f*p(4)).^(1-p(5))) ...
                 + p(3)./(1i*2*pi*f*eps0);
end

%% --- integrate over 0..pi/2 (1000 pts) ---------------------------------
Ntheta = 1000;
thetaV = linspace(0, pi/2, Ntheta);
dTheta = thetaV(2) - thetaV(1);

sumR = zeros(F,1);
sumT = zeros(F,1);

dFull = [0 t 0];                 % include exterior air layers (1 x (L+2))
epsAir = ones(F,1);

omega = 2*pi*f;

for theta = thetaV
    % full permittivity stack (F x (L+2))
    epsFull = [epsAir eps_r epsAir];

    % propagation constant gamma (strictly without angle correction)
    gamma = 1i * omega .* sqrt(mu0 * eps0 .* epsFull);

    % Snell's law: compute angles theta_j
    theta_j = zeros(F, L+2);
    theta_j(:,1) = theta;
    for n = 2:(L+2)
        theta_j(:,n) = asin(sqrt(epsFull(:,n-1))./sqrt(epsFull(:,n)) .* sin(theta_j(:,n-1)));
    end

    % phase thickness psi_j = d_j * gamma_j * cos(theta_j)
    phi = gamma .* (dFull .* cos(theta_j));

    % initialise ABCD (region Nst = L+2 is outer air)
    Nst = L + 2;
    A = ones(F,Nst);  B = zeros(F,Nst);
    C = ones(F,Nst);  D = zeros(F,Nst);

    for n = Nst-1:-1:1
        Y = (cos(theta_j(:,n+1))./cos(theta_j(:,n))) .* sqrt(epsFull(:,n+1)./epsFull(:,n));
        Z = (cos(theta_j(:,n+1))./cos(theta_j(:,n))) .* sqrt(epsFull(:,n)./epsFull(:,n+1));

        eP = exp(phi(:,n));
        eM = exp(-phi(:,n));

        A(:,n) = 0.5 .* eP .* ( A(:,n+1).*(1+Y) + B(:,n+1).*(1-Y) );
        B(:,n) = 0.5 .* eM .* ( A(:,n+1).*(1-Y) + B(:,n+1).*(1+Y) );
        C(:,n) = 0.5 .* eP .* ( C(:,n+1).*(1+Z) + D(:,n+1).*(1-Z) );
        D(:,n) = 0.5 .* eM .* ( C(:,n+1).*(1-Z) + D(:,n+1).*(1+Z) );
    end

    % Reflection / Transmission (Balanis equations 5-90, 5-91) 
    % Note:  Equation 5-91 in book is wrong.
    RV = B(:,1)./A(:,1);   TV = 1./A(:,1);   % perpendicular
    RH = D(:,1)./C(:,1);   TH = 1./C(:,1);   % parallel

    w = cos(theta)*sin(theta)*dTheta;
    sumR = sumR + (1 - 0.5*(abs(RV).^2 + abs(RH).^2)).*w;
    sumT = sumT + 0.5*(abs(TV).^2 + abs(TH).^2).*w;
end

result.SUM_Solid_Angle_mean_R   = sumR;
result.SUM_Solid_Angle_mean_T   = sumT;
result.SUM_Solid_Angle_mean_R_T = sumR - sumT;

end
