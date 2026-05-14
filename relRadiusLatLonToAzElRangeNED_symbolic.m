function out = relRadiusLatLonToAzElRangeNED_symbolic()
% relRadiusLatLonToAzElRangeNED_symbolic
%
% Symbolically converts relative spherical-coordinate differences
% [Delta r, Delta latitude, Delta longitude] into local NED components,
% azimuth, elevation, and range, then computes the symbolic linearization.
%
% Requires Symbolic Math Toolbox.
%
% Angles are in radians.
%
% State:
%   x = [dr; dphi; dlambda]
%
% Measurement:
%   h(x) = [azimuth; elevation; range]
%
% Output:
%   out.NED          = exact symbolic [N; E; D]
%   out.h            = exact symbolic [az; el; range]
%   out.H            = exact symbolic Jacobian dh/dx
%   out.G_NED        = exact symbolic Jacobian d[N,E,D]/dx
%   out.G_NED_zero   = small-offset NED linearization matrix
%   out.NED_lin_zero = first-order NED approximation near zero offset
%   out.Hbar         = symbolic Jacobian evaluated about xbar
%   out.hbar         = measurement evaluated about xbar
%   out.h_lin        = first-order measurement linearization about xbar

    %% Symbolic variables
    syms r1 phi1 dr dphi dlambda real
    syms dr_bar dphi_bar dlambda_bar real

    assumeAlso(r1 > 0);

    % State vector
    x = [dr; dphi; dlambda];

    % Linearization point
    xbar = [dr_bar; dphi_bar; dlambda_bar];

    %% Target spherical coordinates relative to observer
    r2 = r1 + dr;
    phi2 = phi1 + dphi;

    % Only longitude difference matters
    dlam = dlambda;

    %% Exact NED line-of-sight components
    N = r2 * ( ...
        cos(phi1)*sin(phi2) ...
        - sin(phi1)*cos(phi2)*cos(dlam) );

    E = r2 * cos(phi2) * sin(dlam);

    D = r1 - r2 * ( ...
        cos(phi1)*cos(phi2)*cos(dlam) ...
        + sin(phi1)*sin(phi2) );

    NED = simplify([N; E; D]);

    %% Azimuth, elevation, range
    horizontalRange = sqrt(N^2 + E^2);

    range = sqrt(N^2 + E^2 + D^2);

    az = atan2(E, N);

    % In NED, D is positive downward, so elevation uses -D
    el = atan2(-D, horizontalRange);

    h = simplify([az; el; range]);

    %% Symbolic Jacobians
    G_NED = simplify(jacobian(NED, x));

    H = simplify(jacobian(h, x));

    %% Small-offset NED linearization about dr = dphi = dlambda = 0
    G_NED_zero = simplify(subs(G_NED, ...
        [dr, dphi, dlambda], ...
        [0, 0, 0]));

    NED_lin_zero = simplify(G_NED_zero * x);

    %% General symbolic linearization about xbar
    hbar = simplify(subs(h, ...
        [dr, dphi, dlambda], ...
        [dr_bar, dphi_bar, dlambda_bar]));

    Hbar = simplify(subs(H, ...
        [dr, dphi, dlambda], ...
        [dr_bar, dphi_bar, dlambda_bar]));

    h_lin = simplify(hbar + Hbar * (x - xbar));

    %% Return output structure
    out.symbols.r1 = r1;
    out.symbols.phi1 = phi1;
    out.symbols.dr = dr;
    out.symbols.dphi = dphi;
    out.symbols.dlambda = dlambda;

    out.x = x;
    out.xbar = xbar;

    out.NED = NED;
    out.h = h;

    out.G_NED = G_NED;
    out.H = H;

    out.G_NED_zero = G_NED_zero;
    out.NED_lin_zero = NED_lin_zero;

    out.hbar = hbar;
    out.Hbar = Hbar;
    out.h_lin = h_lin;

    %% Display key results
    fprintf('\nExact NED vector [N; E; D]:\n');
    pretty(NED)

    fprintf('\nExact measurement vector h = [az; el; range]:\n');
    pretty(h)

    fprintf('\nExact symbolic Jacobian H = dh/d[dr, dphi, dlambda]:\n');
    pretty(H)

    fprintf('\nSmall-offset NED linearization matrix about zero offset:\n');
    pretty(G_NED_zero)

    fprintf('\nSmall-offset NED approximation:\n');
    pretty(NED_lin_zero)

    fprintf('\nGeneral linearized measurement model about xbar:\n');
    pretty(h_lin)
end