clc;
clear;
close all;
warning('off', 'all');
 
ENA_TROPO_ERR_CORR = true; % Tropospheric error estimation by Saastamoinen model
 
%% Read given data
eph = importdata(fullfile('Data','eph.dat'));
eph = sortrows(eph,[1, 2]);
 
% receiver
rcvr = importdata(fullfile('Data','rcvr.dat'));
rcvr = sortrows(rcvr,[1, 2]);
 
% select common satellite
[~, ephIdx] = ismember(eph(:,2), rcvr(:,2));
[~, rcvrIdx] = ismember(rcvr(:,2), eph(:,2));
eph = eph(ephIdx, :);
rcvr = rcvr(rcvrIdx, :);
 
% Given 
init_pos = [-2694685.473; -4293642.366; 3857878.924]; % initial position
voflight = 299792458.0; % speed of light, m/s
wedot = 7.2921151467e-5; % earth's rotation rate, wgs84, r/s
mu = 3.986005e+14; % earth's universal gravitation, wgs84, m^3/s^2
F = -4.442807633e-10; % relativistic correction term
tar_XR = [-2700400; -4292560; 3855270];  % target receiver position
 
[wlat, wlon, walt] = wgsxyz2lla(tar_XR);
tar_XR_lla = [wlat, wlon, walt]'; % target receiver position
 
% Tropospheric constant
% Standard atmosphere - Berg, 1948 (Bernese)
% pressure [mbar]
% temperature [K]
% numerical constants for the algorithm [-] [m] [mbar]
% temperature at sea-level
Pr = 1013.25;
Tr = 291.15;
Hr = 50.0;
temp0 = 15;
 
 
% Main program
% Calculate satellite ECEF position vector
XS = nan(size(eph,1), 3); 
dtS = nan(size(eph,1), 1); 
for i = 1:size(eph,1)
    
% assign requried variables 
    pr      =   rcvr(i, 3);
    t       =   eph(i, 1) - pr/voflight; % transmission time
    svid    =   eph(i, 2);
    toc     =   eph(i, 3);
    toe     =   eph(i, 4);
    af0     =   eph(i, 5);
    af1     =   eph(i, 6);
    af2     =   eph(i, 7);
    ura     =   eph(i,8);
    e       =   eph(i, 9);
    sqrta   =   eph(i, 10);
    dn      =   eph(i, 11);
    m0      =   eph(i, 12);
    omega   =   eph(i, 13);
    omg0    =   eph(i, 14);
    i0      =   eph(i, 15);
    odot    =   eph(i, 16);
    idot    =   eph(i, 17);
    cus     =   eph(i, 18);
    cuc     =   eph(i, 19);
    cis     =   eph(i, 20);
    cic     =   eph(i, 21);
    crs     =   eph(i, 22);
    crc     =   eph(i, 23);  
    iod     =   eph(i,24);  
     
    % clokck offset 
    dt_sv = af0 + af1 * (t - toc) + af2 * (t - toc).^2;
    dtS(i,1) = dt_sv;
    
    % ECEF position 
    a = sqrta^2;   % semimajor axis
    
    n = sqrt(mu / a^3) + dn;   % corrected mean motion (rad/sec)
    
    tk = t - toe;   % time from ephemeris epoch
 
      if(tk >= 302400)
        tk = tk - 604800;
    elseif(tk <= -302400)
        tk = tk + 604800;
    end;
    
    mk = m0 + n*(tk);   
    
    % eccentrix anomaly
    max_iter = 4; % set maximum number of iterations
    iter = 0;
    EK = mk;
    EK = EK;
    while iter < max_iter || abs(EK-EK)>1e-12
        EK = EK;
        EK = EK - (EK - e * sin(EK) - mk) / (1 - e * cos(EK));
        iter = iter + 1;
    end
    
        
    fk = atan2(sqrt(1 - e^2) * sin(EK), cos(EK) - e);   % true anomaly
    
    phik = fk + omega;   
    
    duk = cuc * cos(2 * phik) + cus * sin(2 * phik);   % argument of latitude correction
    drk = crc * cos(2 * phik) + crs * sin(2 * phik);   % radius correction
    dik = cic * cos(2 * phik) + cis * sin(2 * phik);   % inclination correction
    
    uk = phik + duk;   % corrected argument of latitude
    rk = a * (1 - e * cos(EK)) + drk;   % corrected radius
    ik = i0 + dik + idot * tk;   % corrected inclination 
    
    omgk = omg0 + (odot - wedot) * tk - wedot * toe; % corrected longitude of node
    
    % in-plane x and y position
    xp = rk * cos(uk);   
    yp = rk * sin(uk);   
    
    % ECEF x,y,z-coordinate
    xs = xp * cos(omgk) - yp * cos(ik) * sin(omgk);   
    ys = xp * sin(omgk) + yp * cos(ik) * cos(omgk);   
    zs = yp * sin(ik);   
 
    
    % earth rotation correction
    % rotation angle
    omegatau = wedot * pr / voflight;
    
    % rotation matrix
    R3 = [ cos(omegatau)    sin(omegatau)   0;
          -sin(omegatau)    cos(omegatau)   0;
           0                0               1];
    XSpot = R3 * [xs, ys, zs]';
    
    XS(i, :) = XSpot';  % store corresponding satellite position
end
 
% least squares to estimate receiver location 
XR = init_pos;  
dtR = 0;  
[wlat, wlon, walt] = wgsxyz2lla(XR); 
x = [XR; dtR]; 
iter = 0;  
posErr = norm(XR-tar_XR); 
 
% printing
dryRunTbl = [iter, nan, nan, nan, nan, nan, XR', dtR, dtR*voflight, wlat, wlon, walt, nan]; 
fprintf('Iter#%d (initial):\n\n', iter);
fprintf('Initial position :\n'  );
fprintf('ECEF(m): %.3fm, %.3fm, %.3fm \n', XR);
fprintf('(WGS84 LLA): %.9f°, %.9f°, %.3fm\n', wlat, wlon, walt);
fprintf('\n  Total position error is %.3fm.\n\n', posErr);
fprintf('\n------------------------------------------------------------------\n');
 
while 1
    % resolve the troposheric error 
    if ENA_TROPO_ERR_CORR
        humi = 1.0;
        [wlat, wlon, walt] = wgsxyz2lla(XR);
        tropo_err = zeros(size(XS, 1), 1);
        for i = 1:size(XS, 1)
            enu = wgsxyz2enu(XS(i,:)', wlat, wlon, walt);
            el = asin(enu(3)/sqrt(enu(1)^2 + enu(2)^2));
            if walt < 0
                hgt = 0;
            else
                hgt = walt;
            end
            pres = Pr * (1 - 2.2557e-5 * hgt)^5.2568;
            temp = temp0 - 6.5e-3 * hgt + 273.16;
            e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));
            
            % saastamoninen model
            z = pi / 2.0 - el;
            trph = 0.0022768 * pres / (1.0-0.00266*cos(2.0 * wlat) - 0.00028 * hgt / 1e3) / cos(z);
            trpw = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z);
            tropo_err(i, 1) = trph + trpw;
        end
    end
    % troposheric error (saastamoinen model)
    
    pr = rcvr(:, 3); 
    
    dist = sqrt(sum((XS - XR').^2, 2)); % geometric distance
    b = dist + voflight .* (dtR - dtS) + tropo_err; % A. pseudorange
    
    H = [[XR' - XS]./dist, ones(size(XS, 1),1).*voflight];  % design matrix
    
    dx = inv(H' * H) * H' * (pr - b); % LS 
    
    residual = pr - (b + H*dx); % residual
    residualSE = residual' * residual; % residual, squared error
    
    % new solution
    x = x + dx; 
    XR = XR + dx(1:3);
    [wlat, wlon, walt] = wgsxyz2lla(XR); 
    dtR = dtR + dx(4); 
    posErr = norm(XR-tar_XR);     
    iter = iter + 1; 
    
    % printing
    dryRunTbl = [dryRunTbl; iter, dx', residualSE, XR', dtR, dtR*voflight, wlat, wlon, walt, posErr];
    fprintf('Iter #%d:\n\n', iter);
    fprintf('Updated receiver position: \n\n');
    fprintf('ECEF(m): %.3fm, %.3fm, %.3fm \n', XR);
    fprintf('WGS84 LLA: %.9f°, %.9f°, %.3fm \n', wlat, wlon, walt);
    fprintf('  Updated receiver clock offset is %.7fs (%.3fm)\n', dtR, dtR*voflight);
    fprintf('\n  Total position error is %.3fm\n', posErr);
    fprintf('  LS residual & squared error is %.3fm^2\n\n', residualSE);
    fprintf('\n------------------------------------------------------------------\n');
    
    if norm(dx(1:3)) < 1e-2 || ...  
            iter > 2  
       break; 
    end
end
 
% summary
fprintf('Summary: \n');
fprintf('\nInitial position ECEF are %.3fm, %.3fm, %.3fm', init_pos);
fprintf('\nFinal position ECEF are %.3fm, %.3fm, %.3fm\n',  XR);
fprintf('\nTarget position ECEF are %.3fm, %.3fm, %.3fm, positioning error: %.3fm\n', tar_XR, norm([tar_XR-XR]));
fprintf('\nReceiver clock offset are %.7fs, %.3fm\n', dtR, dtR*voflight);
fprintf('\n\n 19093616r AAE6102');
 
 
