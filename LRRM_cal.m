function [output] = LRRM_cal(varargin)

% Classical LRRM calibration
% Ref.: [1990 ARFTG] by A. Davidson and [2006 ARFTG] by L. Hayden

files                                           = varargin{1,1};
params                                          = varargin{1,2};

path                                            = files.path;

% load in data for debug purposes
% [thru, freq]                                    = readTouchStone_statistiCal([path files.thru]);
% [short,d]                                       = readTouchStone_statistiCal([path files.short]);
% [open,d]                                        = readTouchStone_statistiCal([path files.open]);
% [load,d]                                        = readTouchStone_statistiCal([path files.load]);

% [thru, freq]                                    = readTouchStone([path files.thru]);
% [short,d]                                       = readTouchStone([path files.short]);
% [open,d]                                        = readTouchStone([path files.open]);
% [load,d]                                        = readTouchStone([path files.load]);


if strcmp(files.thru(end-2:end), 'txt')
    [thru, freq]                                = readTouchStone_statistiCal([path files.thru]);
else
    [thru, freq]                                = readTouchStone([path files.thru]);
end

if strcmp(files.short(end-2:end), 'txt')
    [short, d]                                  = readTouchStone_statistiCal([path files.short]);
else
    [short, d]                                  = readTouchStone([path files.short]);
end

if strcmp(files.open(end-2:end), 'txt')
    [open, d]                                   = readTouchStone_statistiCal([path files.open]);
else
    [open, d]                                   = readTouchStone([path files.open]);
end

if strcmp(files.load(end-2:end), 'txt')
    [load, d]                                   = readTouchStone_statistiCal([path files.load]);
else
    [load, d]                                   = readTouchStone([path files.load]);
end

numpts                                          = length(freq);

if isfield(files, 'sw')
    file_sw                                     = files.sw;
    if strcmp(file_sw(end-2:end), 'cti')
        [sw, d]                                 = readCti_switchTerm([path files.sw]); % assume GammaF GammaR
        sw                                      = [sw(:,2) sw(:,1) zeros(numpts, 2)];
    else
        [sw, d]                                 = readTouchStone_statistiCal([path files.sw]);
    end    
else
    sw                                          = zeros(numpts,4);
end
    
Rdc                                             = params.Rdc;    

% correct for switch-terms
thru_sw                                         = correct_switchTerms(thru, sw);
short_sw                                        = correct_switchTerms(short, sw);
open_sw                                         = correct_switchTerms(open, sw);
load_sw                                         = correct_switchTerms(load, sw);

Tthru                                           = S2ABCD_Zc(thru_sw, 50*ones(numpts,1));

zshort                                          = S2Z(short_sw);
zopen                                           = S2Z(open_sw);
zload                                           = S2Z(load_sw);
% thru measurement
t1                                              = squeeze(Tthru(1,1,:));
t2                                              = squeeze(Tthru(1,2,:));
t3                                              = squeeze(Tthru(2,1,:));
t4                                              = squeeze(Tthru(2,2,:));

zxshort                                         = squeeze(zshort(1,1,:));
zyshort                                         = squeeze(zshort(2,2,:));

zxopen                                          = squeeze(zopen(1,1,:));
zyopen                                          = squeeze(zopen(2,2,:));

zxload                                          = squeeze(zload(1,1,:));
zyload                                          = squeeze(zload(2,2,:));

if strcmp(params.reflect, 'Load')
    % % a1, a2, and Va using match
    a1                                          = t1.*zyload - t2 + t3.*zxload.*zyload - t4.*zxload;
    a2                                          = 2*t4 - 2*t3.*zyload;
    Va                                          = 2*t1.*zxload.*zyload - 2*t2.*zxload;
elseif strcmp(params.reflect, 'Short')
    %eq(24) in Leonard 2006 ARFTG
    % a1, a2, and Va using short
    a1                                          = t1.*zyshort - t2 + t3.*zxshort.*zyshort - t4.*zxshort;
    a2                                          = 2*t4 - 2*t3.*zyshort;
    Va                                          = 2*t1.*zxshort.*zyshort - 2*t2.*zxshort;
end

% b1, b2, and Vb using open
b1                                              = t1.*zyopen - t2 + t3.*zxopen.*zyopen - t4.*zxopen;
b2                                              = 2*t4 - 2*t3.*zyopen;
Vb                                              = 2*t1.*zxopen.*zyopen - 2*t2.*zxopen;

P1                                              = (Va.*b2 - Vb.*a2)./(a1.*b2 - a2.*b1);
P2                                              = (Vb.*a1 - Va.*b1)./(a1.*b2 - a2.*b1);
% two solutions
A_over_C1                                       = (P1 + sqrt(P1.^2 - 4*P2))/2;
B1                                              = (P1 - sqrt(P1.^2 - 4*P2))/2;

A_over_C2                                       = (P1 - sqrt(P1.^2 - 4*P2))/2;
B2                                              = (P1 + sqrt(P1.^2 - 4*P2))/2;

z_load_over_short1                              = ((B1 - zxload)./(zxload - A_over_C1))./( (B1 - zxshort)./(zxshort - A_over_C1)  );
z_load_over_open1                               = ((B1 - zxload)./(zxload - A_over_C1))./( (B1 - zxopen)./(zxopen - A_over_C1)  );

z_load_over_short2                              = ((B2 - zxload)./(zxload - A_over_C2))./( (B2 - zxshort)./(zxshort - A_over_C2)  );
z_load_over_open2                               = ((B2 - zxload)./(zxload - A_over_C2))./( (B2 - zxopen)./(zxopen - A_over_C2)  );

% traditional LRRM for the root choice
Ls1                                             = (imag(z_load_over_open1/sqrt(-1))./real(z_load_over_open1/sqrt(-1)))./( (2*pi*freq)./Rdc );
Ls2                                             = (imag(z_load_over_open2/sqrt(-1))./real(z_load_over_open2/sqrt(-1)))./( (2*pi*freq)./Rdc );

zaload1                                         = Rdc + sqrt(-1)*2*pi*freq.*Ls1;
zaload2                                         = Rdc + sqrt(-1)*2*pi*freq.*Ls2;

zashort1                                        = zaload1./z_load_over_short1;
zashort2                                        = zaload2./z_load_over_short2;

zaopen1                                         = zaload1./z_load_over_open1;
zaopen2                                         = zaload2./z_load_over_open2;
for i = 1:numpts
    if imag(zashort1(i))<0 % assumes the short is equivalent to a negative inductance 
        A_over_C(i,1)                           = A_over_C1(i);
        B(i,1)                                  = B1(i);
        Ls(i,1)                                 = Ls1(i);                
        zaload(i,1)                             = zaload1(i);                
        zashort(i,1)                            = zashort1(i);
        zaopen(i,1)                             = zaopen1(i);
        z_load_over_open(i,1)                   = z_load_over_open1(i,1);
        z_load_over_short(i,1)                  = z_load_over_short1(i,1);        
    else
        A_over_C(i,1)                           = A_over_C2(i);
        B(i,1)                                  = B2(i);
        Ls(i,1)                                 = Ls2(i);                
        zaload(i,1)                             = zaload2(i);                
        zashort(i,1)                            = zashort2(i);
        zaopen(i,1)                             = zaopen2(i);        
        z_load_over_open(i,1)                   = z_load_over_open2(i,1);        
        z_load_over_short(i,1)                  = z_load_over_short2(i,1);        
    end        
end 
    
z_short_over_load                               = 1./z_load_over_short;


if params.Ls == 0 % Ls_mean from the measurements
    Ls_mean                                     = mean(Ls(find_index(freq, params.Ls_start):find_index(freq, params.Ls_end)));
else % Ls_mean given by the user
    Ls_mean                                     = params.Ls;
end

% calculate the load impedance
if length(Ls_mean) > 1
    zaload_mean                                 = Rdc + sqrt(-1)*2*pi*freq.*Ls_mean;
else    
    zaload_mean                                 = Rdc + sqrt(-1)*2*pi*freq*Ls_mean;
end

% knowing the load impedance, calculate the rest of the calibration
% coefficients
% using dispersive Ls measured by LRRM. "A_over_C", "zaload" contain root
% choice. 
C                                               = (B - zxload)./((zxload - A_over_C).*zaload);
A                                               = A_over_C.*C;
D                                               = -1./sqrt(A - B.*C);
Ex                                              = [A.*D C.*D B.*D D];
Sx                                              = ABCD2S_Zc(Ex, 50*ones(numpts,1));
[d, Ey]                                         = DeEmbedd_v2(Ex, Tthru, repmat([1 0 0 1], numpts,1), 'RRR');
Sy                                              = ABCD2S_Zc(Ey, 50*ones(numpts));
Rx                                              = S2R(Sx,1);
Ry                                              = S2R(Sy,1);  
% Using Ls_mean averaged over frequency. "A_over_C", "zaload_mean" contain
% root choice
C_mean                                          = (B - zxload)./((zxload - A_over_C).*zaload_mean);
A_mean                                          = A_over_C.*C_mean;
D_mean                                          = 1./sqrt(A_mean - B.*C_mean);

Ex_mean                                         = [A_mean.*D_mean C_mean.*D_mean B.*D_mean D_mean];
Sx_mean                                         = ABCD2S_Zc(Ex_mean, 50*ones(numpts,1));
[d, Ey_mean]                                    = DeEmbedd_v2(Ex_mean, Tthru, repmat([1 0 0 1], numpts,1), 'RRR');
Sy_mean                                         = ABCD2S_Zc(Ey_mean, 50*ones(numpts));
Rx_mean                                         = S2R(Sx_mean,1);
Ry_mean                                         = S2R(Sy_mean,1); 
% For debug: first solution due to root choice of "A_over_C1", using Ls_mean. Only "zaload_mean" contains root choice.
C1                                              = (B1 - zxload)./((zxload - A_over_C1).*zaload_mean);
A1                                              = A_over_C1.*C1;
D1                                              = -1./sqrt(A1 - B1.*C1);
Ex1                                             = [A1.*D1 C1.*D1 B1.*D1 D1];
Sx1                                             = ABCD2S_Zc(Ex1, 50*ones(numpts,1));
[d, Ey1]                                        = DeEmbedd_v2(Ex1, Tthru, repmat([1 0 0 1], numpts,1), 'RRR');
Sy1                                             = ABCD2S_Zc(Ey1, 50*ones(numpts));

% For debug: second solution due to root choice of "A_over_C2", using Ls_mean. Only "zaload_mean" contains root choice
C2                                              = (B2 - zxload)./((zxload - A_over_C2).*zaload_mean);
A2                                              = A_over_C2.*C2;
D2                                              = -1./sqrt(A2 - B2.*C2);
Ex2                                             = [A2.*D2 C2.*D2 B2.*D2 D2];
Sx2                                             = ABCD2S_Zc(Ex2, 50*ones(numpts,1));
[d, Ey2]                                        = DeEmbedd_v2(Ex2, Tthru, repmat([1 0 0 1], numpts,1), 'RRR');
Sy2                                             = ABCD2S_Zc(Ey2, 50*ones(numpts));

%%
% main output
if strcmp(params.output, 'Ls') % Use the dispersive Ls extracted during LRRM as the load inductance value at each frequency
    output.S1                                   = Sx;
    output.S2                                   = Sy;      
    output.R1                                   = Rx;
    output.R2                                   = Ry;
elseif strcmp(params.output, 'Ls_mean') % use Ls_mean as load inductance value                         
    if isfield(params, 'S1ref') && isfield(params, 'S2ref') % if there is a reference calibration, do the calibration comparison
        S1ref                                   = params.S1ref;
        S2ref                                   = params.S2ref;
        R1ref                                   = S2R(S1ref,1);
        R2ref                                   = S2R(S2ref,1);    
        [d, emax]                               = calcomp_Song_v2(Rx_mean, Ry_mean, R1ref, R2ref);
        % flip the sign of some frequency in LRM cal
        for i = 1:numpts
            if emax(i) > 0.5
                Sx_mean(i,[2 3])                = -Sx_mean(i,[2 3]);
                Sy_mean(i,[2 3])                = -Sy_mean(i,[2 3]);        
            end
        end
        Rx_mean                                 = S2R(Sx_mean, 1);
        Ry_mean                                 = S2R(Sy_mean, 1);    
        [d, emax]                               = calcomp_Song_v2(Rx_mean, Ry_mean, R1ref, R2ref);
        output.emax                             = emax;      
    end
        
    output.S1                                   = Sx_mean;
    output.S2                                   = Sy_mean;    
    output.R1                                   = Rx_mean;
    output.R2                                   = Ry_mean;       
end
    
% Two debug calibrations due to root choice.
output.S11                                      = Sx1;
output.S21                                      = Sy1;
output.R11                                      = S2R(Sx1,1);
output.R21                                      = S2R(Sy1,1);

output.S12                                      = Sx2;
output.S22                                      = Sy2;
output.R12                                      = S2R(Sx2,1);
output.R22                                      = S2R(Sy2,1);

% summarizing outputs
meas.freq                                       = freq;
meas.thru                                       = thru; % raw data without any corrections
meas.short                                      = short;
meas.open                                       = open;
meas.load                                       = load;
meas.thru_sw                                    = thru_sw; % switch-term corrected data
meas.short_sw                                   = short_sw;
meas.open_sw                                    = open_sw;
meas.load_sw                                    = load_sw;
meas.sw                                         = sw; % switch terms
meas.Ls                                         = Ls; % raw load inductance values over frequency measured by the LRRM method.
meas.Ls_mean                                    = Ls_mean; 
meas.z_load_over_open                           = z_load_over_open; % measured impedance ratio
meas.z_short_over_load                          = z_short_over_load;


output.input_params                             = params;
output.meas                                     = meas;