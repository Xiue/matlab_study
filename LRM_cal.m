function [output] = LRM_cal(varargin)

% A general LRM calibration.
% treats unequal but known load impedances. 

files                                           = varargin{1,1};
params                                          = varargin{1,2};

path                                            = files.path;
% load in data
[thru, freq]                                    = readTouchStone_statistiCal([path files.thru]);
[reflect,d]                                     = readTouchStone_statistiCal([path files.reflect]);
[load,d]                                        = readTouchStone_statistiCal([path files.load]);
numpts                                          = length(freq);

% calculate the load impedances
% if the user gives Rdc 
if isfield(params, 'Rdc1')
    if length(params.Rdc1) == 1
        Rdc1                                        = params.Rdc1*ones(numpts,1);
    else
        Rdc1                                        = params.Rdc1;
    end
end

if isfield(params, 'Rdc2')
    if length(params.Rdc2) == 1
        Rdc2                                        = params.Rdc2*ones(numpts,1);
    else
        Rdc2                                        = params.Rdc2;
    end
end

% if the user gives load inductance
if isfield(params, 'Ls1')
    if length(params.Ls1) == 1
        Ls1                                         = params.Ls1*ones(numpts,1);
    else
        Ls1                                         = params.Ls1;
    end
end

if isfield(params, 'Ls2')
    if length(params.Ls2) == 1
        Ls2                                         = params.Ls2*ones(numpts,1);
    else
        Ls2                                         = params.Ls2;
    end
end

% If the user directly gives load impedances
if isfield(params, 'Zx')
    yax                                         = 1./params.Zx;
else
    yax                                         = 1./(Rdc1 + sqrt(-1)*2*pi*freq.*Ls1);
end

if isfield(params, 'Zy')
    yay                                         = 1./params.Zy;
else
    yay                                         = 1./(Rdc2 + sqrt(-1)*2*pi*freq.*Ls2);
end

if isfield(files, 'sw')
    [sw, d]                                     = readTouchStone_statistiCal([path files.sw]);
else
    fprintf('no switch-terms are supplied. \n');
    sw                                          = zeros(numpts,4);
end
    
% correct for switch-terms
thru_sw                                         = correct_switchTerms(thru, sw);
reflect_sw                                      = correct_switchTerms(reflect, sw);
load_sw                                         = correct_switchTerms(load, sw);

Tthru                                           = S2ABCD_Zc_v2(thru_sw, 50*ones(numpts,1));

zreflect                                        = S2Z(reflect_sw);
zload                                           = S2Z(load_sw);
% thru measurement
t1                                              = Tthru(:,1); % A
t2                                              = Tthru(:,3); % B
t3                                              = Tthru(:,2); % C
t4                                              = Tthru(:,4); % D

zxreflect                                       = squeeze(zreflect(1,1,:));
zyreflect                                       = squeeze(zreflect(2,2,:));

% measured load impedances at the VNA measurement plane
zmx                                             = squeeze(zload(1,1,:));
zmy                                             = squeeze(zload(2,2,:));
% solve B_X
% measurement of the loads gives: 
% m1*A_X + m2*B_X + m3*C_X + m4 = 0 and
% n1*A_X + n2*B_X + n3*C_X + n4 = 0.
m1                                              = 1;
m2                                              = yax;
m3                                              = -zmx;
m4                                              = -yax.*zmx;

n1                                              = -zmy.*t3 + t4;
n2                                              = yay.*zmy.*t3 - yay.*t4;
n3                                              = zmy.*t1 - t2;
n4                                              = -yay.*zmy.*t1 + yay.*t2;
% equation A_X / C_X = k * (o1*B_X + o2)/(q1*B_X + q2)
k                                               = (m3.*n1 - n3.*m1)./(m1.*n3 - n1.*m3);
o1                                              = n2.*m3 - m2.*n3;
o2                                              = n4.*m3 - m4.*n3;
q1                                              = n2.*m1 - m2.*n1;
q2                                              = n4.*m1 - m4.*n1;
%eq(24) in Leonard 2006 ARFTG
% a1, a2, and Va using a reflect
a1                                              = t1.*zyreflect - t2 + t3.*zxreflect.*zyreflect - t4.*zxreflect;
a2                                              = 2*t4 - 2*t3.*zyreflect;
Va                                              = 2*t1.*zxreflect.*zyreflect - 2*t2.*zxreflect;
% equation A*B_X^2 + B*B_X + C = 0
A                                               = a1.*q1 + a2.*k.*o1;
B                                               = a1.*q2 + a2.*k.*o2 - Va.*q1 + a1.*k.*o1;
C                                               = -Va.*q2 + a1.*k.*o2;
B_X1                                            = (-B + sqrt(B.^2 - 4*A.*C))./(2*A);
B_X2                                            = (-B - sqrt(B.^2 - 4*A.*C))./(2*A);

A_X1                                            = ((n2.*m3 - m2.*n3).*B_X1 + n4.*m3 - m4.*n3)./(m1.*n3 - n1.*m3);
A_X2                                            = ((n2.*m3 - m2.*n3).*B_X2 + n4.*m3 - m4.*n3)./(m1.*n3 - n1.*m3);

C_X1                                            = ((n2.*m1 - m2.*n1).*B_X1 + n4.*m1 - m4.*n1)./(m3.*n1 - n3.*m1);
C_X2                                            = ((n2.*m1 - m2.*n1).*B_X2 + n4.*m1 - m4.*n1)./(m3.*n1 - n3.*m1);
  
yaref1                                          = C_X1.*( (zxreflect - A_X1./C_X1)./(B_X1 - zxreflect) );
yaref2                                          = C_X2.*( (zxreflect - A_X2./C_X2)./(B_X2 - zxreflect) );

zaref1                                          = 1./yaref1;
zaref2                                          = 1./yaref2;
%%
for i = 1:numpts
    if strcmp(params.reflectType, 'Open') % negative capacitance for open
        
        if isfield(params, 'reflectGuess')
            copen_guess                         = params.reflectGuess;
        else
            copen_guess                         = 0;
        end
        
        if isfield(params, 'deltaOpenc') % tuning parameter in making the root choice using the capacitance
            delta_openc                         = params.deltaOpenc;
        else
            delta_openc                         = 2*1e-15;
        end
        
        copen1(i,1)                             = imag(yaref1(i))/(2*pi*freq(i));
        copen2(i,1)                             = imag(yaref2(i))/(2*pi*freq(i));        

        if i == 1 
            if abs(copen1(i) - copen_guess) < abs(copen2(i) - copen_guess)                
                select                          = 1;
            else 
                select                          = 2;
            end
        else % based on observations, capacitance can only decrease. select the one which decreases slower
            deltac1                             = copen1(i) - copen(i-1);
            deltac2                             = copen2(i) - copen(i-1);
            
            if abs(deltac1) <= delta_openc  && abs(deltac2) >= delta_openc                 
                select                          = 1;
            elseif  abs(deltac1) <= delta_openc  &&  abs(deltac2) < delta_openc 
                if deltac1 < 0 && deltac2 > 0
                    select                      = 1;
                elseif deltac1 < 0 && deltac2 <= 0 
                    if abs(deltac1) < abs(deltac2)
                        select                  = 1;
                    else
                        select                  = 2;
                    end                                                                
                elseif deltac1 >= 0 && deltac2 < 0
                    select                      = 2;
                elseif deltac1 >= 0 && deltac2 >= 0
                    if abs(deltac1) <= abs(deltac2)
                        select                  = 1;
                    else
                        select                  = 2;
                    end                                        
                else
                    select                      = 2;
                end
            elseif abs(deltac1) > delta_openc && abs(deltac2) <= delta_openc                
                select                          = 2;
            elseif abs(deltac1) > delta_openc && abs(deltac2) > delta_openc
                if abs(deltac1) < abs(deltac2)
                    select                      = 1;
                else
                    select                      = 2;
                end                                
            end
        end
                         
        if select == 1                  
            A_X(i,1)                            = A_X1(i,1);
            B_X(i,1)                            = B_X1(i,1);
            C_X(i,1)                            = C_X1(i,1);
            yaref(i,1)                          = yaref1(i,1);
            copen(i,1)                          = copen1(i,1);                    
        elseif select == 2
            A_X(i,1)                            = A_X2(i,1);
            B_X(i,1)                            = B_X2(i,1);
            C_X(i,1)                            = C_X2(i,1);
            yaref(i,1)                          = yaref2(i,1);
            copen(i,1)                          = copen2(i,1);                 
        end                                
    elseif strcmp(params.reflectType, 'Short') % negative inductance for short
        if imag(zaref1(i)) < 0
            A_X(i,1)                            = A_X1(i,1);
            B_X(i,1)                            = B_X1(i,1);
            C_X(i,1)                            = C_X1(i,1);
            yaref(i,1)                          = yaref1(i,1);            
        else
            A_X(i,1)                            = A_X2(i,1);
            B_X(i,1)                            = B_X2(i,1);
            C_X(i,1)                            = C_X2(i,1);
            yaref(i,1)                          = yaref2(i,1);            
        end
    else
            disp('please specify reflect type...');
            A_X(i,1)                            = A_X1(i,1);
            B_X(i,1)                            = B_X1(i,1);
            C_X(i,1)                            = C_X1(i,1);
            yaref(i,1)                          = yaref1(i,1);            
    end        
end

if isfield(params, 'debug')
    if strcmp(params.debug, 'on')
        figure1                                     = figure('Position',[100 100 800 800]);
        if strcmp(params.reflectType, 'Open')        
            plot(freq, imag(yaref1)./(2*pi*freq), 'blue');
            hold on;
            plot(freq, imag(yaref2)./(2*pi*freq), 'red');
            plot(freq, imag(yaref)./(2*pi*freq), 'greeno');
            axis([0 110e9 -150e-15 10e-15]);
        end
    end
    
end
%% 
D_X                                             = 1./sqrt(A_X - B_X.*C_X);
Ex                                              = [A_X.*D_X C_X.*D_X B_X.*D_X D_X];
Sx                                              = ABCD2S_Zc(Ex, 50*ones(numpts,1));
[d, Ey]                                         = DeEmbedd_v2(Ex, Tthru, repmat([1 0 0 1], numpts,1), 'RRR');
Sy                                              = ABCD2S_Zc(Ey, 50*ones(numpts));

Rx                                              = S2R(Sx,1);
Ry                                              = S2R(Sy,1);

if isfield(params, 'S1ref') && isfield(params, 'S2ref')
    S1ref                                       = params.S1ref;
    S2ref                                       = params.S2ref;
    
    R1ref                                       = S2R(S1ref,1);
    R2ref                                       = S2R(S2ref,1);    
    % calibration comparison with reference calibrations
    [d, emax0]                                  = calcomp_Song_v2(Rx, Ry, R1ref, R2ref);
    flip_table                                  = zeros(numpts,1);
    % flip the sign of some frequency in LRM cal
    for i = 1:numpts
        if emax0(i) > 0.5
            Sx(i,[2 3])                         = -1*Sx(i,[2 3]);
            Sy(i,[2 3])                         = -1*Sy(i,[2 3]);   
            flip_table(i,1)                     = 1;
        end
    end
    Rx                                          = S2R(Sx, 1);
    Ry                                          = S2R(Sy, 1);    
    [d, emax1]                                  = calcomp_Song_v2(Rx, Ry, R1ref, R2ref);
    % save the calibration comparison result
    output.emax                                 = emax1;
end

output.S1                                       = Sx;
output.S2                                       = Sy;
output.R1                                       = Rx;
output.R2                                       = Ry;