function [output] = SORefMat_cal(varargin)

% A general LRM calibration.
% treats unequal but known load impedances. 

files                                           = varargin{1,1};
params                                          = varargin{1,2};

path                                            = files.path;
% load in data
[short, freq]                                   = readTouchStone_s1p([path files.short]);
[open,d]                                        = readTouchStone_s1p([path files.open]);
[load,d]                                        = readTouchStone_s1p([path files.load]);

plotSparS11(sqrt(real(short).^2+imag(short).^2), freq)
plotSparS11(sqrt(real(open).^2+imag(open).^2), freq)
plotSparS11(sqrt(real(load).^2+imag(load).^2), freq)
% calculate the load impedance of the Open and Reference material standards

D                                               = params.D;                                               
d                                               = params.d;                                               
er_sub                                          = params.er_sub; 
c0                                              = params.c0;
cf                                              = params.cf;

Z0                                              = 138*log10(D/d)/sqrt(er_sub);
Er_open                                         = ones(length(freq),1) - sqrt(-1)*zeros(length(freq),1);
Er_ellison                                      = water_Er_Ellison(20, freq);
x                                               = [c0 cf Z0];
[Z_L2,Gamma_L2]                                 = E_com2Load(x,freq,Er_open);
[Z_L3,Gamma_L3]                                 = E_com2Load(x,freq,Er_ellison);

gamma_2                                         = Gamma_L2;
gamma_3                                         = Gamma_L3;

% 
Rho1                                            = short;
Rho2                                            = open;
Rho3                                            = load;

% defination 
Gamma1                                          = -1;

numpts                                          = length(freq);
S11                                             = zeros(numpts, 1);
S22                                             = zeros(numpts, 1);
S12S21                                          = zeros(numpts, 1);
for i = 1:numpts
    Gamma2                                      = gamma_2(i,1);
    Gamma3                                      = gamma_3(i,1);
    
    T12                                         = Gamma1*Gamma2*(Rho1(i)-Rho2(i));
    T23                                         = Gamma2*Gamma3*(Rho2(i)-Rho3(i));
    T31                                         = Gamma3*Gamma1*(Rho3(i)-Rho1(i));
    
    s11                                         = (Rho3(i)*T12 + Rho2(i)*T31 + Rho1(i)*T23)/(T12 + T31 + T23);
    s22                                         = (Gamma1*(Rho2(i) - s11) + Gamma2*(s11 - Rho1(i)))/(Gamma1*Gamma2*(Rho2(i) - Rho1(i)));
    s12s21                                      = (Rho1(i) - s11)*(1 - s22*Gamma1)/Gamma1;
    
    S11(i,1)                                    = s11;
    S22(i,1)                                    = s22;
    S12S21(i,1)                                 = s12s21;
end

output.S11                                      = S11;
output.S22                                      = S22;
output.S12S21                                   = S12S21;

end

