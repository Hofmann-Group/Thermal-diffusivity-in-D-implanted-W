function y=decay_inc_ampl(X,x)
% This function is called by leasqr
% x is a vector which contains the coefficients of the
% equation. X and Y are the option data sets that were
% passed to leasqr.
% A-peak height  B-centre  C-width  D-background

q = 2*pi/4* 10^-6;

Ap=x(1); % amplitude of phase grating decay
Bp=x(2); % thermal diffusivity
alpha = Bp^2/q^2;
%alpha=x(2);
Cp=x(3); % amplitude of saw signal
Dp=x(4); % frequency of saw signal
Ep=x(5); % phase of saw signal
Fp=x(6); % time constant of saw signal decay
Gp=x(7); % background
%Hp=x(8); % amplitude of amplitude grating

% mar
%ph_grating  = Ap.*erfc(q.*sqrt(alpha).*sqrt(X));
ph_grating  = Ap.*erfc(Bp.*sqrt(X));

saw         = Cp.*cos((Dp*2*pi).*X + Ep).*exp(-X./Fp);
%amp_grating = Hp.*(alpha.*X).^(-0.5).*exp(-alpha.*q^2.*X);
bkg = Gp;

%y=ph_grating - amp_grating + saw + bkg;
y=ph_grating  + saw + bkg;