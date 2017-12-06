function Y = sigm_system_fg(x,fmu,x0,x1,T,fT,P)

fsig= x(1);
fbar= x(2);

% unit area under the curve
Y(1) = -1 + fbar*trapz(linspace(x0,x1,100),sigmf(linspace(x0,x1,100),[fsig,fmu])); %-1/fbar + x1 - x0 + log((1+exp(-fsig*(x1-fmu)))/(1+exp(-fsig*(x0-fmu))));

% area before T equal to P
Y(2)= -P + fbar*trapz(linspace(x0,T,100),sigmf(linspace(x0,T,100),[fsig,fmu])); %fbar*(T-x0+log((1+exp(-fsig*(T-fmu)))/(1+exp(-fsig*(x0-fmu))))) - P;

% fixing fmu
Y(3) = fmu - fT;
