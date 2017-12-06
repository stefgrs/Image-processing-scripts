function Y = sigm_system_bg(x,gmu,x0,x1,T,gT,P)

gsig= x(1);
gbar= x(2);

% unit area under the curve
Y(1) = -1 + gbar*trapz(linspace(x0,x1,100),sigmf(linspace(x0,x1,100),[gsig,gmu])); %-1/fbar + x1 - x0 + log((1+exp(-fsig*(x1-fmu)))/(1+exp(-fsig*(x0-fmu))));

% area before T equal to P
Y(2)= -P + gbar*trapz(linspace(T,x1,100),sigmf(linspace(T,x1,100),[gsig,gmu])); %fbar*(T-x0+log((1+exp(-fsig*(T-fmu)))/(1+exp(-fsig*(x0-fmu))))) - P;

% fixing fmu (??? better constraints ???)
Y(3) = gmu - gT;
