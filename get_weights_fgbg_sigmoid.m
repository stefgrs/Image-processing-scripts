function [w_t_fg,w_t_bg,sigma] = get_weights_fgbg_sigmoid(I,T,fT,gT,commA)

% here I assume sigmoidal probability distributions for foreground and
% background pixels, f(x), g(x) respectively, such that:
% i. each of them is given by a sigmoidal function: 
% f(x) = fbar*1/(1+exp(-fsig*(x-fmu)))
% g(x) = gbar*1/(1+exp(-gsig*(x-gmu)))
% ii. they go from x0= min(I(:)) to x1= max(I(:))
% iii. they intersect at T: f(T) = g(T)
% iv. fmu= T1
% v. gmu= T2
% area(f)= 1;
% area(g)= 1;
% area(f < T) = commA/2
% area(g > T) = commA/2
% the parameters are given by solving a system of nonlinear equations

Is=I;
x0= min(Is(~isnan(Is(:))));
x1= max(Is(~isnan(Is(:))));

options = optimoptions('fsolve','Display','off','Algorithm','levenberg-marquardt');
xf= fsolve(@(x)sigm_system_fg(x,fT,x0,x1,T,fT,commA/2),[10,1],options);
xg= fsolve(@(x)sigm_system_bg(x,gT,x0,x1,T,gT,commA/2),[-10,1],options);

fsig= xf(1);
fbar= xf(2);
gsig= xg(1);
gbar= xg(2);

I(isnan(Is(:)))= x0;

% calculate weights
% foreground: first probability, then turn into costs
p_fg= fbar*sigmf(I(:),[fsig,fT]);
w_t_fg= -log(abs(p_fg));
% background
p_bg= gbar*sigmf(I(:),[gsig,gT]);
w_t_bg= -log(abs(p_bg));

% compute sigma, that is scale of separation
sigma= std(I(:));
