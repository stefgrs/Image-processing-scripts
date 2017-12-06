function [Th,nEl,cent]= T_point_thresh(Data,Bins)

% input
% D: data from which you want to compute the threshold, given as a vector
% Bins: number of bins for the histogram

[M,L,nEl,cent]= get_descend_bins(Data,Bins);

n= 0;
Sh=0;
Shh= 0;
Sg= 0;
Sgg= 0;
Sgh= 0;
eps1= zeros(1,L-M);
% for kk=L, it would be like having one slope only and epsilon2 would be
% undefined, therefore it's better to stop at L-1, also because it's
% practically impossible that having only one slope minimizes the error.
for kk=M:L-1
    n= n+1;
    Sh= Sh + nEl(kk);
    Shh= Shh + nEl(kk).^2;
    Sg= Sg + cent(kk);
    Sgg= Sgg + cent(kk).^2;
    Sgh= Sgh + nEl(kk).*cent(kk);
    eps1(kk-M+1) = Shh - Sh^2/n - ((n*Sgh - Sg * Sh)^2)/(n * (n*Sgg-Sg^2));
end
n= 0;
Sh=0;
Shh= 0;
Sg= 0;
Sgg= 0;
Sgh= 0;
eps2= zeros(1,L-M);
% imagine that it's again kk=L-1:-1:M, then change variable kk=kk+1
for kk=L:-1:(M+1)
    n=n+1;
    Sh= Sh + nEl(kk);
    Shh= Shh + nEl(kk).^2;
    Sg= Sg + cent(kk);
    Sgg= Sgg + cent(kk).^2;
    Sgh= Sgh + nEl(kk).*cent(kk);
    eps2(L-M-n+1) = Shh - Sh^2/n - ((n*Sgh - Sg*Sh)^2)/(n * (n*Sgg-Sg^2));
end

eps_tot= eps1 + eps2;
[~,idx]= min(eps_tot);

idx= idx+M-1;
if idx<=M
    idx= M+1;
end
Th= cent(idx);
    

