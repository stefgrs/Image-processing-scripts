function [Bw,flow]= graph_cut_segV3(I,modS)
% description: binary segmentation of image I with graph cut algorithm
% See book "Optimization for computer vision", chapter 6
% Input:
% I: 2D-image to segment
% modS: multiplier for the parameter sigma. The higher the multiplier, the
% less pixels will be assigned to the foreground (generally)
% Output:
% Bw= binary image. Pixels with label 0 are background and pixels with
% label 1 are foreground
% flow= the calculated maximum flow value

% if isempty(Ishow) %nargin<3
%     Ishow= I;
% end

if isempty(modS)
    modS= 1;
end

% number of nodes, size and list of the indices
N= numel(I);
Sz= size(I);
idxAll=1:N;

% Terminal weights. Note: here they indicates the cost of assigning each
% node to each label (bg or fg), as required in the 'maxflow' function,
% and not the weight of the edge that connect the node to either the source
% or the sink node as described in the book
sig= 1;
Is= I; %imfilter(I,fspecial('gaussian',[4*ceil(sig),4*ceil(sig)],sig),'symmetric');
% Assume linear probability distribution and estimate weights accordingly
% Th: threshold between foreground and background --> unimodal threshold?
% commA: percentage of overlap between fg and bg
% Th and commA determine how "conservative" I want to be (CHECK!, ?????)
Th= modS*mean(Is(and(Is(:)>0,~isnan(Is(:)))));
commA= 2*.05;
[w_t_fg,w_t_bg,sigmaI]= get_weights_fgbg_sigmoid(Is(:),Th,Th*1.1,Th*0.9,commA);
sigma= max(sigmaI/10,min(sigmaI,Th*(1.1-1)));
%[w_t_fg,w_t_bg,sigma]= get_weights_fgbg_lin(Is(:),Th,commA);
% obs: the terminal weights must be higher enough with respect to the 
% internal ones such that they do NOT form the mincut.
T= [w_t_bg,w_t_fg];

I(isnan(I))= min(I(~isnan(I(:))));
% Create graph adiacency matrix, 4 neighbours connections (A basically represents
% a rectangular grid)
% A= sparse(N,N);
% connection between i and i+1 (for not-border nodes), so the connection to
% the node immediately below
idxD=(idxAll+1).*repmat([ones(1,Sz(1)-1), 0],[1,Sz(2)]);
% connection between i and i-1 (for not-border nodes), so the connection to
% the node immediately above
idxU=(idxAll-1).*repmat([0,ones(1,Sz(1)-1)],[1,Sz(2)]);
% connection to the node on the right, which means i --> i+Sz(1) and every
% node from 1 to N-Sz(1) has it
idxR= idxAll(1:end-Sz(1))+Sz(1);
% connection to the node on the left, which means i --> i-Sz(1) and every
% node from Sz(1)+1 to N has it
idxL= idxAll(Sz(1)+1:end)-Sz(1);
% get the proper weights
coeff= 1;
w_nD= (exp(-((I(idxAll(idxD>0))-I(idxD(idxD>0))).^2)/2/(coeff*sigma)^2));
w_nU= (exp(-((I(idxAll(idxU>0))-I(idxU(idxU>0))).^2)/2/(coeff*sigma)^2));
w_nR= (exp(-((I(idxAll(1:end-Sz(1)))-I(idxR)).^2)/2/(coeff*.5*sigma)^2));
w_nL= (exp(-((I(idxAll(Sz(1)+1:end))-I(idxL)).^2)/2/(coeff*.5*sigma)^2));
% directed weights
% w_nD= (exp(-((I(idxAll(idxD>0))-I(idxD(idxD>0))))/2/sigma));
% w_nU= (exp(-((I(idxAll(idxU>0))-I(idxU(idxU>0))))/2/sigma));
% w_nR= (exp(-((I(idxAll(1:end-Sz(1)))-I(idxR)))/2/sigma));
% w_nL= (exp(-((I(idxAll(Sz(1)+1:end))-I(idxL)))/2/sigma));
% Create the sparse matrix
A= sparse(idxAll(idxD>0),idxD(idxD>0),w_nD,N,N);
A= A + sparse(idxAll(idxU>0),idxU(idxU>0),w_nU,N,N); % I can sum since the indices are not overlapping
A= A + sparse(idxAll(1:end-Sz(1)),idxR,w_nR,N,N); 
A= A + sparse(idxAll(Sz(1)+1:end),idxL,w_nL,N,N); 
%A(sub2ind([N,N],idxAll,idxU))= w_nU;


% % to plot A with the right coordinates
% coord(:,1)= reshape(repmat(1:Sz(2),[Sz(1),1]),prod(Sz),[]);
% coord(:,2)=repmat((Sz(1):-1:1)',[Sz(2),1]);
% coord(prod(Sz)+1,:)= [(Sz(2)+1)/2, Sz(1)+1];
% coord(prod(Sz)+2,:)= [(Sz(2)+1)/2, 0];
% figure
% gplot(A,coord)


% Compute maxflow
[flow,labels] = maxflow(A,sparse(T));

% Reshape labels to have binary mask
Bw=reshape(~labels,Sz);

% % partition stability algorithm
% T = 10.^[-1:.1:1];
% [S, N, VI, C] = stability(A, T,'noVI');

% %
% Bwbg= 0;
% Bwfg= 0;


% Loose Description:
% Basically, for each intensity value in the image there are two associated 
% probability values that indicate how likely it is for that intensity value 
% to belong to a foreground pixel or to a background pixel. Generally, the 
% higher the intensity value, the higher the foreground probability and the
% lower the background one. I modelled these probability distributions as 
% sigmoid functions (increasing for the foreground probability and 
% decreasing for the background), that intersect at a threshold T, reach 
% the half activation point at certain thresholds (two different ones, 
% given by 1.1*T and 0.9*T respectively) and have a certain percentage of 
% area in common (the parameter commA in the graph cut function). I opted 
% for the sigmoid functions having in mind some sort of “activation process” 
% (in a very loose sense of the word). The parameter T is given by the 
% average of the positive intensity values of the image multiplied by the
% modifier modS: the higher modS is, the less probable it is for even high 
% intensity values to be associated with the foreground, therefore the 
% algorithm performs more conservative choices and more pixels are assigned
% to the background. The default value for modS is 1 and, for now, I never 
% needed large deviation from unity to achieve good results (good at least 
% in the sense of simple visual inspection).