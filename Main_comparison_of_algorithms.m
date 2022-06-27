clc;
clear all;close all;
global n1 n2 n mk m  S2 q  kk;
%filenames=  {'Pincat.mat','brain_T2_T1.mat','speech_seq.mat','Cardiac_ocmr_data.mat','lowres_speech.mat','FB_ungated.mat'};
%filenames={'brain_T2_T1.mat'};




under_k_space=load('undersampled_kspace.mat')
under_k_space=double(cell2mat(struct2cell(under_k_space)));
[n1,n2,q]=size(under_k_space);
n=n1*n2;
uks=reshape(under_k_space,[n1*n2,q]);

%%%%% Mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask=load('mask.mat');
mask=double(cell2mat(struct2cell(mask)));
mask4 = fftshift(mask);
mask3=reshape(mask,[n1*n2,q]);

mk=[];
for i=1:1:q
    mk(i)=length(find(logical(mask3(:,i))));
    S2(1:mk(i),i)=double(find(logical(mask3(:,i))));
end
m=max(mk);
Y=zeros(m,q);
for k=1:1:q
    ksc=uks(:,k);
    Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ktslr %%%%%%%%%%%%%%%%%%%%%%%%%%
S=find(mask~=0);
Akt = @(z)A_fhp3D(z,S,n1,n2,q);
Aktt=@(z)At_fhp3D(z,S,n1,n2,q);
step_size = [1,1,1];
[D,Dt] = defDDt(step_size);
b = uks(S);
tic;
x_init = Aktt(b);
mu1 =1e-10;
mu2 =4e-9;
opts.mu1 = mu1;
opts.mu2 = mu2;
opts.p=0.1;
[~,sq,~]=givefastSVD(reshape(x_init, n1*n2,q));
opts.beta1=10./max(sq(:));
opts.beta2=10./max(abs(x_init(:)));
opts.beta1rate = 50;
opts.beta2rate = 25;
opts.outer_iter =15;
opts.inner_iter = 50;
[Xhat_ktslr,cost,opts] = minSNandTV(Akt,Aktt,D,Dt,x_init,b,1,opts);
Time_Ktslr=toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L+S-Otazo %%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=find(mask~=0);
A = @(z)A_fhp3D(z, S,n1,n2,q);
At = @(z) At_fhp3D(z, S, n1,n2,q);
param.A=A;
param.At=At;
param.d = uks(S);
param.T=TempFFT(3);
tic;
param.lambda_L=0.01;
param.lambda_S=0.01;
param.nite=50;
param.tol=0.0025;
M=At(param.d);

M=reshape(M,[n1*n2,q]);
Lpre=M;
S=zeros(n1*n2,q);
ite=0;
while(1)
    ite=ite+1;
    M0=M;
    [Ut,St,Vt]=svd(M-S,0);
    St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
    L=Ut*St*Vt';
    S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
    resk=param.A(reshape(L+S,[n1,n2,q]))-param.d;
    M=L+S-reshape(param.At(resk),[n1*n2,q]);
    Lpre=L;
    tmp2=param.T*reshape(S,[n1,n2,q]);
    if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
end
Xhat_LpS1=L+S;
Xhat_LpS=reshape(Xhat_LpS1,[n1,n2,q]);
Time_LSparse= toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% L+S-Lin %%%%%%%%%%%%%%%%%%%%%%%%%%%%


under_k_space_lin=load('undersampled_kspace_lin.mat')
under_k_space_lin=double(cell2mat(struct2cell(under_k_space_lin)));
smap = ones(n1,n2,1);
E=getE_sc(smap,q,'samp',mask);
tic;
d = under_k_space_lin;
param.d = d;
opt.d = param.d;
L = E'*param.d;
res = E*L-param.d;
[~,St,~]=svd(reshape(L,[n,q])-reshape(E'*res,[n,q]),0);

param.E = getE(smap,q,'samp',mask);
param.T = getT(n1,n2,q);
param.nite = 10;
param.scaleL = St(1);
param.scaleS = 1/1.887;
param.lambda_L=0.01;
param.lambda_S=0.05*param.scaleS;

[L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
L = L_pogm;S = S_pogm;
LplusS=L+S;
Time_LplusS_jeff=toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  altGDmin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=max(mk);
Y=Y(1:m,1:q);
T=70;
tic;
[U0]=initAltGDMin(Y);
[Uhat, Bhat]=AltGDmin(T,U0,Y);
X_hat=reshape(Uhat*Bhat,[n1, n2,q]);
Time_GD=  toc;
T=70;
tic;
L=[];
[zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
Ytemp=reshape(Afft(zbar_hat),[m,q]);
Ybar=Y-Ytemp;
[U0]=initAltGDMin(Ybar);
[Uhat, Bhat]=AltGDmin(T,U0,Ybar);
X_hat=Uhat*Bhat;

S=find(mask~=0);
param.Samp_loc=S;
A = @(z)A_fhp3D(z, S,n1,n2,q);
At = @(z) At_fhp3D(z, S, n1,n2,q);
param.A=A;
param.At=At;
param.d = uks(S);
param.T=TempFFT(3);
param.lambda_L=0.01;

param.nite=10;
param.tol=0.0025;
M=At(param.d);
M=reshape(M,[n1*n2,q]);
Lpre=M;
Ehat=zeros(n1*n2,q);
L(:,1:q)=X_hat+zbar_hat;
param.lambda_S=0.001*max(max(abs(M-L)));
ite=0;

while(1)
    ite=ite+1;
    M0=M;
    Ehat=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
    resk=param.A(reshape(L+Ehat,[n1,n2,q]))-param.d;
    M=L+Ehat-reshape(param.At(resk),[n1*n2,q]);
    Lpre=L;
    tmp2=param.T*reshape(Ehat,[n1,n2,q]);
    if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
end
Zhat=L+Ehat;
Zhat_MRI2=reshape(Zhat,n1,n2,q);
Time_MRI2=  toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% altGDmin-MRI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=max(mk);
L=[];
T=70;
tic;


[zbar_hat,flag,resNE,iter] = cgls(@Afft,@Att, Y,0,1e-36,10);
Ytemp=reshape(Afft(zbar_hat),[m,q]);
Ybar=Y-Ytemp;
[U0]=initAltGDMin(Ybar);
[Uhat, Bhat]=AltGDmin(T,U0,Ybar);
X_hat=Uhat*Bhat;

Yhat_hat=Y-Afft(X_hat+zbar_hat);
Ehat=[];
for kk=1:1:q
    Ehat(:,kk)=cgls_modi(@Afft_modi,@At_modi, Yhat_hat(:,kk) ,0,1e-36,3);
end
Zhat=X_hat+zbar_hat+Ehat;
Zhat_MRI=reshape(Zhat,[n1, n2,q]);

Time_MRI=  toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1=max(max(max(abs(Xhat_ktslr))));subplot(1,5,1); imagesc(abs(Xhat_ktslr(60:196,:,20)),[0,0.5*m1]);colormap(gray); axis off;
m2=max(max(max(abs(Xhat_LpS))));subplot(1,5,2); imagesc(abs(Xhat_LpS(60:196,:,20)),[0,0.5*m2]);colormap(gray); axis off;
m3=max(max(max(abs(LplusS))));LplusS1=imrotate(LplusS(:,:,20),-180);subplot(1,5,3); imagesc(abs(LplusS1(60:196,:,1)),[0,0.5*m3]);colormap(gray); axis off;
m4=max(max(max(abs(Zhat_MRI))));subplot(1,5,4); imagesc(abs(Zhat_MRI(60:196,:,20)),[0,0.5*m4]);colormap(gray); axis off;
m5=max(max(max(abs(Zhat_MRI2))));subplot(1,5,5); imagesc(abs(Zhat_MRI2(60:196,:,20)),[0,0.5*m5]);colormap(gray); axis off;

function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end
