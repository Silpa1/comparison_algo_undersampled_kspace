% L+S reconstruction of undersampled multicoil PINCAT phantom
% Claire Lin, 06/05/2018
clear all;close all;
filenames={'brain_T2_T1.mat'}
%new=load('FB_ungated.mat')
%filenames=  {'Pincat.mat','brain_T2_T1.mat','speech_seq.mat','Cardiac_ocmr_data.mat','FB_ungated.mat'}%,'FB_ungated.mat',,'image_xyt'}  {'Cardiac_ocmr_data.mat'}%; %Cardiac'FB_ungated.mat',
%filenames=  {'lowres_speech.mat'};%,'lowres_speech.mat','FB_ungated.mat'}%,'FB_ungated.mat',,'image_xyt'}  {'Cardiac_ocmr_data.mat'}%; %Cardiac'FB_ungated.mat',
fid = fopen('Comparison_Jeff.txt','wt');
fprintf(fid, '%s(%s) & %s    \n','Dataset','Radial','LplusS_jeff');

for jj = 1:1:numel(filenames)
    S = load(filenames{jj});
    new=cell2mat(struct2cell(S));
    [~,name,~] = fileparts(filenames{jj});% Best to load into an output variable.
    
    %load('Xinf.mat')
    tmp = max(new(:));
    Xtrue = div0(new,tmp);
    [nx,ny,nt] = size(Xtrue);
    %im(Xtrue)
    
    
    smap = ones(nx,ny,1);
    
    radialline =[4];
    for ii=1:1:length(radialline)
        tic;
        
        samp = goldencart(nx,ny,nt,radialline(ii));
        %im(samp)
        %% data
        E=getE_sc(smap,nt,'samp',samp);
        d = E*Xtrue;
        % add noise
        rng(0)
        dn = randn(size(d)) + 1i * randn(size(d));
        param.snr_db =23;
        param.scale_noise = norm(d(:)) / norm(dn(:)) / 10.^(param.snr_db / 20);
        param.d = d;% + param.scale_noise * dn;
        %printm('data rmse = %g, snr = %g', rms(param.d(:)-d(:)), ...
            %20*log10(norm(d(:)) / norm(param.d(:)-d(:))))
        opt.d = param.d;
        % prepare for regularization scaling
        L = E'*param.d;
        res = E*L-param.d;
        [~,St,~]=svd(reshape(L,[nx*ny,nt])-reshape(E'*res,[nx*ny,nt]),0);
        %semilogy(diag(St),'ok')
        %% prepare for AL: opt
       
        %% prepare for PGM: param
        param.E = getE(smap,nt,'samp',samp);
        param.T = getT(nx,ny,nt);
        param.nite = 10;
        param.scaleL = St(1);
        param.scaleS = 1/1.887; %1 / b1 constant squared in middle of image
        param.lambda_L=0.01;
        param.lambda_S=0.05*param.scaleS;
        %param.Xinf = reshape(Xinf.pincat,nx*ny,nt);
        %% ISTA
        [L_ista,S_ista,xdiff_ista,cost_ista,time_ista,rankL_ista] = PGM(param);
        %% FISTA
        [L_fista,S_fista,xdiff_fista,cost_fista,time_fista,rankL_fista] = PGM(param,'fistaL',1,'fistaS',1);
        %% POGMmask2
        [L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
        %% Display: 1 frame
        %figure; tpick = 21;
        L = L_pogm;S = S_pogm;
        LplusS=L+S;
     Time_LplusS_jeff=toc;
        Error_LplusS_jeff=RMSE_modi(LplusS,Xtrue);
        similarity_index=[];
        for i =1:1:nt
            mssim=ssim(abs(LplusS(:,:,i)/max(max(LplusS(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
            similarity_index(i)=mssim;
        end
        save('C:\Users\sbabu\Desktop\Results\low_res_8\Xhat_LpS_lin.mat', 'LplusS');
        sim_LplusS_jeff=min(similarity_index);
        %filename(ii)=name;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fid, '%s(%d) & %8.4f (%5.2f,%8.4f) \n', name, radialline(ii),Error_LplusS_jeff,Time_LplusS_jeff,sim_LplusS_jeff);

    end
end
fclose(fid);
        imshow(abs(Xtrue(:,:,tpick)),[0,1])
        title({'$X_{true}$'},'Interpreter','latex')
        subplot(1,5,2)
        imshow(abs(LplusS(:,:,tpick)),[0,1])
        title({'$\tilde{X}$'},'Interpreter','latex')
        subplot(1,5,3)
        imshow(abs(L(:,:,tpick)),[0,1])
        title({'$L_\infty$'},'Interpreter','latex')
        subplot(1,5,4)
        imshow(abs(S(:,:,tpick)),[0,1])
        title({'$S_\infty$'},'Interpreter','latex')
        subplot(1,5,5)
        imagesc(abs(LplusS(:,:,21))-abs(Xtrue(:,:,tpick)),[0,0.2])
        title({'$|X_{true}-\tilde{X}|$'},'Interpreter','latex')
        axis off;axis square;