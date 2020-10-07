%{
Description: This code compares nmse performance of UNTF and sparse
gaussian matrices using BP.

Measurement matrix with different dimensions are considered.
save nmse for every Measurement matrix

%}
clc,clear;
close all;

addpath('gen_random_untf');
addpath('matrix_generation');
addpath('SRA');


addpath('/home/pradipsasmal/Documents/');
defpaper;% setting for figures


%% Code starts from here
 pBuf = [11,13,17,19]';
%pBuf = [17,19,23,29]';

nLoops = length(pBuf);

colors = 'rkb';
markers1 = {'x','<','p','+'};
markers2 = {'*','>','s','o'};
markers3 = {'*','>','s','o'};



if length(markers1) ~= nLoops
    warning('Make same length for pBuf and markers');
end

legendNames = {'Embedd\_UNTF','Sparse\_Gaussian','Gaussian'};
legend1 = cell(nLoops*3,1);


nmseFull1= cell(nLoops,1);
nmseFull2= cell(nLoops,1);
nmseFull3= cell(nLoops,1);

flg_run = 1;% 1- runs one algorithm at a time otherwise run all the three algorithms
tic

for ii = 1:nLoops
    %% Parameters for UNTF
    ii
    p = pBuf(ii);
    %p = 19;
    b = 1;
    r = 1;
    %z = 0.5*(p-1);%p^b-17;
    z = 2;
    q = z;%no of ones in each column
    
    n = p^2;%no fo columns of A matrix
    k = q;
    %% UNTF matrix generation
    Ainit = uob_binarymats(p,b,r,z);
    flg_untf = 2; % 1-dft,2-dct,other for hadamard
    
    Auntf = untf_mod(Ainit',flg_untf);
    
    %% Construct sparse random measurement matrix with normalized columns
    [M,N] = size(Auntf);
    BrandnSp= sparse_gau_mat_gen(M,N,p);
    
    %% Gaussian Matrix generation
    
    Brandn = colNormalise(randn(size(Auntf)));
    
    L = 1000;
    
    max_sparsity = 0.5*p*(p-1)
    
    nmse1 = zeros(max_sparsity,1);
    nmse2 = zeros(max_sparsity,1);
    nmse3 = zeros(max_sparsity,1);
    %# for l1 recovery initialisation
    pinvAuntf = pinv(Auntf);
    pinvBsp= pinv(BrandnSp);
    pinvB  = pinv(Brandn);
    
    for i= 1:max_sparsity
        
        sparsity = i;
        
        err1 = 0;
        err2 = 0;
        err3 = 0;
       
        for j=1:L
            
            %# Sparse vector generation
            x = zeros(N,1);
            indx = randperm(N,sparsity);
            x(indx) = randn(sparsity,1);
            
            %# Compressed measurements
            y1 = Auntf*x;
            y2 = BrandnSp*x;
            y3 = Brandn*x;
            
            %# Sparse recovery
            xhat1 = l1eq_pd(pinvAuntf*y1,Auntf,[],y1);
            xhat2 = l1eq_pd(pinvBsp*y2,BrandnSp,[],y2);
            xhat3 = l1eq_pd(pinvB*y3,Brandn,[],y3);
            
            %# nmse computation
            ratio1 = norm(xhat1-x,2)/norm(x);
            ratio2 = norm(xhat2-x,2)/norm(x);
            ratio3 = norm(xhat3-x,2)/norm(x);
            
            err1 = err1+ratio1.^2;
            err2 = err2+ratio2.^2;
            err3 = err3+ratio3.^2;
            
        end
        
        %# nmse computation
        nmse1(i) = err1/L;
        nmse2(i) = err2/L;
        nmse3(i) = err3/L;
        
    end
    toc
    nmseFull1{ii} = nmse1;
    nmseFull2{ii} = nmse2;
    nmseFull3{ii} = nmse3;
    
    %% save file with date
    time = datestr(now, 'yyyy_mm_dd');
    filename = sprintf('errMatsparseRand_%s.mat',time);
    save( fullfile('./matfiles', filename), 'ii', 'nmseFull1','nmseFull2','nmseFull3');
    
    
    figure(1);
    plot(nmse1,strcat('-',colors(1),markers1{ii}),'MarkerIndices', 1:2:length(nmse1));hold on;
    plot(nmse2,strcat('-',colors(2),markers2{ii}),'MarkerIndices', 1:2:length(nmse2));
    plot(nmse3,strcat('--',colors(3),markers3{ii}),'MarkerIndices', 1:2:length(nmse3));
    
    legend1{3*(ii-1)+1} = strcat(legendNames{1},'-',num2str(p^2),'x',num2str(2*p^2));
    legend1{3*(ii-1)+2} = strcat(legendNames{2},'-',num2str(p^2),'x',num2str(2*p^2));
    legend1{3*(ii-1)+3} = strcat(legendNames{3},'-',num2str(p^2),'x',num2str(2*p^2));

end
toc
figure(1);
ylabel('NMSE')
xlabel('Sparsity')
title('BP performance');
ylim([0,0.5])
grid on;
legend(legend1,'Location','Best');







