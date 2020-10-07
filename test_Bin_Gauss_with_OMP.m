%{
Description: This code compares nmse performance of binary and gaussian
matrices using OMP.

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



%% Code begins
% pBuf = [11,13,17,19]';
pBuf = [17,19,23,29]';


nLoops = length(pBuf);

colors = 'rk';
markers1 = {'x','<','p','+'};
markers2 = {'*','>','s','o'};


if length(markers1) ~= nLoops
    warning('Make same length for pBuf and markers');
end

legendNames = {'Binary','Gaussian'};
legend1 = cell(nLoops*2,1);

nmseFull1= cell(nLoops,1);
nmseFull3= cell(nLoops,1);


tic

for ii = 1:nLoops
    %% Parameters for binary matrix
    
    p = pBuf(ii);
    b = 1;
    r = 1;
    z = 0.5*(p-1);%p^b-17;
    q = z;%no of ones in each column
    
    n = p^2;%no fo columns of A matrix
    k = q;
    
    %% Binary matrix generation using euler square
    Ainit = uob_binarymats(p,b,r,z);
    
    
    %% construct column normalized gaussian matrix
    
    Brandn = colNormalise(randn(size(Ainit)));
    [M,N] = size(Ainit);
    
    L = 1000;
    
    max_sparsity = 0.5*p*(p-1)
    
    nmse11 = zeros(max_sparsity,1);
    nmse33 = zeros(max_sparsity,1);
    
    
    for i= 1:max_sparsity
        
        sparsity = i;
        err11 = 0;
        err33 = 0;
        
        for j=1:L
            x = zeros(N,1);
            indx = randperm(N,sparsity);
            x(indx) = randn(sparsity,1);
            
            y1 = Ainit*x;
            y3 = Brandn*x;
            
            xhat11 = omp(Ainit,y1,sparsity);
            
            xhat33 = omp(Brandn,y3,sparsity);
            
            ratio11 = norm(xhat11-x,2)/norm(x);
            ratio33 = norm(xhat33-x,2)/norm(x);
            
            err11 = err11+ratio11.^2;
            err33 = err33+ratio33.^2;
            
        end
        
        nmse11(i) = err11/L;
        nmse33(i) = err33/L;
        
    end
    
    nmseFull1{ii} = nmse11;
    nmseFull3{ii} = nmse33;
    
    %% save file with date
    time = datestr(now, 'yyyy_mm_dd');
    filename = sprintf('errMatBinaryRandCombined_%s.mat',time);
    save( fullfile('./matfiles', filename), 'ii', 'nmseFull1', 'nmseFull3', 'Brandn');
    
    
    figure(1);
    plot(nmse11,strcat('-',colors(1),markers1{ii}),'MarkerIndices', 1:2:length(nmse33));hold on;
    plot(nmse33,strcat('--',colors(2),markers2{ii}),'MarkerIndices', 1:2:length(nmse33));
    
    legend1{2*(ii-1)+1} = strcat(legendNames{1},'-',num2str(p^2),'x',num2str(2*p^2));
    legend1{2*(ii-1)+2} = strcat(legendNames{2},'-',num2str(p^2),'x',num2str(2*p^2));
    
end
toc

figure(1);
ylabel('NMSE')
xlabel('Sparsity')
title('OMP performance');
ylim([0,0.5])
grid on;
legend(legend1,'Location','Best');







