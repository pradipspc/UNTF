%{
Construct sparse random measurement matrix with normalized columns
No of non zeros per column remains same but support is different.

%}
function aSparseMat = sparse_gau_mat_gen(nrows,ncols,nNonZeros_per_col)
%     [nrows,ncols] = size(Auntf);
%     nNonZeros_per_col = p;
    aSparseMat = zeros(nrows,ncols);
    
    for i=1:ncols
        temp = zeros(nrows,1);
        indx = randperm(nrows,nNonZeros_per_col);
        temp(indx) = randn(nNonZeros_per_col,1);
        
        aSparseMat(:,i) = temp;
    end
    
    aSparseMat = normc(aSparseMat);
end