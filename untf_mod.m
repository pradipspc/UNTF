

function B = untf_mod(A,flg)
%{
A is binary matrix with constant number of ones in each column and rows

%}
%this function generates untf from A

% default is dftmatrix
if nargin ~=2
    flg=1;
end

[m,n] = size(A);
s = nnz(A(:,1));

if flg == 1 % normalised dft matrix
    pdtmat = (1/sqrt(s))*fft(eye(s));
elseif flg == 2 % normalised dct matrix
    pdtmat = dct(eye(s));
else % normalised hadamard matrix (make sure s is 2^k form) 
    pdtmat = (1/sqrt(s))*hadamard(s);
end

temp2 = zeros(m,s,n);
for i=1:n
    indx = abs(A(:,i)) > 0;
    temp = zeros(m,s);
    temp(indx,:) = pdtmat;
    temp2(:,:,i) = temp;
end

B = reshape(temp2,[m,s*n]);

end