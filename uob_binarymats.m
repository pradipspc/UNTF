function Phi = uob_binarymats(q,n,r,k)

% here q is a prime, n is positive integer, k is less or equal to q^n, r is less than k. %this function generates binary matrix of size kq^n x (q^n)^(r+1) with k ones and (q^n)^r ones in each column and row, respectively.


p = q^n;
% elements of the field Fp
field = gftuple((-1:p-2)',n,q);

% constructing the matrix by using different polynomials for different columns
Phi  = zeros(p*k, p^(r+1));
for Poly_ind = 1 : p^(r+1)
    % the coefficients of the polynomials, starting from the highest
    % degree, i.e. [f0 f1 ..] represents f0x^r + f1x^(r-1) etc.
    PolyCoeff   = de2bi(Poly_ind - 1 ,  r + 1  ,  p, 'left-msb') -1;
    % the -ones represent the zero value in the field
    PolyVal     = -ones(1 , k);
    % just add the various powers to evaluate the polynomial
    for coef_ind = 1 : r+1
        PolyVal = gfmul(PolyVal , (-1 : k - 2) , field);
        PolyVal = gfadd(PolyVal , repmat(PolyCoeff(coef_ind) , 1 , k) , field);
    end
    % build the matrix
    PolyVal(PolyVal == -inf)    = -1;
    Phi(2 + PolyVal + p * (0 : k - 1)  ,  Poly_ind) = 1;
end