function Anorm = colNormalise(A)
Anorm = A./sqrt(sum(abs(A).^2));

end