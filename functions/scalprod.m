function p = scalprod(f1,f2,xk,ak)
% p = scalprod(f1,f2,xk,ak) computes the scalar product between f1 and f2
% using the GLL nodes xk and the associate weights ak.
% (If f1 and f2 contain the values computed over xk, leave the field xk
% empty; if f1 and f2 are function handles, xk must be specified)

if isa(f1, 'function_handle'), fk1 = f1(xk); else, fk1 = f1; end
if isa(f2, 'function_handle'), fk2 = f2(xk); else, fk2 = f2; end

if isrow(fk1), fk1 = fk1'; end
if isrow(fk2), fk2 = fk2'; end
if iscolumn(ak), ak = ak'; end

p = ak * (fk1.*fk2);


