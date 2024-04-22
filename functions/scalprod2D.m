function p = scalprod2D(f1x,f1y,f2x,f2y,xk,yk,ajk)
% p = scalprod(f1x,f1y,f2x,f2y,xk,yk,ajk) computes the scalar product 
% between the vectorial functions f1 and f2 using the GLL nodes (xk,yk) and 
% the associated weights ak. xk and yk must me matrices of appropriate
% dimensions.
% If f1 and f2 contain the values computed over (xk,yk), leave the fields 
% xk and yk empty.

N = size(xk,1);
if isa(f1x, 'function_handle'), fk1x = f1x(xk,yk); else, fk1x = f1x; end
if isa(f1y, 'function_handle'), fk1y = f1x(xk,yk); else, fk1y = f1y; end
if isa(f2x, 'function_handle'), fk2x = f2x(xk,yk); else, fk2x = f2x; end
if isa(f2y, 'function_handle'), fk2y = f2x(xk,yk); else, fk2y = f2y; end

p = 0;
for m = 1 : N+1
  for n = 1 : N+1
    p = p + ajk(m,n).*(fk1x(m,n).*fk2x(m,n) + fk1y(m,n).*fk2yx(m,n));
  end
end


