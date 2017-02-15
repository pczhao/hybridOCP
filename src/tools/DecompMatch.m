function coeff = DecompMatch( p, mon, v )
% 
% coeff = DecompMatch(p, mon)
% coeff = DecompMatch(p, mon, v)
% 
% where:
% p -- n-by-1 msspoly
% mon -- m-by-1 msspoly, may not contain all possible monomials of p.
% v -- l-by-1 free msspoly, optional
% 
% Returns:
% When nargin = 2,
%   coeff -- n-by-m array of doubles
% When nargin = 3,
%   coeff -- n-by-m msspoly that can depend on v
% 
% Satisfying:
% p = coeff * mon
% 
% 

if isempty(p)
    coeff = zeros( 1, length(mon) );
    return;
end

if nargin == 2
    [var,pow,M] = decomp(p);
else
    [var,pow,M] = decomp(p,v);
end
basis = recomp(var,pow,speye(size(pow,1)));
mtch = zeros(size(pow,1), 1);
coeff = msspoly(zeros( size(M, 1), length(mon) ));
if (numel(basis) == 0)
    return;
end

for j = 1 : size( mon, 1 )
    t = basis - mon(j);
    [~,~,b3] = decomp(t);
    mtch( sum(abs(b3),2)==0 ) = j;
end

if prod(mtch) == 0
    warning('Some term(s) of the polynomial is lost in the matching process.');
end

for i = 1 : length(mtch)
    if mtch(i) ~= 0
        coeff(:, mtch(i)) = M(:, i);
    end
end

if nargin == 2
    coeff = double( coeff );
end