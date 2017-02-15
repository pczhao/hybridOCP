function [out] = Stance_f_Approx( x, params )
l0 = params.l0;
g = params.g;
k = params.k;
m = params.m;
umax = params.umax;

if size(x,1) == 1
    x = x';
end
q1 = x(1,:);
q2 = x(2,:);
q3 = x(3,:);
q4 = x(4,:);

out = ...
    [ q2;
      (1/2).*m.^(-1).*(m.*(g.*((-2)+q3.^2)+2.*q1.*q4.^2)+ ...
            k.*(2.*l0+(-2).*q1));
      q4;
      (1/6).*l0.^(-3).*(g.*q3.*((-18).*l0.*q1+6.* ...
            q1.^2+(-1).*l0.^2.*((-18)+q3.^2))+12.*l0.*((-2).*l0+q1).*q2.*q4);
      (-1).*q2.*q3+(-1).*q1.*q4+(1/2).*l0.*q3.^2.*q4];