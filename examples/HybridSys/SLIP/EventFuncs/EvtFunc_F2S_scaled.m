function [value, isterminal, direction] = EvtFunc_F2S_scaled( ~, x, params )
% Event function of Flight 2, scaled version
isterminal = ones(8,1);
direction = -ones(8,1);
domain = params.domain{3};
yR = params.yR;

value = [ x(3) - yR;                    % y = yR
          domain(1:2,2) - x(1:2);
          x(1:2) - domain(1:2,1);
          domain(3,2) - x(3);
          domain(4,2) - x(4);
          x(4) - domain(4,1) ];
