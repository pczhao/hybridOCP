function [value, isterminal, direction] = EvtFunc_F1_scaled( ~, x, params )
% Event function of Flight 1, scaled version
isterminal = ones(8,1);
direction = -ones(8,1);
domain = params.domain{2};

value = [ x(4);                    % y_dot = 0
          domain(1:3,2) - x(1:3);
          x(1:3) - domain(1:3,1);
          domain(4,2) - x(4); ];
