function [value, isterminal, direction] = EvtFunc_S2F_scaled( ~, x, params )
% Event function of Stance 1, scaled version
% keyboard;
% isterminal = ones(10,1);
isterminal = [1; zeros(9,1)];
direction = -ones(10,1);
l0 = params.l0;
domain = params.domain{1};

value = [ l0 - x(1);     % transition to Flight 1
          x(1) - domain(1,1);
          domain(2:end,2) - x(2:end);
          x(2:end) - domain(2:end,1)];
