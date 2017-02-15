function [value,isterminal,direction] = EventFcn(~,x)

value = zeros(2,1);
direction = zeros(2,1);
isterminal = zeros(2,1);

value(1) = x(2) - 0.3;
direction(1) = 1;
isterminal(1) = 1;

value(2) = x(1)^2 + x(2)^2 - 0.3;
direction(2) = 0;
isterminal(2) = 0;