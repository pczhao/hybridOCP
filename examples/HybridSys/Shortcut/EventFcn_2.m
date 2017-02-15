function [value,isterminal,direction] = EventFcn_2(~,x)

value(1) = x(1) - 0.8;
direction(1) = 0;
isterminal(1) = 1;