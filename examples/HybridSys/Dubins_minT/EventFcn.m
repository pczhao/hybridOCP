function [value,isterminal,direction] = EventFcn(~,x)

value(1) = x(1);
value(2) = x(2);

direction = zeros(1,2);
isterminal = zeros(1,2);
