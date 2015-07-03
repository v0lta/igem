function [ mu ] = muFunction( h,DRhoMin,DRhoMax )
% The mu-cell mobility function.

Kh = 4*10^8;
%m = 20; <- leads to trouble....
m = 5;

mu = (DRhoMax + DRhoMin .* (h./Kh).^m)./(1 + (h./Kh).^m);
end

