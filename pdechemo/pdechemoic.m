function u0 = pdechemoic(x)

u10 = 15*(x - x.^2).*exp(-50 .*((x - 0.5).^2));
u20 = 1;
u0 = [u10; u20];
