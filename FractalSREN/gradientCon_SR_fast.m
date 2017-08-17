function out=gradientCon_SR_fast(in, ker, dx, dy, lambda)


otfk = psf2otf(ker,size(in));
Denom1 = conj(otfk).*otfk;
Nomin1 = conj(otfk).*fft2(in);
clear otfk;

Denom2 = abs(psf2otf([-1 1],size(in))).^2+abs(psf2otf([-1 1]',size(in))).^2;

Fx = conj(psf2otf([-1 1],size(in)));
Fy = conj(psf2otf([-1 1]',size(in)));
Fx_t = fft2(dx);
Fy_t = fft2(dy);
F_yout = (Nomin1+lambda*(Fx.*Fx_t+Fy.*Fy_t))./(Denom1+lambda*Denom2);
clear Nomin1;
clear Fx;
clear Fx_t;
clear Fy;
clear Fy_t;
clear Denom1;
clear Denom2;
out = real(ifft2(F_yout));
clear F_yout;