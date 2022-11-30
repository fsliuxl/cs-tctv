function psnr = PSNR_me(Xfull,Xrecover)

maxP = max(abs(Xfull(:)));
Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);
MSE = norm(Xfull(:)-Xrecover(:))^2/(numel(Xrecover));
psnr = 10*log10(maxP^2/MSE);