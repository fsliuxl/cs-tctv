function X_rec = tctv_cs(A, y, lambda, opts)
%% Problem formulation
%  min_x  1/2 || A x - y||^2 + lambda\sum_{i=1}^3  ||\nabla_i X||_{*L}
%  A: measurement operator
%  y: observed signal
%  lambda: penalty parameters for TCTV norm
%  X: true image
% version 1.0 - 30/11/2022
%
% Written by Xinling Liu (fsliuxl@163.com)
%% Main
% Initialization
dim = opts.dim;
X0 = zeros(dim);
mu = opts.mu;
rho = opts.rho;
shift_dim = opts.shift_dim;
transform = opts.transform;
MaxIter = opts.MaxIter;
tol = opts.tol;
X = opts.X; % used for calculate quality indices
Z1 = X0; Z2 = X0; Z3 = X0;
DualZ1 = Z1; DualZ2 = Z2; DualZ3 = Z3;
% Loop
for k = 1:MaxIter
    tic;
    % Update X
    df1 = diff_1T(Z1 + DualZ1/mu);
    df2 = diff_2T(Z2 + DualZ2/mu);
    df3 = diff_3T(Z3 + DualZ3/mu);
    df = df1 + df2 + df3;
    rhs = (A'*y) + mu * df(:);
    x = tctv_pcg(X0(:),A,rhs(:),mu,dim);
    X1 = reshape(x,dim);
    % Update auxiliary variables Z_i
    Z11 = permute(diff_1(X1) - 1/mu*DualZ1,shift_dim);
    Z22 = permute(diff_2(X1) - 1/mu*DualZ2,shift_dim);
    Z33 = permute(diff_3(X1) - 1/mu*DualZ3,shift_dim);
    Z1 = prox_tnn(Z11,lambda/mu,transform);
    Z2 = prox_tnn(Z22,lambda/mu,transform);
    Z3 = prox_tnn(Z33,lambda/mu,transform);
    Z1 = permute(Z1,shift_dim);
    Z2 = permute(Z2,shift_dim);
    Z3 = permute(Z3,shift_dim);
    % Update dual variables DualZ_i
    DualZ1 = DualZ1 + mu * (Z1 - diff_1(X1)); 
    DualZ2 = DualZ2 + mu * (Z2 - diff_2(X1)); 
    DualZ3 = DualZ3 + mu * (Z3 - diff_3(X1)); 
    relerr = min(norm(X0(:)-X1(:),'fro')/max(norm(X0(:),'fro'),1),1); % relative error
    [psnr1,ssim1] = msqia(X,X1);
    fprintf('Iter = %g, mu1=%.3f, RelErr = %.6f, PSNR = %.4f, SSIM = %.4f\n',k, mu,relerr,psnr1,ssim1);
    if relerr <= tol
        fprintf('<<<<<<<<<<<<<<<<<< ADMM convergenced >>>>>>>>>>>>>>>>\n');
        break;
    end
    X0 = X1;
    mu = min(mu * rho, 1e6); 
end
X_rec = X1;
end