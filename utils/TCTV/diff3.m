function diff_x = diff3(x,dim)

tenX   = reshape(x,dim);

dfx1     = diff(tenX, 1, 1);
dfy1     = diff(tenX, 1, 2);
dfz1     = diff(tenX, 1, 3);

dfx      = zeros(dim);
dfy      = zeros(dim);
dfz      = zeros(dim);
dfx(1:end-1,:,:) = dfx1;
dfx(end,:,:)     =  tenX(1,:,:) - tenX(end,:,:);
dfy(:,1:end-1,:) = dfy1;
dfy(:,end,:)     = tenX(:,1,:) - tenX(:,end,:);
dfz(:,:,1:end-1) = dfz1;
dfz(:,:,end) = tenX(:,:,1) - tenX(:,:,end);
diff_x = [dfx(:); dfy(:); dfz(:)];
end

