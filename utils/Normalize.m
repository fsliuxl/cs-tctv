function norm_tensor = Normalize(orginal_tensor)
[m,n,p] = size(orginal_tensor);
norm_tensor = zeros([m,n,p]);
for band =1:p
    tmp = orginal_tensor(:,:,band);
    if max(tmp(:)) - min(tmp(:))~=0
        tmp = (tmp-min(tmp(:)))/(max(tmp(:)) - min(tmp(:)));
    else
        tmp = max(tmp(:));
    end
    norm_tensor(:,:,band) = tmp;
end