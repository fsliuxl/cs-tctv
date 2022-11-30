function x= tctv_pcg(x,A,T,mu,dim)
  temp =T(:);
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
        x1 = reshape(x,dim);
        df = diff_1T(diff_1(x1)) + diff_2T(diff_2(x1)) + diff_3T(diff_3(x1));
        y = A'*(A*x) + mu*df(:);
    end
end