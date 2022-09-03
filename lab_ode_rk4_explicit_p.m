function X = lab_ode_rk4_explicit_p(fnc, T, X0)
    N=length(T);
    h=T(2)-T(1);
    X = X0;
    for i =2:N
        X_k = X;
        t_k = T(i-1);
        
        K1=fnc(X_k, t_k);
        K2=fnc(X_k+h/2*K1, t_k+h/2);
        K3=fnc(X_k+h/2*K2, t_k+h/2);
        K4=fnc(X_k+h*K3, t_k+h);
        
        X=X_k + h/6*(K1+2*K2+2*K3+K4);
    end
end

