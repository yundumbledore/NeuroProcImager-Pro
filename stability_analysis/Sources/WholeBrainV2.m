function [eq,H_diagonal] = WholeBrainV2(options, alpha_ip, alpha_pi, alpha_pe, alpha_ep, varsigma, v0, tau_in, tau_ex, x0, n_channels, w, mu_avg, H_diagonal, scale)
    fun = @coupled_nmm;
    eq = fsolve(fun,x0,options);
    
    function F = coupled_nmm(x)
%         vp = H_diagonal*x + mu_avg;
%         gvp = g(vp/scale,v0,varsigma);
        
        for i = 1:n_channels
%             mu = w(i,:)*gvp;
            v_ip = x(1+(i-1)*9);
            z_ip = x(2+(i-1)*9);
            v_pi = x(3+(i-1)*9);
            z_pi = x(4+(i-1)*9);
            v_pe = x(5+(i-1)*9);
            z_pe = x(6+(i-1)*9);
            v_ep = x(7+(i-1)*9);
            z_ep = x(8+(i-1)*9);
            v_p = x(9+(i-1)*9);
            w_i = w(i,:);


            F(9*(i-1)+1) = scale*z_ip; % v_ip
            F(9*(i-1)+2) = alpha_ip(i)*g(v_pi/scale, v0, varsigma) - 2*z_ip/tau_in - v_ip/(scale*tau_in^2); % z_ip

            F(9*(i-1)+3) = scale*z_pi; % v_pi
            F(8*(i-1)+4) = alpha_pi(i)*g(v_p/scale, v0, varsigma) - 2*z_pi/tau_ex - v_pi/(scale*tau_ex^2); % z_pi

            F(9*(i-1)+5) = scale*z_pe; % v_pe
            F(9*(i-1)+6) = alpha_pe(i)*g(v_p/scale, v0, varsigma) - 2*z_pe/tau_ex - v_pe/(scale*tau_ex^2); % z_pe

            F(9*(i-1)+7) = scale*z_ep; % v_ep
            F(9*(i-1)+8) = alpha_ep(i)*g(v_pe/scale, v0, varsigma) - 2*z_ep/tau_ex - v_ep/(scale*tau_ex^2); % z_ep
            
            F(9*(i-1)+9) = sum(w.*g(x(9:9:end)));
        end
    end

    function out = g(v,v0,varsigma)
        out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
    end
end