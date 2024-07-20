function eq = SingleNode(options, x0, alpha_ip, alpha_pi, alpha_pe, alpha_ep, mu, varsigma, v0, tau_in, tau_ex, scale)
    fun = @nmm;
    eq = fsolve(fun,x0,options);
    
    function F = nmm(x)
    % x(1) v_ip (2) z_ip (3) v_pi (4) z_pi (5) v_pe (6) z_pe (7) v_ep (8) z_ep
    
        F(1) = scale*x(2); % v_ip
        v_i = x(3);
        F(2) = alpha_ip*g(v_i/scale, v0, varsigma) - 2*x(2)/tau_in - x(1)/(scale*tau_in^2); % z_ip
        
        F(3) = scale*x(4); % v_pi
        v_p = x(1) + x(7) + mu;
        F(4) = alpha_pi*g(v_p/scale, v0, varsigma) - 2*x(4)/tau_ex - x(3)/(scale*tau_ex^2); % z_pi
        
        F(5) = scale*x(6); % v_pe
        v_p = x(1) + x(7) + mu;
        F(6) = alpha_pe*g(v_p/scale, v0, varsigma) - 2*x(6)/tau_ex - x(5)/(scale*tau_ex^2); % z_pe
        
        F(7) = scale*x(8); % v_ep
        v_e = x(5);
        F(8) = alpha_ep*g(v_e/scale, v0, varsigma) - 2*x(8)/tau_ex - x(7)/(scale*tau_ex^2); % z_ep
    end

    function out = g(v,v0,varsigma)
        out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
    end
end

% function eq = SingleNode(options, x0, alpha_ip, alpha_pi, alpha_pe, alpha_ep, mu, varsigma, v0, tau_in, tau_ex)
%     fun = @nmm;
%     eq = fsolve(fun,x0,options);
%     
%     function F = nmm(x)
%     % x(1) v_ip (2) z_ip (3) v_pi (4) z_pi (5) v_pe (6) z_pe (7) v_ep (8) z_ep
%         F(1) = x(2); % v_ip
%         F(2) = alpha_ip*0.5*(erf((x(3)-v0)/varsigma)+1)-2*x(2)/tau_in-x(1)/(tau_in^2); % z_ip
%         
%         F(3) = x(4); % v_pi
%         F(4) = alpha_pi*0.5*(erf(((x(1)+x(7))+mu-v0)/varsigma)+1)-2*x(4)/tau_ex-x(3)/(tau_ex^2); % z_pi
%         
%         F(5) = x(6); % v_pe
%         F(6) = alpha_pe*0.5*(erf(((x(1)+x(7))+mu-v0)/varsigma)+1)-2*x(6)/tau_ex-x(5)/(tau_ex^2); % z_pe
%         
%         F(7) = x(8); % v_ep
%         F(8) = alpha_ep*0.5*(erf((x(5)-v0)/varsigma)+1)-2*x(8)/tau_ex-x(7)/(tau_ex^2); % z_ep
%     end
% end