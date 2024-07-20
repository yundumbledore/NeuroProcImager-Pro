function [eq] = WholeBrain(options, alpha_ip, alpha_pi, alpha_pe, alpha_ep, varsigma, v0, tau_in, tau_ex, x0, n_channels, mu_avg, scale)
    fun = @coupled_nmm;
    eq = fsolve(fun,x0,options);
    
    function F = coupled_nmm(x)
%         vp = H_diagonal*x + mu_avg;
%         gvp = g(vp/scale,v0,varsigma);
        
        for i = 1:n_channels
%             mu = w(i,:)*gvp;
            
            F(8*(i-1)+1) = scale*x(8*(i-1)+2); % v_ip
            v_i = x(8*(i-1)+3);
            F(8*(i-1)+2) = alpha_ip(i)*g(v_i/scale, v0, varsigma) - 2*x(8*(i-1)+2)/tau_in - x(8*(i-1)+1)/(scale*tau_in^2); % z_ip

            F(8*(i-1)+3) = scale*x(8*(i-1)+4); % v_pi
            v_p = x(8*(i-1)+1) + x(8*(i-1)+7) + mu_avg(i);
            F(8*(i-1)+4) = alpha_pi(i)*g(v_p/scale, v0, varsigma) - 2*x(8*(i-1)+4)/tau_ex - x(8*(i-1)+3)/(scale*tau_ex^2); % z_pi

            F(8*(i-1)+5) = scale*x(8*(i-1)+6); % v_pe
            v_p = x(8*(i-1)+1) + x(8*(i-1)+7) + mu_avg(i);
            F(8*(i-1)+6) = alpha_pe(i)*g(v_p/scale, v0, varsigma) - 2*x(8*(i-1)+6)/tau_ex - x(8*(i-1)+5)/(scale*tau_ex^2); % z_pe

            F(8*(i-1)+7) = scale*x(8*(i-1)+8); % v_ep
            v_e = x(8*(i-1)+5);
            F(8*(i-1)+8) = alpha_ep(i)*g(v_e/scale, v0, varsigma) - 2*x(8*(i-1)+8)/tau_ex - x(8*(i-1)+7)/(scale*tau_ex^2); % z_ep
        end
    end

    function out = g(v,v0,varsigma)
        out = 0.5*erf((v - v0) / (sqrt(2)*varsigma)) + 0.5;
    end
end

% function [eq,H_diagonal] = WholeBrain(options, alpha_ip, alpha_pi, alpha_pe, alpha_ep, varsigma, v0, tau_in, tau_ex, x0, n_channels, w, mu_avg, H_diagonal)
%     fun = @coupled_nmm;
%     eq = fsolve(fun,x0,options);
%     
%     function F = coupled_nmm(x)
%         vp = H_diagonal*x' + mu_avg;
%         gvp = g_function(vp,v0,varsigma);
%         
%         for i = 1:n_channels
%             mu = w(i,:)*gvp;
%             
%             F(8*(i-1)+1) = x(8*(i-1)+2); % v_ip
%             F(8*(i-1)+2) = alpha_ip(i)*0.5*(erf((x(8*(i-1)+3)-v0)/varsigma)+1)-2*x(8*(i-1)+2)/tau_in-x(8*(i-1)+1)/(tau_in^2); % z_ip
% 
%             F(8*(i-1)+3) = x(8*(i-1)+4); % v_pi
%             F(8*(i-1)+4) = alpha_pi(i)*0.5*(erf(((x(8*(i-1)+1)+x(8*(i-1)+7))+mu-v0)/varsigma)+1)-2*x(8*(i-1)+4)/tau_ex-x(8*(i-1)+3)/(tau_ex^2); % z_pi
% 
%             F(8*(i-1)+5) = x(8*(i-1)+6); % v_pe
%             F(8*(i-1)+6) = alpha_pe(i)*0.5*(erf(((x(8*(i-1)+1)+x(8*(i-1)+7))+mu-v0)/varsigma)+1)-2*x(8*(i-1)+6)/tau_ex-x(8*(i-1)+5)/(tau_ex^2); % z_pe
% 
%             F(8*(i-1)+7) = x(8*(i-1)+8); % v_ep
%             F(8*(i-1)+8) = alpha_ep(i)*0.5*(erf((x(8*(i-1)+5)-v0)/varsigma)+1)-2*x(8*(i-1)+8)/tau_ex-x(8*(i-1)+7)/(tau_ex^2); % z_ep
%         end
%     end
% end