%% This function calculates dvip/dvp or dvep/dvp for each intrinsic cortical column
function dv_xp_dv_p = diagonal_from_x(alpha_xp, alpha_px, tau_xp, tau_px, v0, varsigma, vp)
    numerator = 0.1592*alpha_xp*alpha_px*tau_xp*tau_px*exp((0.5*(v0+alpha_px*tau_px*(0.5*erf(0.7071*(v0-vp)/varsigma)-0.5))^2)/(-varsigma^2))*exp((0.5*(v0-vp)^2)/(-varsigma^2));
    dv_xp_dv_p = numerator/varsigma^2;
end