%% This function calculates dmu/dvp
function dmu_dv_p = offdiagonalelements(w, v0, vp, varsigma)
    dmu_dv_p = (0.3989*w*exp((0.5*(v0-vp)^2)/(-varsigma^2)))/varsigma;
end