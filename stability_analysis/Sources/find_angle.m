function [CosTheta]= find_angle(u, v)
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
end