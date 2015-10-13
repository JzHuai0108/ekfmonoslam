function [diffx, diff_hor2]=horizondist(llhref, llhpos)
Ce2n0=llh2dcm_v000(llhref(1:2),[0;1]);
xyz_pos=ecef2geo_v000(llhpos,1);
xyz_ref=ecef2geo_v000(llhref,1);
diffx=Ce2n0*(xyz_pos-xyz_ref);
diff_hor2=sqrt(sum(diffx(1:2).^2));
end