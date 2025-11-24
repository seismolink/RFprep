function seis_out = depthmapping(t_in,seis_in,z_out,ldepth,velp,vels,slow)
% Map time to depth

t_out = z2t(z_out,ldepth,velp,vels,slow);
seis_out = interp1(t_in,seis_in,t_out);

end