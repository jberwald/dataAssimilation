function a = linear_albedo( ai, ao, E )

shift = 20; 
scale= 14.06;
if E < -shift
  gamma = -E - shift;
else
  gamma = 0;
end
a = ao + (ai - ao) * min( gamma/14.06, 1 );
