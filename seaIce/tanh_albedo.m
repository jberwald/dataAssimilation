function a = tanh_albedo( ai, ao, E )
a = ao + (ai - ao) * -min( tanh(E), 0 );
