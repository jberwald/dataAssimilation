function a = original_albedo( ai, ao, E, tanha )

if tanha>0 % using gradual albedo transition
    a=ao+(ai-ao)/2*(1-tanh(E/tanha));
else
    a=ai*(E<0)+ao*(E>=0);
end
