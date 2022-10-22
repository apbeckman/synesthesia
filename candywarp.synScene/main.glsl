bool invert = (invert_bool > 0.5); 


////////////////////////////////////////////////////////////
// CandyWarp  by mojovideotech
//
// based on :  
// glslsandbox.com/e#38710.0
// Posted by Trisomie21
// modified by @hintz
//
// Creative Commons Attribution-NonCommercial-ShareAlike 3.0
////////////////////////////////////////////////////////////




vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);

	float s = RENDERSIZE.y / scale;
	float radius = RENDERSIZE.x / cycle;
	float gap = s * (1.0 - thickness);
	vec2 pos = _xy.xy - RENDERSIZE.xy * 0.5;
	float d = length(pos);
	float T = ((smoothTime*0.5) + rate);
	d += warp * (sin(pos.y * 0.25 / s + T) * sin(pos.x * 0.25 / s + T * 0.5)) * s * 5.0;
	float v = mod(d + radius / (loops * 2.0), radius / loops);
	v = abs(v - radius / (loops * 4.0));
	v = clamp(v - gap, 0.0, 1.0);
	d /= radius - T;
	vec3 m = fract((d - 1.0) * vec3(loops * hue, -loops, loops * tint) * 0.5);
	if (invert) 	out_FragColor = vec4(m / v, 1.0);
	else out_FragColor = vec4(m * v, 1.0);

return out_FragColor; 
 } 
