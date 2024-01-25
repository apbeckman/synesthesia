vec4 renderMain(void)
{
    vec4 fragColor = vec4(0.0);
	vec2 position = _uv;
	position -= _mouse.xy / RENDERSIZE;
    position.x *= 1.667;
	float circle = length(position);

	for (int i = 0; i < 4; i++){
	fragColor += smoothstep(.4, .7, circle)*_pulse(circle, fract(TIME*0.1+i), .01);
	fragColor += smoothstep(.3, .5, circle)*_pulse(1.0-circle, fract(TIME*0.1+i), .01);

	}
	return fragColor;
   
}
