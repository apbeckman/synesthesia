

			//******** BuffA Code Begins ********

vec2 offset [9];
float kernel [9];
//float highhits = pow(syn_HighLevel*0.25 + (syn_Level*0.25+syn_HighHits*0.35), 2.0)*syn_Intensity;
//float bass = pow(syn_BassLevel*0.25 + (syn_Level*0.25+syn_BassHits*0.5), 2.0)*syn_Intensity;
vec2 rotate(vec2 coords, float angle){
	float sin_factor = sin(angle );
    float cos_factor = cos(angle );
    coords = vec2((coords.x - 0.5) , coords.y - 0.5) * mat2(cos_factor, sin_factor, -sin_factor, cos_factor);
    coords += 0.5;
    return coords;
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 res = RENDERSIZE.xy;
    float t = smoothTime*0.25;
    vec2 tc = fragCoord / res;
    vec2 uv = -1.0 + 2.0 * tc;
    
    //zoom 
    uv *= 0.995+(Zoom-sign(Zoom)*(basshits*0.0125*sign(Zoom)));
    uv = uv * 0.5 + 0.5;
    
    //rotation
    uv = rotate(uv, 0.00+rotat);
    
    vec2 step = (.50+StepSize-(basshits*0.5)) / res;
    
    offset[0] = vec2(-step.x, -step.y);
    offset[1] = vec2(0.0, -step.y);
    offset[2] = vec2(step.x, -step.y);
    
    offset[3] = vec2(-step.x, 0.0);
    offset[4] = vec2(0.0, 0.0);
    offset[5] = vec2(step.x, 0.0);
    
    offset[6] = vec2(-step.x, step.y);
    offset[7] = vec2(0.0, step.y);
    offset[8] = vec2(step.x, step.y);
    
    //convolution values
    float outer = -0.25;
    float inner = 0.125;
    
    kernel[0] = 0.0; kernel[1] = 0.0; kernel[2] = inner+test;
    kernel[3] = 0.0; kernel[4] = inner+test; kernel[5] = 0.0;
    kernel[6] = outer-test; kernel[7] = 0.0; kernel[8] = 0.0;
    
    vec4 sum = texture(BuffA, uv);
    vec4 last = sum+test;
    last.gb = -1.0 + 2.0 * last.gb;
    
    for (int i = 0; i < 9; i++) {
        vec4 color = texture(BuffA, (uv + offset[i]) - last.gb*0.0125);
        sum += color * kernel[i];
    }
    
    vec4 src = texture(syn_UserImage, tc);
    sum.rgb = mix(sum.rgb, src.rgb, src.rgb*(0.1+hits*1.));
    
    if(FRAMECOUNT < 10 || _mouse.z > 0.0){
        fragColor = src;
    } else {
    	fragColor = vec4(clamp(vec3(sum.rgb), vec3(0.0), vec3(1.0)), 1.0);
    }
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
	fragColor = texture(BuffA, uv);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}