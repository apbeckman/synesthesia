

			//******** BuffA Code Begins ********
float dist(vec2 p0, vec2 pf){return sqrt((pf.x-p0.x)*(pf.x-p0.x)+(pf.y-p0.y)*(pf.y-p0.y));}
float d2 = dist(RENDERSIZE.xy*0.5,_xy.xy)*0.005;
#define loops 15

vec2 offset [loops];
float kernel [loops];
//float highhits = pow(syn_HighLevel*0.25 + (syn_Level*0.25+syn_HighHits*0.35), 2.0)*syn_Intensity;
//float bass = pow(syn_BassLevel*0.25 + (syn_Level*0.25+syn_BassHits*0.5), 2.0)*syn_Intensity;
    float Noisey = Noise*low;
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
    //fragCoord += (syn_BassLevel*0.5+syn_Presence*0.5) *vec2(_noise(TIME)*syn_BassLevel, _noise(TIME)*syn_BassLevel);

    vec2 res = RENDERSIZE.xy;
    float t = smoothTime*0.25;
    vec2 tc = fragCoord / res;
    vec2 uv = -1.0 + 2.0 * tc;
    float punch_in_out = punch_out - punch_in;
    //zoom 
    uv *= 0.995+(Zoom*0.5-sign(Zoom)*(basshits*0.0125*sign(Zoom)));
    uv += uv*sin(mix(_noise(_rotate(_uvc, TIME*0.1)), _noise(TIME), 0.25))*warp*0.01;
    uv = uv * 0.5 + 0.5;
    uv += (d2)*_uvc*(Fisheye+punch_in_out*0.1)* (1.0 + 2.0 * low + syn_Intensity * 1.5)*0.025;
    //rotation
    float twist = spin_right - spin_left;
    uv = rotate(uv, 0.00+(rotat+twist*0.025)*d2)+_uvc*twist*0.01;
    
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
    float outer = -0.25 * Convolution;
    float inner = 0.125 * Convolution;
    
    kernel[0] = 0.0; kernel[1] = 0.0; kernel[2] = inner+test;
    kernel[3] = 0.0; kernel[4] = inner+test; kernel[5] = 0.0;
    kernel[6] = outer-test; kernel[7] = 0.0; kernel[8] = 0.0;
    
    vec4 sum = texture(BuffB, uv);
    vec4 last = sum;
    last.gb = -1.0 + 2.0 * last.gb;

    for (int i = 0; i < loops; i++) {
        vec4 color = texture(BuffA, (uv + offset[i]+.05*_noise(offset[i])) - last.gb*0.25);
        sum.rgb -= vec3(_fbm(_uv))*Noisey*0.00125;
        sum += color * kernel[i];
       // sum.rgb -= _noise(sin(-1*color.rgb))*0.0001*Noisey;
    }
    
    vec4 src = texture(media_pass, tc);
    // src *= 0.75 + _edgeDetectSobel(syn_UserImage)*0.5;
    // src = mix(src, texture(feedback, _uv), .15*(syn_HighLevel)) ;
    vec4 mediaEdges = _edgeDetectSobelMedia();
    sum.rgb = mix(sum.rgb, src.rgb-(sum.rgb)*0.001*syn_BassLevel, src.rgb*(0.65)*Media);
    sum.rgb -= .0125*(_noise(texture(feedback, _uv).rgb*0.01));
    if(FRAMECOUNT < 10 || _mouse.z > 0.0){
        fragColor = src;
    } else {
    	fragColor = vec4(clamp(vec3(sum.rgb), vec3(0.0), vec3(1.0)), 1.0);
    }
	return fragColor; 
 } 


vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
	fragColor = texture(BuffA, uv);
    fragColor.rgb -= _fbm(fragColor.rgb)*0.001*Noisey*syn_BassLevel*d2;
    fragColor.rgb += _edgeDetectSobel(BuffA, _uv).rgb*.0000001;
	return fragColor; 
 } 

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord.xy / RENDERSIZE.xy;
	fragColor = texture(BuffA, uv);
	return fragColor; 
 } 
vec4 mediaPass() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    vec4 media = _loadMedia();

	return media;
 } 
vec4 Feedback() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    vec4 fB = texture(syn_FinalPass, _uv);
	return fB;
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return mediaPass();
	}
	if(PASSINDEX == 3){
		return Feedback();
	}
	if(PASSINDEX == 4   ){
		return renderMainImage();
	}
}