

			//******** Common Code Begins ********

const int R = 360;

//const vec2 C = vec2(0.5, 0.0);


			//******** BuffA Code Begins ********


vec4 process(vec2 c) {
    vec2 sz = RENDERSIZE.xy;
    vec4 v0 = texture(BuffA, (c + vec2(-1,0)) / sz);
    vec4 v1 = texture(BuffA, (c + vec2( 1,0)) / sz);
    vec4 v2 = texture(BuffA, (c + vec2(0,-1)) / sz);
    vec4 v3 = texture(BuffA, (c + vec2(0, 1)) / sz);    
    vec4 v4 = texture(BuffA, c / sz);
    float w = ((int(FRAMECOUNT) % 2) == 0)?0.367879:3.0;
    float k = (1.0 - w) / 4.0;    
    return w * v4 + k * (v0 + v1 + v2 + v3); 
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 C = _mouse.xy / RENDERSIZE.xy;
    if (int(FRAMECOUNT) <= 1) {
        vec2 uv = (fragCoord / RENDERSIZE.xy) * 2.0 - 1.0;
                
        uv.x *= RENDERSIZE.x / RENDERSIZE.y;
        float c0 = length(uv) - 0.5;
        float c1 = length(uv - vec2(0.0,0.25));
        float c2 = length(uv - vec2(0.0,-0.25));
        float d = max(min(max(c0, -uv.x),c1-0.25),-c2+0.25);        
        d = max(d, -c1+0.07);
        d = min(d, c2-0.07);
        d = min(d, abs(c0)-0.01);
        float w = sign(d) * sin(c0*80.0);
        fragColor = vec4(w);
    } else if ((int(FRAMECOUNT) % R) == 0) {
        vec2 c = vec2(fragCoord / RENDERSIZE.xy);
        fragColor = texture(BuffA, (c - C) / 2.0 + C);
    } else {
        vec4 s = process(fragCoord);
        fragColor = clamp(s, -1.0, 1.0);
    }    
	return fragColor; 
 } 


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 C = _mouse.xy / RENDERSIZE.xy;
    vec2 uv = fragCoord / RENDERSIZE.xy;
    float m = float(int(FRAMECOUNT) % R) / float(R);  
    uv -= C;
    uv.x *= RENDERSIZE.x / RENDERSIZE.y;
    uv = uv / (1.0 + 0.1*dot(uv,uv));
    uv.x /= RENDERSIZE.x / RENDERSIZE.y;
    uv *= exp(mix(log(1.0),log(0.5),m));
    float f = 0.004;
    vec2 uv3 = uv * exp(0.0) + C;
    vec2 uv2 = uv * exp(-f) + C;
    vec2 uv1 = uv * exp(-f*2.0) + C;
    float r = texture(BuffA, uv1).r*0.5+0.5;
    float g = texture(BuffA, uv2).r*0.5+0.5;
    float b = texture(BuffA, uv3).r*0.5+0.5;    
    fragColor = vec4(pow(vec3(r,g,b), vec3(0.5)),1.0);
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