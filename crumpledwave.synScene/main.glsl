

			//******** Common Code Begins ********


#define R RENDERSIZE.xy

// Dave Hoskins https://www.shadertoy.com/view/4djSRW
vec2 hash23(vec3 p3)
{
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yzx+33.33);
    return fract((p3.xx+p3.yz)*p3.zy);
}

float gyroid (vec3 seed)
{
    return dot(sin(seed),cos(seed.yzx));
}

float fbm (vec3 seed)
{
    float result = 0.;
    float a = .5;
    for (int i = 0; i < 3; ++i) {
        seed += result / 2.;
        result += gyroid(seed/a)*a;
        a /= 2.;
    }
    return result;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv =  (2.0 * fragCoord - RENDERSIZE.xy) / min(RENDERSIZE.x, RENDERSIZE.y);
   
    for(float i = 1.0; i < 8.0; i++){
    uv.y += i * 0.1 / i * 
      sin(uv.x * i * i + smoothTimeC * 0.5) * sin(uv.y * i * i + smoothTimeC * 0.5);
  }
    
   vec3 col;
   col.r  = uv.y - 0.1;
   col.g = uv.y + 0.3;
   col.b = uv.y + 0.95;
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}