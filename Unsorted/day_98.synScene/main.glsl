//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********


#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))

float xor(float a, float b){
	return float(int(a)^ int(b));
}

float sdBox(vec2 p, vec2 s){
	p = abs(p) - s;
    return max(p.y,p.x);
}
vec4 noise(float t){return texture(image30,vec2(floor(t), floor(t))/256.);}
vec4 valueNoise(vec2 t, float w){
    vec2 fr = fract(t);
	return 
        mix(
            mix( 
                texture(image30,vec2(floor(t.x), floor(t.y))/256.),
                texture(image30,vec2(floor(t.x), floor(t.y) + 1.)/256.),
            	smoothstep(0.,1.,fr.y)
            ),
            mix( 
                texture(image30,vec2(floor(t.x) + 1.,floor(t.y))/256.),
                texture(image30,vec2(floor(t.x) + 1.,floor(t.y) + 1.)/256.),
            	smoothstep(0.,1.,fr.y)
            ),
            smoothstep(0.,1.,pow(fr.x, w)));
}
vec4 fbm(vec2 uv){
	vec4 n = vec4(0);
    n += valueNoise(uv*800.,0.1);
    n += valueNoise(uv*1700.,0.1)*0.5;
    n -= valueNoise(uv*10.,1.)*1.;
    n -= valueNoise(uv*20.,0.5)*0.5;
    n = max(n, 0.);
    return n;
}
vec3 get( vec2 fragCoord )
{
    vec2 p =( fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
	
    
    vec3 col = vec3(0);
    
    
    p *= 1. - dot(p,p)*0.05+test ;
    
    
    vec2 b = p;
    float d = 10e6;
    
    p *= 1.;
    
    p = vec2(9.*atan(p.x,p.y)/6.14, length(p));
    
    p.y = log(p.y)*1. - smoothTime*0.16 + valueNoise(vec2(smoothTime*0.15),41.).x*1.;
    vec2 j = p;
    
    float id = floor(p.y);
    vec4 n = noise(id);
    p.xy *= rot(0. + sin(p.y*0.1)*0.);
    
    //vec2 id =
    b *= 1.5;
    //p = cos(p*(1. - exp(-length(b)*length(b)*1000.)*0.2));
    p = cos(p*(1. - exp(-length(b)*length(b)*1000.)*0.));
    
    for(float i = 0.; i < 9. + float(min(FRAMECOUNT,0))     ; i++){
    		p = abs(p);
        
        	float sc = (5. + n.y*6.);
        	vec2 q = p;
        
        	float x = xor(q.y*sc,q.x*sc);
        	p *= 1. - x*(0.06 - sin(n.y*4.)*0.02);
        	p -= 0.02 + n.y*0.045;
        	p *= rot(0.125*3.14);
        {
        if(mod(i, 2.) == 0.)
    		d = min(d, sdBox(p, vec2(0.04 + x*0.005)));
        else
    		d = max(d, -sdBox(p, vec2(0.02 + x*0.005))); 
        
        }
    }

    
    col += smoothstep(0.001,0.,d);
    
    
    float f = smoothstep(0.,1.,fbm(j*80.).x);
    
    f = pow(f, 6.)*2.4;
    if(col.x < 0.4)
        col += f*6.5;
    
    
    col -= f;
    
    
    
    //col += fbm(p).x*10.;
    //col.r += sin(id)*3.5;
	return col;
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 col = vec3(0);
    
    float aa = 3.;
    
    for(float i =0.; i < aa*aa + float(min(FRAMECOUNT,0)); i++){
    	col += get(fragCoord + vec2(mod(i,aa),floor(i/aa))/aa);
    }
    col /= aa*aa;
    
    
    col = max(col, 0.);
	col = pow(col, vec3(0.4545));
    
    col = 1. - col;
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 




// radiual blur in this buffer

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
	fragColor = vec4(0);
    // Radial blur
    float steps = 7.;
    float scale = 0.00 + pow(dot(uvn,uvn),1.1)*0.04;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*1.2;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
        scale *= 0.97;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1.4/RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    fragColor += radial;
    
    fragColor.b *= 0.97 + dot(uvn,uvn)*0.4;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.6);
    
    fragColor.t *= 1.  - smoothstep(0.,1.,dot(uvn,uvn))*0.;
    
    
    //fragColor.b *= 1. + uv.x*0.02;
    //fragColor.g *= 1. + uv.t*0.05;
    
    fragColor = max(fragColor, 0.);
    fragColor = pow(fragColor, vec4(0.4545 ));
    fragColor *= 1. - dot(uvn,uvn)*1.   ;
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