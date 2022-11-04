vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

// Fork of "Day 94" by jeyko. https://shadertoy.com/view/WsXcDM
// 2020-03-22 08:59:44

#define pmod(p,x) mod(p,x) - 0.5*x
#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))


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
    
    
    n += valueNoise(uv*500.,0.1);
    //n += valueNoise(uv*1000., 0.1)*0.5;
    n += valueNoise(uv*1700.,0.1)*0.5;
    
    n -= valueNoise(uv*10.,1.)*1.;
    n -= valueNoise(uv*20.,0.5)*0.5;
    
    
    n = max(n, 0.);
    return n;
}


vec3 get(vec2 fragCoord )
{
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
	vec2 uvn = uv;
    uv *= 1. - dot(uv,uv)*0.3;
    
    uv *= 2. + valueNoise(vec2(TIME), 0.2).x*1.5;
   
    vec3 col = vec3(0);
    
    
    
    for(int i = 0; i<3; i++){
        if(uv.x < uv.y){uv.xy = uv.yx;}
        uv *= rot(0.5);
        uv.x -= 0.9;
        
    }
    
    float m = 1.;
    for(float i = 0.; i< 1.; i++){
        //uv = uv.yx;
    	float dpp = dot(uv,uv);
        dpp = clamp(dpp, 0.1,1.);
		
        vec2 z1 = vec2(0.4);
        vec2 z2 = vec2(-0.1);
		z1 = uv - z1; uv -= z2;        
    	uv = mat2(z1, z1.y, -z1.x)*(uv*1. + 1.*vec2(-0.1, 0.7))/ dpp; 
        
        uv -= vec2(0.1,0.);
    }
    
    uv = sin(uv*3.14*(0.6));
    
    float dC = log(length(uv) + 90.)+ TIME*0.001;
    
    
    float modD = 0.00099;
    float id = floor(dC/modD);
    dC = pmod(dC, modD);
    
    
    vec3 c = pal(0.5,0.5,vec3(0.1,0.6,0.9), 0.6 + sin(id*0.4)*0.4, id*0.9 + TIME*0.5);
    
    c = max(c, 0.);
    c = smoothstep(0.,1.,c);
    col += c;
    
    
    
    //col -= fbm(uv).x*0.001;
    //col += fbm(uv + 4.).x*0.001;
    
    
    col = pow(col, vec3(1. + dot(uvn,uvn)));
    
    col *= 1.  - dot(uvn,uvn);
        
    //col += smoothstep(0.01);

    return col;
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 col = vec3(0);
    
    // float(min(FRAMECOUNT,0)) hack for faster compilation (stops loop unrolling)
    for(float i =0.; i < 9. + float(min(FRAMECOUNT,0)); i++){
    	col += get(fragCoord + vec2(mod(i,3.),floor(i/3.))*0.25);
    }
    col /= 9.;
    
	col = pow(col, vec3(0.4545));
    
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
    float steps = 16.;
    float scale = 0.00 + pow(dot(uvn,uvn),3.1)*3.93;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*2.3;
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
    
    fragColor.b *= 0.97 + dot(uvn,uvn)*1.9;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.4);
    
    fragColor.r *= 1.  - smoothstep(0.,1.,dot(uvn,uvn))*0.9;
    
    
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