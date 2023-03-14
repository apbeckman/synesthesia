

			//******** BuffA Code Begins ********


vec3 glow = vec3(0);
vec3 glowF = vec3(0);
    
#define TIME (TIME + 10.)
#define mx (bass_time*1.9)
#define my (20.*_mouse.y/RENDERSIZE.x)
#define pal(a,b,c,d,e) (a + b*sin(c*d + e))
    
#define dmin(a, b) a.x < b.x ? a : b
#define pmod(a,x) mod(a,x) - x*0.5

#define pi acos(-1.)


vec2 NOISE = vec2(0.);
vec2 valueNoise(float p){
	vec2 a = texture(image30, vec2(floor(p))/1024.).xy;
	vec2 b = texture(image30, vec2(floor(p) + 1.)/1024.).xy;
    return mix(a,b,smoothstep(0.,1.,fract(p)));
}

vec3 path(float z){
    z *= 0.44;
	return vec3(
    	sin(z + cos(z)*0.016),
    	cos(z + sin(z*0.8)*0.025),
    	0.
    )*0.376;
}
float opSmoothUnion(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1)/k,0.,1.);
    return mix(d2, d1, h) - k*h*(1. - h);  
}
// Hex code from BigWings! He has a tutorial on them.
float HexDist(vec2 p) {
	p = abs(p);
    float c = dot(p, normalize(vec2(1,1.73)));
    c = max(c, p.x);
    return c;
}


vec4 HexCoords(vec2 uv) {
	vec2 r = vec2(1, 1.73);
    vec2 h = r*.5;
    vec2 a = mod(uv, r)-h;
    vec2 b = mod(uv-h, r)-h;
    vec2 gv = dot(a, a) < dot(b,b) ? a : b;
    float x = atan(gv.x, gv.y);
    float y = .5-HexDist(gv);
    vec2 id = uv-gv;
    return vec4(x, y, id.x,id.y);
}



float modu;
float moduB;
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
vec2 map(vec3 p){
	vec2 d = vec2(10e6);
	p -= path(p.z);
    vec2 n = normalize(p.xy);
    #define modDist 1.
    #define tunnW 0.7
    #define pipeW 0.04
    
    vec3 k = p;
    float id = floor(p.z/modDist);
    vec3 g = p;
    p.z = pmod(p.z, modDist);
    
    vec3 o = p;
    //p.xy *= rot(0.4 + p.z*(0.1 + sin(TIME*0.1) )+ TIME*0.3);
    vec2 pC = vec2(atan(p.y,p.x), length(p.xy));
    
    vec3 q = vec3(pC, p.z);

    vec4 hc = HexCoords(vec2(4.*pC.x/pi, p.z*2.)*1.);
    
    float dHex = hc.y - 0.02 + sin(o.z*100.)*0.01;
    dHex = max(dHex, -length(o.xy*1.) + tunnW*0.8);
    //dHex = max(dHex, length(o.xy*1.) - tunnW*0.97);
    //d = dmin(d, vec2(dHex, 2.));    
    
    k = pmod(k, modDist*10.);
    
    float hcy = pmod(hc.y, 0.5);
    float dThings = max(hcy - 0.06,  -hcy - 0.1);
    dThings = max(dThings, -length(o.xy*1.1) + tunnW*0.8);
    
    dThings = max(dThings, -k.z - 0.25*modDist);
    dThings = max(dThings, k.z - 0.25*modDist);
    
    //dThings = max(pmod());
    if(mod(floor(4.*pC.x/ pi ), 2.) != 1.)
    	d = dmin(d, vec2(dThings, 2.)); 
    
    //hc.y += 0.03;
    float dThingsB = max(hc.y - 0.03,  -hc.y - 0.1);
    dThingsB = max(dThingsB, -length(o.xy*1.1) + tunnW*0.8);
        
        
    pC.x += 0.25*pi;
    if(mod(floor(4.*pC.x/ pi ), 2.) != 1.)
    	d = dmin(d, vec2(dThingsB, 2.)); 
    
    
    // dots
    
    
    // mod
    float mm = sin(TIME*0.5 + g.z*0.5 + p.z);
    modu = (mm/sqrt(0.02 + mm*mm ))*0.5 + 0.5;
    float mmB = sin(TIME*1.25 + g.z*0.25 + p.z*0.8 + p.y);
    moduB = (mmB/sqrt(0.01 + mmB*mmB ))*0.5 + 0.5;
    
    d = dmin(d, vec2(dHex, 2.));
    glow += exp(-d.x*100.)*pal(1.3,0.7,vec3(1.8+modu*0.5,0.4,0.8), 3.9 +modu*0.2 + sin(p.z)*0.5,2.)*2.;
    // tunnel
    float dTunn = -length(o.xy*1.) + tunnW;
    dTunn = max(length(o.xy*1.) - tunnW - 0.02, dTunn);
    d = dmin(d, vec2(dTunn, 10.));
    
    
    glowF += exp(-d.x*(10. + NOISE.x*200.9))*pal(1.39,0.7+ exp(-dThingsB*20.)*0.6,vec3(1.1+modu*0.5 ,0.4,0.8), 3.9 +modu*0.1 + sin(p.z)*0.5,2.)*2.;
    
    
    
    
    pC.x += pi/3.5;
    pC.x = pmod(pC.x, pi/1.);
    pC.y -= 0.6;
    pC = abs(pC) - vec2(0.05,0.15);
    //d = dmin(d, vec2( max(pC.x, pC.y), 2.));
    d.x *= 0.5;
	return d;
}

vec2 march(vec3 ro,vec3 rd,inout vec3 p,inout float t,inout bool hit){
	hit = false;
    p = ro;
    t = 0.;
    vec2 d;
    for(int i = 0; i < 180 ; i++){
    	d = map(p);
                    
        
        if(d.x < 0.001){
        	hit = true;
            break;
        }
        t += d.x;;
    	p = ro + rd*t;
    }
	return d;
}

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = cross(vec3(0,1,0), dir);
	vec3 up = cross(dir, right);
    float fov = 0.9;
	return normalize(dir + right*uv.x*fov + up*uv.y*fov);
}
vec3 getNormal(vec3 p){
	vec2 t = vec2(0.001, 0.);
    return normalize(map(p).x - vec3(
    	map(p-t.xyy).x,
    	map(p-t.yxy).x,
    	map(p-t.yyx).x
    ));
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. - dot(uv,uv)*0.3;
    
    vec3 col = vec3(0);

    vec3 ro = vec3(0);
    ro.z += mx*1. ;
    ro += path(ro.z);
    vec3 lookAt = ro + vec3(0,0,4);
    lookAt += path(lookAt.z);
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    ro.xyz += rd * texture(image30, (uv+ TIME)*20.).xyz*0.7;
    
    NOISE = valueNoise(smoothTime*.5);
    ro.xy += NOISE*(0.025+syn_BassLevel*0.125)*Shake;
    float t; bool hit;
    vec3 p;
    
    vec2 d = march(ro, rd, p, t, hit);
    
    if(hit){
        vec3 n = getNormal(p);
        float fres = pow(1.0 - max(dot(n, -rd), 0.), 4.);
        
        col += fres*1.;
        
        
    	col *= glow*(0.021 + sin(smoothTimeB + t*0.6)*0.01);
    }
    
    
    col += glowF*0.1*smoothstep(0.,1., t*(0.01));
        
    col *= 0.5+syn_HighLevel;
    //col = mix(col, vec3(0.1,0.4,0.3),pow(smoothstep(0.,1.,t*0.05), 2.)*0.56);
    col = max(col, 0.);
    col = pow(col,vec3(0.4545));
    
    uv *= 0.8;
    col *= 1. - dot(uv,uv);
    //col *= 1. - t*0.2;    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// Fork of "Day 48" by jeyko. https://shadertoy.com/view/WtcSW4
// 2020-02-28 11:33:26

// HEX FUNCTION FROM BigWIngs !! he has a tutorial on this kind of tiling

// There's too much noise on fullscreen

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 100.;
    float scale = 0.00 + dot(uvn,uvn)*0.5;
    float chromAb = dot(uvn,uvn)*4.*(2.0+syn_BassLevel*2.);
    vec2 offs = vec2(0) + texture(image30, uv + TIME*4.).xz*0.001;
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.99;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.x += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.y += texture(BuffA, target).y;
    	radial.z += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    
    
    
    radial /= steps;
    
    fragColor = texture(BuffA, uv)*0.1 + radial*1.5;
    fragColor *= 0.97;
    //fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.5);
    //fragColor *= 1. - dot(uvn,uvn)*2.;
    
    
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