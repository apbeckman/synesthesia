float smin(float a, float b, float k) {
    float h = max(k - abs(a-b), 0.) / k;
    return min(a, b) - h*h*h*k*1./6.;
}



			//******** BuffA Code Begins ********


vec3 glow = vec3(0);
#define dmin(a, b) a.x < b.x ? a : b
#define PI acos(-1.)
#define tau (2.*PI)
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
//#define TIME (TIME + 3.6)
#define pal(a,b,c,d,e) (a + b*sin(c*d + e))

vec3 att = vec3(1);

float pModPolar(inout vec2 p, float repetitions) {
	float angle = 2.*PI/repetitions;
	float a = atan(p.y, p.x) + angle/2.;
	float r = length(p);
	float c = floor(a/angle);
	a = mod(a,angle) - angle/2.;
	p = vec2(cos(a), sin(a))*r;
	// For an odd number of repetitions, fix cell index of the cell in -x direction
	// (cell index would be e.g. -5 and 5 in the two halves of the cell):
	if (abs(c) >= (repetitions/2.)) c = abs(c);
	return c;
}
#define pmod(p,x) mod(p,x) - 0.5*x
vec4 valueNoise(float t){
	return mix(texture(image30,vec2(floor(t)/256.)),texture(image30,vec2(floor(t) + 1.)/256.), smoothstep(0.,1.,fract(t)));
}

vec2 map(vec3 p){
	vec2 d = vec2(10e7);
    float modD = 3.;
    float idz = floor(p.z/modD);
    p.z = pmod(p.z, modD);
    
    
    vec3 q = p;
    vec3 b = p;
    float idb = pModPolar(q.xy, polar_repeat);
    //b.xy *= rot(idz*0.5);
    pModPolar(b.xy, polar_repeat);
    
    b.x -= 0.8;
    // b.xy = mix(b.xy, _rotate(b.xy, TIME), 0.5);
    vec3 u = p;
    
    //u.xy *= rot(idz);
    float o = pModPolar(u.xy, 6.);
    u.x -= 1.;
    
    q.x -= 0.8;
    // q.xy = abs(q.xy);
    float dG = -u.x;
    d = dmin(d, vec2(dG, 3.));
    

    
    
    //u -= 0.1;
    u.y = abs(u.y);
    u.y -= 0.7;
    dG = length(u.xy) - 0.02;
    d = dmin(d, vec2(dG, 8.));
    
        //glow += 0.2/(0.01 + dG*dG*2.)*att;
	//q = abs(q);
    //q.xy *= rot(-1.5);
    //q.x -= 0.2;
    vec3 z = q;
    z = abs(z) - vec3(0.01,0.4,0.1);
    float dC = max(z.z, max(z.y, z.x));
    d = dmin(d, vec2(dC, 0.));
    z = q;
    z.x += 0.02;
    z = abs(z) ;
    z -= vec3(0.01,0.3,0.02);
    float dCb = max(z.z, max(z.y, z.x));
    d = dmin(d, vec2(dCb, 1.));
    
    
    
    z = b;
    z.y -= 0.2;
    z.x += 0.3;
    z.z += modD*0.10;
    z.xy *= rot(0.7 + sin(TIME*0.2 + idz*0.5));
    z = abs(z);
    z.zx *= rot(-0.1);
    z = abs(z) - vec3(0.01,0.5,0.04);
    float dD = max(z.z,max(z.x, z.y));
    
    d = dmin(d, vec2(dD, 5.));
    
    
    z = u;
    z.y -= 0.2;
    z.x += 0.3;
    z.z -= modD*0.25;
    
    z = abs(z);
    z.xz *=rot(0.25*PI);
    z.x -= 0.4;
    z.xy *=rot(0.25*PI);
    
    //z.x -= 0.08;
    z = abs(z) - vec3(0.07,0.1,0.04);
    float dDd = max(z.z,max(z.x, z.y));
    //d = dmin(d, vec2(dDd, 5.));
    
    
    vec4 a = valueNoise((idb + TIME*3. + idz*3.));
    
    //vec3 c = max(pal(0.7,1., vec3(3.7,0.3,0.6), 0.6,4.4 + sin(TIME) + sin(idz * idb)*0.2), 0.);
    vec3 c = max(pal(0.7,1., vec3(3.,0.3,0.1), 0.6,4.4 + sin(smoothTimeC*0.3) + idz + sin(idz * idb)*0.2), 0.1);
    
    glow += pow(smoothstep(0.,1.,a.z*1.5), 20.)*1.5/(0.005 + dCb*dCb*(90. - a.x*20.))*att*c* pow(smoothstep(1.,0.,length(q.y*1.)), 5.);
    //glow += pow(smoothstep(0.,1.,a.z*1.5), 20.)*1.5/(0.005 + dCb*dCb*(100. - a.x*20.))*att*c;
    //glow += pow(smoothstep(0.,1.,a.z*1.5), 20.)*1.5/(0.005 + dCb*dCb*(90. - a.x*20.))*att*c* (smoothstep(1.,0.,length(q.y*1.6)));
    //glow += pow(smoothstep(0.,1.,a.z*1.5), 20.)*0.01/(0.0004 + dCb*dCb*dCb*dCb*(50. - a.x*20.))*att*c;
    
    //glow += pow(smoothstep(0.,1.,a.z*1.5), 20.)*exp(-dCb*30.)*att*c*30.;
    //glow += pow(smoothstep(0.,1.,a.z*1.5), 20.)*0.6/(0.004 + dCb*dCb*(20. - a.x*20.))*att*c;
    
    return d;
}
float dith;
vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e7);

    p = ro; t = 0.; hit = false;
    for(int i = 0; i < 230 ; i++){
    	d = map(p);
        d.x *= dith;
        
    	//glow += exp(-d.x*20.);
        if(d.x < 0.002){
        	hit = true;
            break;
        }
        
        t += d.x;
        p = ro + rd*t;
    }
    
    
    return d;
}

vec3 getNormal(vec3 p){
	vec2 t= vec2(0.001,0);
	return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
    float inv_screen = mix(1.0, -1.0, invert_screen);

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    //
    // uv = mix(uv, uv+_uvc, fov);
	uv *= 1. - dot(uv,uv)*0.14;
    uv = mix(uv, uv+_uvc*PI, fov);
    float x_mirror_screen = smin(uv.x * -1, uv.x, 0.0125);
    float y_mirror_screen = smin(uv.y * -1, uv.y, 0.0125) ;
    vec2 smxy = vec2(x_mirror_screen, y_mirror_screen)*inv_screen;
    vec2 sm = vec2(screen_mirror_x, screen_mirror_y);
    //rd.yz *= rot(TIME);
    uv.xy = vec2(mix(uv.x, smxy.x, sm.x), mix(uv.y, smxy.y, sm.y));
    
    uv.xy *= rot((TIME - 3.6)*0.1);
    
    vec3 col = vec3(0);

    dith = mix(0.8,1., texture(image30, 20.*uv*256.).x);
    vec3 ro = vec3(0);
    
    ro.z += bass_time*1.5;
    
    vec3 rd = normalize(vec3(uv,2.));
    rd.xz = _rotate(rd.xz, lookxy.x*PI);
    rd.yz = _rotate(rd.yz, lookxy.y*PI);
    vec3 p; float t; bool hit;
    float side = 1.;
    float tA;
    
    for(int i = 0; i < 3; i ++){
    	vec2 d = march(ro, rd, p, t, hit);
    	vec3 n = getNormal(p);
        
        vec3 ld = normalize(vec3(1));
        vec3 h = normalize(ld - rd);
        
        float diff = max(dot(n, ld), 0.);
        float spec = pow(max(dot(n, -h), 0.), 10.);
        float fres = pow(1. - max(dot(n, -rd), 0.), 5.);
        
        if(i == 0){
        	tA = t;
        }    
        if(d.y == 5.){
        	//col += fres*0.1*att*(glow);
            col += fres*0.06*att*(glow);
        }
        if(d.y == 8.){
        	col += fres*0.02*att*(glow);
            //col += fres*20.*att;
        }
        if (d.y == 3.){
        	rd = reflect(rd, n);
            att *= vec3(0.6,0.8,0.8)*0.2;
            col += spec*0.04*att;
            ro = p + n*0.2;
        } else {
        	break;
        }
    }
    
    glow *= 1.0+highs*10.;
    
    col += glow*0.001;
    
    col = mix(col, vec3(0.4,0.4,0.7)*0.004, pow(smoothstep(0.,1.,tA*0.013), 1.6));
    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// radial blur and chromatic abberation in this buffer
// thx iq for pallette and hg-sdf for polarMod


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 30.;
    float scale = 0.00 + pow(length(uv - 0.5),4.)*0.25;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.)*3.7;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.97;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    float ss = smoothstep(0.,1.,dot(uvn,uvn)*3.);
    fragColor = radial*1.; 
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.4);
    fragColor *= 18.;
    fragColor = pow(fragColor, vec4(0.4545));
    fragColor *= 1. - dot(uvn,uvn)*2.;
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