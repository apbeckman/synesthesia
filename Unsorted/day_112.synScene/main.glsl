

			//******** BuffA Code Begins ********
float smin(float a, float b, float k) {
    float h = max(k - abs(a-b), 0.) / k;
    return min(a, b) - h*h*h*k*1./6.;
}

#define pi acos(-1.)

#define tau (2.*pi)

#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define dmin(a,b) a.x < b.x ? a : b

#define pmod(p,z) mod(p,z) - 0.5*z

#define pal(a,b,d,e,f) ((a) + (b)*sin((d) * (e) + (f)))

float sdBox(vec3 p, vec3 s){
	p = abs(p) - s;
    return max(p.x,max(p.y,p.z));
}

vec3 glow = vec3(0);

vec2 map(vec3 p){
	vec2 d = vec2(10e7);
    vec3 j = p;

    p.xy = (mix(p.xy,p.xy+ PI*_uvc*2., fov));
    vec2 m = vec2(x_mirror, y_mirror);
    float mirror_x = smin(-p.x, p.x, 0.25);
    float mirror_y = smin(-p.y, p.y, 0.25);
    vec2 mxy = vec2(mirror_x, mirror_y);
    p.xy = vec2(mix(p.x, mxy.x, m.x), mix(p.y, mxy.y, m.y));

    p.z = pmod(p.z, 10.);
    // p.xy = _rotate(p.xy, rotate*PI);
    p.xz = _rotate(p.xz, rotatexy.x*PI);
    p.yz = _rotate(p.yz, rotatexy.y*PI);
    
    
    for(int i = 0; i < 5; i++){
    	p = abs(p);
        p.x -= 2.;
        p.xy *= rot(0.25*pi);
    //    p.xy = _rotate(p.xy, 0.1*sin(TIME));
        // p.xy = _rotate(abs(0.75-2*p.xy), xtime*0.1);
        // p.xz = _rotate(p.xz, xtime*0.05);
        p.t -= 1.;
        p.z += 0.2;
        // p.xz = _rotate(p.xz, .2*_ncos(xtime*0.125));
        
        
    }
    vec3 glown = vec3(0);
    
    //float dB = length(p) - 0.4;
    
    p = abs(p);
    p.xy -= 0.3;
    
    for(int i = 0; i < 2; i++){
        p.x = smin(p.x, -p.x, 0.12);
        p.y = smin(p.y, -p.y, 0.12);
        // p.xy = _rotate(p.xy, 1.0 + 0.5*sin(TIME));
        p.xz = _rotate(p.xz, PI*0.25+.1*cos(xtime*0.125)+PI*morph_1.x);
        p.yz = _rotate(p.yz, PI*0.25+.1*sin(xtime*0.125+PI*morph_1.y));
        p.z = smin(p.y, -p.z, 0.12);
    }    
    p = abs(p);
    p.y -= 0.2;
    p = abs(p);
    p.x -= 0.5;
    p = abs(p);
    p.x -= 0.5;
    
    float attd = pow(abs(sin(j.z*0.5 + TIME*0.4)), 10.);; 
    float atte = pow(abs(sin(j.z*0.5 + TIME*10.4)), 100.)*attd; 
    
    vec3 q = p;
    for(int i = 0; i < 4; i++){
    q = abs(q);
    q.xy = _rotate(q.xy,.5+.2*sin( -ytime*0.1));
    }    
    q.x -= 0.4;
    q = abs(q);
    q.x -= 0.4;

    q.y -= 0.4;
    
    vec3 sz = vec3(0.02,0.02,20.5);
    float dB = sdBox(q, sz)*1.2;
    
    d = dmin(d, vec2(dB, 2.));
    float att = pow(abs(sin(j.z*0.1 + smoothTimeC*0.0125 + sin(abs(j.x)*4.)*0.2)), 200.);
    
    //glow += 0.9/(0.001 + dB*dB*2000.)*vec3(1.8 + attd*8.,0.9,0.7)*att;
    glown += 0.9/(0.0001 + dB*dB*10.)*vec3(2.8 + attd*8.,0.9,0.7)*att;
    
    p.z -= 1.5;
    
    p = abs(p);
    p.z -= 0.2;
    //p.xz *= rot(0.5);
    float dC = sdBox(p, vec3(0.02,2000.7,0.02));
    d = dmin(d, vec2(dC, 2.));
    float attb = pow(abs(sin(p.x*0.4 + TIME - 1.)), 10.);
    //glown += 7.9/(0.04 + dC*dC*2000.)*vec3(1.,1.,1.7)*attb;
	glown += 2.9/(0.001 - - atte*2. + dC*dC*400.)*vec3(1.,1.,1.7)*attb;

    //glow += 2.9/(0.001 + dC*dC*400.)*vec3(1.,1.,1.7)*attb;
    
    p.xy = _rotate(p.xy, TIME*0.1);
    p.x -= 0.4;
    p = abs(p);
    p.x += 0.1;
    p.y -= 0.2;
    p.xy *= rot(-0.25*pi);
    p.z -= 1.;
    //p -= 0.4;
    p.x -= 0.3;
    
    float dD = sdBox(p, vec3(0.02,1.7,0.02));
    d = dmin(d, vec2(dD, 2.));
    float attc = pow(abs(sin(p.y*0.4 + TIME + 4.)), 10.);
    //glow += 10.9/(0.01 + dD*dD*2000.)*sin(vec3(0.1,0.8,0.7) + vec3(0,0,attd*2.))*attc;
    //glow += 0.7/(0.001 + dD*dD*100.)*sin(vec3(0.1,0.1,0.9) + vec3(0,0,attd*2.))*attc;
    
    
    //glown += 0.7/(0.0001 + dD*dD*100.)*sin(vec3(0.1,1.1,0.9) + vec3(0,0,attd*1.))*attc;
    //glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.,1. - atte*0.6,1.9) + vec3(0,0,attd*1.))*attc*vec3(1,1,1.5)*0.8;
    //glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.,0.4 - atte*0.6,1.9) + vec3(0,0,attd*1.))*attc*vec3(1,1,1.5)*0.8;
    glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.1,0.4 - atte*0.6,1.1) + vec3(0,0,attd*0.))*attc*vec3(1,1,1.)*0.8;
    
    
    d.x = abs(d.x) + 0.002;
    
    
    glow += glown;
    
    d.x *= 0.4;
    
    return d;
}
int it = 0;

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e7);
    
    t = 0.; hit = false; p = ro;
    
    for(it = 0; it< 60 + min(FRAMECOUNT,0); it++){
    	d = map(p);
        
        if(d.x < 0.001){
        	hit = true;
        	break;
        }
    	t += d.x;
        p = ro + rd*t;
    }
    
    
    return d;
}

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
    vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
    return normalize(dir + (right*uv.x + up*uv.y)*0.6);
}

vec3 getNormal(vec3 p){
    vec2 t = vec2(0.001,0.);
    return -normalize(vec3(
    	map(p - t.xyy).x - map(p + t.xyy).x,
    	map(p - t.yxy).x - map(p + t.yxy).x,
    	map(p - t.yyx).x - map(p + t.yyx).x
    ));
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0);
    
	uv *= 1. - dot(uv,uv)*0.06;
    
    vec3 ro = vec3(0);
    
    ro.z += bass_time*3.;
    
    vec3 lookAt = vec3(0);
    
    lookAt.z = ro.z + 2.;
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    vec3 p; float t; bool hit;
    
    vec2 d = march(ro, rd, p, t, hit);
    
    if (d.x < 0.001){
    	vec3 n = getNormal(p);
        vec3 l = normalize(vec3(1));
        float dnh = max(dot(n,normalize(l - rd)),0.);
        float diff = max(dot(n,l),0.);
        float spec = pow(max(dot(n,normalize(l - rd)),0.), 20.);
        float fres = pow(1. - dnh, 5.);
    
        
        //col *= fres;
        
    }
    
    //col += float(it)*0.003;
    
    
    col += glow*0.001;
    col = mix(col, vec3(0.1,0.1,0.16 + (uv.x + 0.5)*0.1)*0.02, smoothstep(0.,1.,t*0.04 ));
    
    
    
    //col += pal(0.5,0.6,vec3(0.7,0.3,0.7), 0.7,0.6)*smoothstep(1.,0.,length(glow)*0.3);;
    col = max(col, 0.);
    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


// Fork of "Day 98" by jeyko. https://shadertoy.com/view/wdXczX
// 2020-04-06 10:18:27



// radiual blur in this buffer

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
	fragColor = vec4(0);
    // Radial blur
    float steps = 26.;
    float scale = 0.00 + pow(dot(uvn,uvn),1.1)*0.04;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.5)*0.4;
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
    
    //fragColor.b *= 0.97 + dot(uvn,uvn)*0.4;
    //fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.6);
    
    fragColor.t *= 1.  - smoothstep(0.,1.,dot(uvn,uvn))*0.;
    
    
    //fragColor.b *= 1. + uv.x*0.02;
    //fragColor.g *= 1. + uv.t*0.05;
    fragColor *= 1. - dot(uvn,uvn)*2.;
    
    fragColor = max(fragColor, 0.);
    fragColor = pow(fragColor, vec4(0.4545 ));
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