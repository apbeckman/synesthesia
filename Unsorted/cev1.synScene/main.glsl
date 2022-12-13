vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

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
    p.z = pmod(p.z, 10.);

    
    
    for(int i = 0; i < 4+iters; i++){
    	p = abs(p);
        p.x -= 2.;
        p.xy *= rot(0.25*pi);
        
        p.t -= 1.;
        p.z += 0.2;
        
        
    }
    vec3 glown = vec3(0);
    
    //float dB = length(p) - 0.4;
    
    p = abs(p);
    p.xy -= 0.3;
    
    p = abs(p);
    p.y -= 0.2+(sin(0.0125*smoothTimeC)*tog);
    p = abs(p);
    p.x -= 0.5+test;
    p = abs(p);
    p.x -= 0.5+(sin(0.0125*smoothTimeC)*tog);
    
    float attd = pow(abs(sin(j.z*0.5 + smoothTimeB*0.4)), 10.);; 
    float atte = pow(abs(sin(j.z*0.5 + smoothTimeB*.124)), 100.)*attd; 
    
    vec3 q = p;
    
    q.x -= 0.4;
    q = abs(q);
    q.x -= 0.4;

    q.y -= 0.4;
    
    vec3 sz = vec3(0.02,0.02+test2,20.5+test2);
    float dB = sdBox(q, sz)*1.2;
    
    d = dmin(d, vec2(dB, 2.));
    float att = pow(abs(sin(j.z*0.1 + smoothTimeB*0.5 + sin(j.x*4.)*0.2)), 400.);
    
    //Beam Colors
    //glow += 0.9/(0.001 + dB*dB*2000.)*vec3(1.8 + attd*8.,0.9,0.7)*att;
    glown += 0.9/(0.0001 + dB*dB*10.)*vec3(.8 + attd*8.,0.9,02.7)*att;
    
    p.z -= 1.5;
    
    p = abs(p);
    p.z -= 0.2;
    //p.xz *= rot(0.5);
    float dC = sdBox(p, vec3(0.02,2000.7,0.02));
    d = dmin(d, vec2(dC, 2.));
    float attb = pow(abs(sin(p.x*0.4 + smoothTimeB - 1.)), 10.);
    //glown += 7.9/(0.04 + dC*dC*2000.)*vec3(1.,1.,1.7)*attb;
	glown += 2.9/(0.001 - - atte*2. + dC*dC*400.)*vec3(1.,1.,1.7)*attb;

    //glow += 2.9/(0.001 + dC*dC*400.)*vec3(1.,1.,1.7)*attb;
    
    p.x -= 0.4;
    p = abs(p);
    p.x += 0.1;
    p.y -= 0.2;
    p.xy *= rot((-0.25+basshits*0.1)*pi);
    p.z -= 1.;
    //p -= 0.4;
    p.x -= 0.3;
    
    float dD = sdBox(p, vec3(0.02,1.7+basshits,0.02+basshits*0.1));
    d = dmin(d, vec2(dD, 2.));
    float attc = pow(abs(sin(p.y*0.24 + smoothTimeB + 4.)), 10.);
    //glow += 10.9/(0.01 + dD*dD*2000.)*sin(vec3(0.1,0.8,0.7) + vec3(0,0,attd*2.))*attc;
    //glow += 0.7/(0.001 + dD*dD*100.)*sin(vec3(0.1,0.1,0.9) + vec3(0,0,attd*2.))*attc;
    
    //Geometry Colors
    //glown += 0.7/(0.0001 + dD*dD*100.)*sin(vec3(0.1,1.1,0.9) + vec3(0,0,attd*1.))*attc;
    //glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.,1. - atte*0.6,1.9) + vec3(0,0,attd*1.))*attc*vec3(1,1,1.5)*0.8;
    //glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.,0.4 - atte*0.6,1.9) + vec3(0,0,attd*1.))*attc*vec3(1,1,1.5)*0.8;
    glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.61,0.4 - atte*0.6,1.1) + vec3(0,0,attd*0.))*attc*vec3(1,1,1.)*0.8;
    //glown += 0.7/(0.001 + dD*dD*(60. - atte*59.))*sin(vec3(0.1,0.4 - atte*0.6,.1) + vec3(0,0,attd*0.))*attc*vec3(1,1,1.)*0.8;
    
    d.x = abs(d.x) + 0.002;
    
    
    
    glow += glown;
    
    d.x *= 0.4;
    
    return d;
}
int it = 0;

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e7);
    
    t = 0.; hit = false; p = ro;
    
    for(it = 0; it< (60+FOV) + min(FRAMECOUNT,0); it++){
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
    lookAt.xy += _uvc*FOV*PI;
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
    
    ro.z += (smoothTime);
    
    vec3 lookAt = vec3(0);
    
    lookAt.z = ro.z + 2.;
    vec3 rd = getRd(ro, lookAt, uv);
        rd.xy = _rotate(rd.xy, PI*Rotate);

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
    float chromAb = pow(length(uv - 0.5),1.5)*0.4+basshits;
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