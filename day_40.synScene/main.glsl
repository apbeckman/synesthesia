//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

// It was going ok, but now I kind of hate it.


// --------- check out CUTEMODE, looks better --------- //
#define CUTEMODE 

#define mx (TIME*3.)


#define dmin(a,b) a.x < b.x ? a : b
#define pmod(p, x) mod(p, x) - x*0.5
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))


vec3 glow = vec3(0);

#define modDist 5.
vec3 u;
vec2 map(vec3 p){
	vec2 d = vec2(10e4);
	p.xz = pmod(p.xz, modDist);
    
    float pR = 4.;
    p.y = abs(p.y);
    p.y -= pR*1.1;
    vec3 q = abs(p) - pR;
    u = q;
    float dPyra = max(q.x, q.z);
    dPyra = max(dPyra, dot(q.xy + 0.5*pR,normalize(vec2(1))));
    dPyra = max(dPyra, dot(q.zy + 0.5*pR,normalize(vec2(1))));
    //d.x = min(d.x, length(p) - 0.5);
    d.x = min(d.x, dPyra);
    //d.x = min(d.x, length(p) - 0.5);
    d.x *= 0.22;
    return d;
}

vec3 getNormal(vec3 p){
    vec2 t = vec2(0.001,0.);
    float d = map(p).x;
	float a = map(p - t.xyy).x;
	float b = map(p - t.yxy).x;
	float c = map(p - t.yyx).x;
    return normalize(d - vec3(a,b,c));
}
void getNormalAndEdge(vec3 p,inout vec3 n,inout float edge){
    vec2 t = vec2(0.05,0.);
    n = getNormal(p);
    vec3 a = getNormal(p + t.xyy);
    vec3 b = getNormal(p + t.yxy);
    vec3 c = getNormal(p + t.yyx);
    vec3 d = getNormal(p - t.xyy);
    vec3 e = getNormal(p - t.yxy);
    vec3 f = getNormal(p - t.yyx);
    edge = length(a - b - c );
    edge = length(a - b - c - d - e - f);
}

#define pal(a,b,c,d,e) (a + b*sin(c*d + e))
vec2 march(vec3 ro, vec3 rd,inout float t,inout bool hit){
	hit = false;
    t = 0.;
    vec2 d;
    vec3 p = ro + rd*1.;
    for(int i = 0; i < 250; i++){
    	d = map(p);
        //glow += exp(-d.x*15.)*pal(0.6 + sin(TIME),vec3(1.3,1.1,1.),vec3(4.7,2.2,3.), d.x*9., 10.);
        //glow += exp(-d.x*10.)*pal(0.5,vec3(1.3 ,1.1 ,1.),vec3(1.7+ sin(TIME*2.)*1.,2.2,3.), d.x*5., 10.);
        
        //glow += exp(-d.x*10.)*pal(0.5,vec3(1.3,1.1,1.),vec3(4.7,2.2,3.), d.x*9., 10.);
        glow += exp(-d.x*10.)*pal(0.5,vec3(1.3 ,1.1 ,1.),vec3(1.7+ sin(TIME*2.)*1.,2.2+highhits,3.), d.x*5., 10.);
        //glow += exp(-d.x*0.)*pal(0.5,vec3(0.5 ,0.5 ,0.5),vec3(1.7+ sin(TIME*1.)*1.,9.2,3.), d.x*5., 10.)*0.1;
        if(d.x < 0.0001){
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
	vec3 right = normalize(cross(vec3(0,1,0),dir));
	vec3 up = normalize(cross(dir, right));
    return normalize(dir + right*uv.x + up*uv.y);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0);

    vec3 ro = vec3(0);
    //ro.z += mx;
    ro.x -= smoothTime*4.;
    ro.y += 0.2;
    vec3 lookAt = ro + vec3(0,0,0);
    lookAt.x += sin(-2.);
    vec3 rd = getRd(ro, lookAt, uv);
    rd.yz *= rot(sin(smoothTime*0.1225)*0.1);
    rd.xy *= 1. - dot(uv,uv)*0.2;

    float t = 0.;
    bool hit = false;
    vec2 d = march(ro, rd, t, hit);
    vec3 n;
    float edge; 
    vec3 p =ro + rd*t;
    vec3 j = p;
    getNormalAndEdge( j*2. - vec3(0,2.5,0), n, edge);
    if(hit){
        float l = 1.;
        //n = 1.;
        //vec3 i = p;
        col += abs(sin(smoothTimeB*0.2)+0.5)*smoothstep(l*2.5,l,pow(edge*0.29,6.))*2.*fract(clamp(sin(smoothTimeB + p.z - cos(p.x) - p.y),0.,1.));
        //#define pal(a,b,c,d,e) (a + b*sin(c*d + e))
#define pi acos(-1.)
        float mR = 6.;
        float w = 0.2;
        u.y = abs(u.y);
        u = abs(u);
        u.xy *= rot(0.25*pi);
        u.zy *= rot(0.225*pi);
        u = abs(u);
        u.x -= 20.;
        float id = floor(u.y*mR);
        vec2 idZX = floor(p.xz/modDist);
        #ifdef CUTEMODE
        float m = abs(mod(u.y*mR,1.)-2.5);
        #else 
        float m = abs(mod(u.y*mR,1.)-.9);
        
        #endif
        col += smoothstep(w,w*.5 + sin(id)*1.4,m)*pal(0.5,vec3(1.9,0.7,0.1),vec3(0.7,0.2,0.8), 2. + idZX.x - idZX.y*2., id + idZX.x + TIME*0.9);
        //col += n;
    } else {
    	
    }
    
    
    col += glow*0.04;
    col = max(col, vec3(0.00));
    
    //col = pow(col, vec3(1.5,1. + sin(TIME)*0.2,1. - sin(TIME)*0.2));
    col *= 0.7;
    //col.b*= 1.2;
    
    col = mix(col, vec3(0.0,0,0.0),clamp(pow(abs(t)*0.11 - 4.9, 1.), 0., 1.));
    col = clamp(col, 0.00,1.);
    col = mix(col,1. - col,0.01);
    //col = smoothstep(0.,1.,col);
    //col = mix(col, vec3(0,0,0),clamp(pow(abs(p.y)*0.5 - 0.4, 1.), 0., 1.));
    //col = mix(col, vec3(0.04,0,0.04),clamp(pow(abs(p.y)*0.5 - 0.9, 1.), 0., 1.));
    //col = smoothstep(-0.05,1.,col);
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 




vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    //uvn.yz = _rotate(uvn.yz, lookXY.y*PI);
    uvn.xy = _rotate(uvn.xy, -1.0*lookXY.x*PI);

    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 20.;
    float scale = 0.00 + pow(length(uv - 0.5)*1.2,3.)*0.1;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.4)*3.1;
    vec2 offs = vec2(0);
    vec4 radial = vec4(0);
    for(float i = 0.; i < steps; i++){
    
        scale *= 0.91;
        vec2 target = uv + offs;
        offs -= normalize(uvn)*scale/steps;
    	radial.r += texture(BuffA, target + chromAb*1./RENDERSIZE.xy).x;
    	radial.g += texture(BuffA, target).y;
    	radial.b += texture(BuffA, target - chromAb*1./RENDERSIZE.xy).z;
    }
    radial /= steps;
    
    
    fragColor = radial*1.; 
    
    fragColor *= 1.3;
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.8);
    
    fragColor = max(fragColor, 0.);
    fragColor.xyz = pow(fragColor.xyz, vec3(0.8,1. + sin(TIME)*0.2,1. - sin(TIME)*0.2));

    //fragColor = pow(fragColor, vec4(0.4545 + dot(uvn,uvn)*1.7));
    fragColor *= 1. - dot(uvn,uvn)*1.8;
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