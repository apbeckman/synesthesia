vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 




#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define mx (10.*iMouse.x/RENDERSIZE.x)
#define my (10.*iMouse.y/RENDERSIZE.x)


vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0,1,0), dir));
	vec3 up = normalize(cross(dir, right));
	return normalize(dir + uv.x*right + uv.y*up);
}
vec3 glowA = vec3(0);
vec3 glowB = vec3(0);
vec3 glowC = vec3(0);
vec3 glowD = vec3(0);


#define ITERS 4
#define ZOOM (1. + my)
#define pmod(p, x) mod(p,x) - 0.5*x
float[] WS = float[ITERS](1.,1.,1.,1.);

float id;
float sd(vec3 q){
	#define repD 3.97
    id = floor(q.z/repD);
    q.z = pmod(q.z, repD);
    q.xy -= 2.;
    q.xy = pmod(q.xy, 4.);
    vec4 p = vec4(q.xyz, 1.);
    vec4 c = vec4(0.9,0.77,0.4,0.2);
    vec4 u = vec4(0.4,0.54,0.7,0.4);
    for(int i = 0; i < ITERS; i++){
        p.xyz = abs(p.xyz) - c.xyz;
        float dpp = dot(p.xyz,p.xyz);
        p=p*(1.5 + u)/clamp(dpp,.4,1.)- mix(p,c,0.9)*0.9;
        if(i==1)
			p.xy += 0.4;
        p.xy *= rot(-0.1 + sin(id)*0.9);
        
        p.xz *= rot(0.9);
        //if(i == 2)
        //	p.y += 0.9;
        WS[i] = p.w;
    }
    //w = p.z;
    p.xyz = abs(p.xyz);
    
    p.xz *= rot(0.2);
    float fr =  max(p.x - 0., max(p.y - 2.4 + sin(id)*0.7, p.z - 2.9))/p.w;
    float d = fr;
        
	return d*0.5;
}

vec2 map(vec3 p){
	vec2 d = vec2(10e5);
    
    d.x = min(d.x, sd(p));
    //d.x = min(d.x, length(p) - 0.6);
    return d;
}

vec2 march(vec3 ro, vec3 rd,inout vec3 p,inout float t,inout bool hit){
	vec2 d = vec2(10e5);
    hit = false;
    p = ro;
    t = 0.;
    for(int i = 0; i < 60; i++){
    	d = map(p);
        glowA += exp(-d.x*20.)*WS[0];
        glowB += exp(-d.x*05.)*WS[1]*0.5;
        glowC += exp(-d.x*10.)*WS[2];
        glowD += exp(-d.x*50.)*WS[3]*0.7;
        if (d.x < 0.001){
        	hit = true;
            break;
        }
        t += d.x;
        p = ro + rd*t;
    }
    
    return d;
}

vec3 getNormal(vec3 p){
	vec2 t = vec2(0.0002, 0);
	return normalize(
    	map(p).x - vec3(
        	map(p - t.xyy).x,
        	map(p - t.yxy).x,
        	map(p - t.yyx).x
        )
    );
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec4 S = texture(image11, vec2(0.01,0.0));
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. + dot(uv,uv)*1.;
    
    vec3 col = vec3(0);
    
    vec3 ro = vec3(0);
    //ro.x += sin(mx)*ZOOM;
    //ro.z += cos(mx)*ZOOM;
    float m = sin(smoothTime*0.125);
    float n = cos(smoothTime*0.125);

    ro.xy += vec2(n*0.125, m*0.125);
    ro.z += smoothTime + 14.;
    vec3 lookAt = vec3(0);
    lookAt.z = ro.z +1.;
    //lookAt.x += m*3.;
    vec3 rd = getRd(ro, lookAt, uv);
    
    vec3 p;float t; bool hit;
    
    vec2 d = march(ro, rd, p, t, hit);
    
    #define pal(a,b,c,d,e) (a + b*sin(c*d + e))
    if(hit){
		vec3 n = getNormal(p);
        vec3 lD = normalize(vec3(1));
        vec3 h = normalize(lD - rd);
        float diff = max(dot(n, lD),0.);
        float spec = pow(max(dot(n, h),0.), 20.);
        float fres = pow(1. - max(dot(n, -rd),0.), 5.);
    	//col += diff*fres*4.;
        
        vec3 CA = WS[0]*vec3(0.6,0.2,0.7)*0.01;
        vec3 CB = WS[1]*vec3(0.6,0.2,0.7)*0.01;
        vec3 CC = WS[2]*vec3(0.6,0.2,0.7)*0.01;
        //vec3 CD = WS[3]*vec3(0.1,0.9,0.4)*0.01;
        //vec3 CD = 0.05*WS[2]*pal(0.2, 0.6, vec3(0.8,0.4,0.7), 5.9 + sin(id)*0.2, 2.4);
        vec3 CD = 0.05*WS[2]*pal(0.2, 0.6, vec3(0.8,0.4,0.7), 5.9 , 0.2 + S.x);
        vec3 C = CA + CB + CC + CD; 
        
        col += diff * 0.2;
        col += pow(CD*0.9,vec3(2.));
    } else {
    
    }
    
    vec3 G = vec3(0);
    G += pow(glowD*0.004,vec3(1.))*pal(0.2, 0.6, vec3(0.8,0.4,0.7), 0.6, 2.4);
    G -= pow(glowC*0.005, vec3(1.1))*pal(0.2, 0.4, vec3(0.8,0.4,0.7), 9.6, 2.4);
    
    G += glowB*0.002 *pal(0.2, 0.6, vec3(0.8,0.4,0.7), 5.99 - sin(smoothTimeB*0.125), 2.4);
    G *= 1. + pow(S.x,5.)*0.2;
    
    col += G;
    uv.y *= 1.5;
    col *= 1. - dot(uv,uv)*0.1;
    col = mix(col, vec3(0.5,0.4,0.35)*0.3, smoothstep(0.,1.,t*0.1 - 0.1));
    col += (glowB*0.004)*(0.5+highhits*0.125);
    //col += glow*0.006;
    col = clamp(col, 0., 1.);
    col = pow(col, vec3(0.7));
    col = smoothstep(0.,0.94, col);

    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}