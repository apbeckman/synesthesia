vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** BuffA Code Begins ********

// Fork of "day 75" by jeyko. https://shadertoy.com/view/3tGXRK
// 2020-03-01 19:53:26

vec3 attenuation = vec3(1);

#define mx (smoothTime*0.5)

#define pmod(p,x) mod(p,x) - x*0.5

#define dmin(a,b) a.x < b.x ? a : b
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))


#define modDist 1.
#define HEIGHT 0.5

#define pi acos(-1.)
#define tau (2.*pi)
#define pal(a,b,c,d,e) ( (a) + (b)*sin(tau*((c)*(d) +(e)) ))
// The "Stairs" flavour produces n-1 steps of a staircase:
// much less stupid version by paniq
float fOpUnionStairs(float a, float b, float r, float n) {
	float s = r/n;
	float u = b-r;
	return min(min(a,b), 0.5 * (u + a + abs ((mod (u - a + s, 2. * s)) - s)));
}


float side = 1.;
vec3 glow = vec3(0);
vec3 glowB = vec3(0);

vec2 map(vec3 p){
	vec2 d= vec2(10e6);

    p.y += _noise((sin(p.z+TIME)) + abs(p.x)*4.  - .1*smoothTime)*0.2*Floor*(1.0+Floor*syn_BassPresence*0.2);
    
    
    vec3 q = p;
    
    vec2 id = floor(p.xz);
    p.xz = pmod(p.xz, modDist);
    vec3 j = p;
    
    float dBall = length(p) - 0.14;
    
    d = dmin(d, vec2(dBall, 0.));
    #define colW 0.1
    
    j.xz = abs(j.xz)- colW;
    vec3 k = j;
    j = abs(j);
    vec3 i = j;
    j.y -= HEIGHT*0.95;
    for (int i =0; i < 3; i ++){
        j = abs(j);
        
        j.x -= 0.05;
        j.xy *= rot(-0.25*pi + sin(smoothTimeC*0.25+ id.y*0.1)*0.75 - 0.34 );
        j.xz *= rot(-0.25*pi);
    }
    
    
    float dCol = min
        (max(j.x, j.z),max(k.x, k.z) );
    
    
    
    d = dmin(d, vec2(dCol, 1.));
    
    float dGoldA = min(length(j.xz),length(k.xz) ) - 0.01 ;
    
    d = dmin(d, vec2(dGoldA, 2.));
    
    i = abs(i);
    i.xy *= rot(0.4);
    i.xz *= rot(0.2);
    i = abs(i);
    i.xz *= rot(0.6*pi);
    
    dGoldA = length(i.xy)- 0.0123 ;
    d = dmin(d, vec2(dGoldA, 2.));
    

    p.xz = pmod(p.xz, 0.2);
    
    p.y -= 0.3;
    p = abs(p) - 0.001;
    
    float dGlow = min(
    	length(p.xy) - 0.1,
    	length(p.zy) - 0.1
    );
    
    //d = dmin(d, vec2(dGlow,9.));
    //float attenC = pow(abs(sin(q.z*0.1  + sin(q.x*0.2 + TIME)*2. + sin(q.y*3.)*0. + TIME*2.2)), 200.);
    float attenC = pow(abs( sin(q.x +smoothTimeB*0.2 + sin(q.z*0.1 + smoothTimeB*0.2) + q.z) ), 40.);

    //float attenC = 0.2;
    
    vec3 col = pal(0.2,0.8 ,vec3(0.1 + pow(abs(sin(smoothTimeB*0.2*1.)), 40. )*0.005,2.2,0.3),0.5 + sin(smoothTimeB*0.2)*0.005,0.5 - attenC*0.6);

    glowB += exp(-dGlow*20.) * col*attenuation*4.*attenC;    
    
    
    q.xz = pmod(q.xz, 0.2);
    
    float dCeilings = q.y + HEIGHT;
    dCeilings = min(dCeilings,-q.y + HEIGHT);
    
    dCeilings += max(exp(-q.x*100.), exp(-q.z*100.))*0.0000001;
    
    d = dmin(d, vec2(dCeilings, 0.));
    
    
    
    //d = dmin(d, vec2(dGlow, 5.));
    
    d.x *= 0.5;
    
    //d.x *= 0.5;
    return d;
}


vec2 march(vec3 ro, vec3 rd,inout vec3 p,inout float t,inout bool hit){
	vec2 d = vec2(10e6);
    p = ro; t = 0.; hit = false;
    
    for(int i = 0; i < 150 ; i++){
    	d = map(p);
        d.x *= side;
        glow += exp(-max(d.x, 0.)*80.)*attenuation*1.;
        
        
        if(d.x < 0.001 || t > 20.){
        	hit = true;
            break;
        }
        t+=d.x;
        p += rd*d.x;
    }
    
    
    return d;
}

vec3 getNormal(vec3 p){
	vec2 t= vec2(0.0001,0);
    return normalize(vec3(
    	map(p - t.xyy).x - map(p + t.xyy).x ,
    	map(p - t.yxy).x - map(p + t.yxy).x ,
    	map(p - t.yyx).x - map(p + t.yyx).x 
    ));
}

vec3 getNormala(vec3 p){
	vec2 t= vec2(0.00001,0);
    return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x 
    ));
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y,
        d;

    uv *= 1. + length(uv)*0.2;
    
    vec3 col = vec3(0),
    	 ro = vec3(0),
    	 rd = normalize(vec3(uv,1.)),
         p;

    ro.y += sin(smoothTime*0.5)*0.0124;    
    ro.z += mx;
    
    
    rd.xy *= rot(sin(smoothTimeC*0.1)*0.0125);
    float t; bool hit;
    float tB = 0.;
    
    for(int i = 0; i < 2 + min(FRAMECOUNT, 0); i++){
    	d = march(ro, rd, p,t, hit);
    	
        if(i == 0)
            tB = t;
        
        vec3 n = getNormal(p)*side;
        
        //float diff = max();
        
        if(d.y == 0.){
        	//col += 0.1;
            /*
            if(i == 0)
                glow *= 0.5;
            attenuation *= 0.5;*/
        }
        if(d.y == 2.){
            vec3 c = vec3(0.7,0.6,0.5)*4.;
            if(i == 0)
                glow *= c*1.;
            
            glow *= c*1.;
        	attenuation *= c*0.6;
        }
        if(d.y == 1.){
        	//col += 0.1;
            //vec3 c = vec3(0.1,0.6,0.9);
            vec3 c = vec3(0.11,0.1,0.1);
            if(i == 0)
                glow *= c;
            glow *= c;
        	attenuation *= c*3.;
            //side *= -1.;
        	rd = refract(rd, n, 0.4 + sin(p.x*40. + p.y*40. + smoothTime*0.2)*0.012);
        	//rd = reflect(rd, n);
        } else {
        	rd = reflect(rd, n);
            
        }
        
        
        ro = p + rd*0.05;
        
        
        attenuation *= 0.7;
    }
    
    
    glow = pow(glow, vec3(1.1))*0.02;
    
	col += glow*0.02;
	col += glowB*0.01;
    
    col = clamp(col, 0., 10.);
    
    
    col = mix(col, vec3(0.7), smoothstep(0.,1.,tB*0.05));
    
    
    uv *= 0.8;
    col *= 1. - dot(uv,uv);
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 



// SHOUTOUT TO LUNA, SHE HAS AWESOME GOLD STUFF shadetoy.com/user/yx

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 10.;
    float scale = 0.00 + dot(uvn,uvn)*0.025;
    float chromAb = dot(uvn,uvn)*4.;
    vec2 offs = vec2(0) + texture(BuffA, uv + smoothTimeC*4.).xz*0.001;
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
    
    fragColor = texture(BuffA, uv)*0.8 + radial*1.5;
    
    fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.7);
    fragColor *= 0.6;
    fragColor = pow(fragColor, vec4(0.45));

    
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