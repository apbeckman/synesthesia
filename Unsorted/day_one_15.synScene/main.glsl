vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 



#define pi acos(-1.)

#define tau (2.*pi)
#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))

#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))

vec3 glow = vec3(0);



float minRadius2 = 0.9;
float fixedRadius2 = 5.7 ;
float foldingLimit = 1.3;

int Iterations = 6;
float Scale = 1.9;

void sphereFold(inout vec3 z, inout float dz) {
	float r2 = dot(z,z);
	if (r2<minRadius2) { 
		// linear inner scaling
		float temp = (fixedRadius2/minRadius2);
		z *= temp;
		dz*= temp;
	} else if (r2<fixedRadius2) { 
		// this is the actual sphere inversion
		float temp =(fixedRadius2/r2);
		z *= temp;
		dz*= temp;
	}
}

void boxFold(inout vec3 z, inout float dz) {
	z = clamp(z, -foldingLimit, foldingLimit) * 2.0 - z;
}


#define pmod(p,j) mod(p,j) - 0.5*j

float map(vec3 z, float t){
	float d = 10e7;
	vec3 p = z;
    z.z = pmod(z.z, 10.);
    
    
    for(int i = 0; i < 4;i ++){
    	z = abs(z);
        
        z.xy *= rot(0.125*pi);
        //z.t -= 0.2;
        //z.z -= 0.3;
    }
    
    
    vec3 q = vec3(z);
    
    vec3 j;
    float jdr;
    
	vec3 offset = z;
	float dr = 1.;
	for (int n = 0; n < Iterations; n++) {
		boxFold(z,dr);       // Reflect
		sphereFold(z,dr);    // Sphere Inversion
        
        if(n == 2){
            j = z;
            jdr = dr;
        }
 		
                z=Scale*z + offset;  // Scale & Translate
                dr = dr*abs(Scale)+1.0;
	}
    
    
    //z = abs(z);
    //z.y -= 10.4;
    
	float r = length(z);
    
    
    
    d = r/abs(dr);
    
    d *= 0.7;
    
    d += smoothstep(1.,0.,t*0.75)*0.15;
    
    float db = length(j)/abs(jdr);;
    
    glow += 0.5/(0.6 + pow( (abs(d) + 0.001)*0.7,2.)*800000.)*0.9;
    
    db += 0.001;
    
    
    float att = pow(abs(sin(p.z + smoothTimeB*0.5 + length(p.xy))),50.);
    glow -= 0.92/(0.04 + pow( (abs(db) + 0.001)*0.7,2.)*16000.)*vec3(0.5,0.9,1.4)*att;
    
    
    
    d = min(d, db);
        
    //d *= 0.7;
    d = abs(d) + 0.00125;
    
    //glow -= 0.01/(0.1 + d*d*10.)*vec3(0.7,0.4,0.8);
    //glow += 0.12/(0.001 + d*d*4000.);
    
    
    return d;
    
    
    //float da = length(z.xyz)/q.w - 0.01;
    
        
    
    float sc = 0.5;
    //d = min(d,da);
    d *= sc;
    d += smoothstep(1.,0.,t*.5)*0.7;
    d = abs(d) + 0.003;
    
    vec3 c = vec3(1,1.,1.);
    return d;
}
vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0.,1.,0), dir));
	vec3 up = normalize(cross( dir, right));
	return normalize(dir + right*uv.x + up*uv.y);
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    
    
    vec3 col = vec3(0.9,0.6,0.4);

    vec3 ro = vec3(0);
    ro.z += smoothTime*0.76;
    
    //ro.xz -= 4.6;
    
    float T = smoothTime*1./tau + pi*0.25;
    ro.xy += _rotate(vec2(cos(T), sin(T))*(0.7+Radius), 0.+smoothTime*0.1);
    
    vec3 lookAt = vec3(0.001);
    
    lookAt.z = ro.z + 4.;
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    float d;
    vec3 p = ro; float t = 0.; bool hit = false;
    
    for(int i = 0; i < 120; i++){
    	d = map(p, t);
        if(d < 0.0001){
        	hit = true;
            //break;
        }
		t += d;
    	p = ro + rd*t;
    }
    
    
    col -= glow*0.07;
    
    col = max(col, 0.);
    
    col = pow(col, vec3(1. + dot(uv,uv)*1.));
    
    //col = smoothstep(0.,1.,col);
    
    col = pow(col, vec3(0.454545*(1.0-(highhits*0.75))));
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}