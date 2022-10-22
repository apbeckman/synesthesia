

			//******** Common Code Begins ********


#define MIPLVL 1.

#define PI acos(-1.)
#define tau (2.*PI)
#define R RENDERSIZE.xy
#define T(uv) texture(iChannel0, (uv)/R)

vec4 gaussian(sampler2D chan, vec2 uv, vec2 RENDERSIZE, float mip){
    float st = 4.;
    vec3 t = vec3(st, 0., -st);
    vec4 C = vec4(0);
	#define TC(uv) texture(chan, (uv)/R, mip) 
    // don't do defines like this lol
    C += 4.*TC(uv);
    C += 2.*TC(uv - t.xy) + 2.*TC(uv + t.xy) + 2.*TC(uv - t.yx) + 2.*TC(uv - t.yx);
    C += 1.*TC(uv - t.xx) + 1.*TC(uv + t.xx) + 2.*TC(uv - t.xz) + 2.*TC(uv - t.xz);

    return C / 16.;
}

float DistributionGGX(vec3 N, vec3 H, float roughness)
{
    float a = roughness*roughness;
    float a2 = a*a;
    float NdotH = max(dot(N, H), 0.0);
    float NdotH2 = NdotH*NdotH;

    float nom   = a2;
    float denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return nom / max(denom, 0.001); // prevent divide by zero for roughness=0.0 and NdotH=1.0
}


float distributionTerm(float roughness, float ndoth) {
	float r2 = roughness * roughness;
	float d	 = (ndoth * r2 - ndoth) * ndoth + 1.0;
	return r2 / (d * d * PI);
}

float D_GGX(float NoH, float roughness)
{
	float a = roughness * roughness;
    float a2 = a * a;
    float nom = a2;
    float denom = (NoH * NoH * (a2 - 1.0) + 1.0);
	denom = PI * denom * denom;
    
    return nom / denom;
}


// ----------------------------------------------------------------------------
float GeometrySchlickGGX(float NdotV, float roughness)
{
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;

    float nom   = NdotV;
    float denom = NdotV * (1.0 - k) + k;

    return nom / denom;
}
// ----------------------------------------------------------------------------
float GeometrySmith(vec3 N, vec3 V, vec3 L, float roughness)
{
    float NdotV = max(dot(N, V), 0.0);
    float NdotL = max(dot(N, L), 0.0);
    float ggx2 = GeometrySchlickGGX(NdotV, roughness);
    float ggx1 = GeometrySchlickGGX(NdotL, roughness);

    return ggx1 * ggx2;
}
// ----------------------------------------------------------------------------
vec3 fresnelSchlick(float cosTheta, vec3 F0)
{
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}


float random(vec2 u){
	return fract(sin(u.y*4125.1 + u.x *125.625)*225.5235);
} 

float noise(vec2 p) {
	vec2 i = ceil(p);
    vec2 f = fract(p);
    vec2 u = f * f * (3. - 2. * f);
   	float a = random(i);
    float b = random(i + vec2(1., 0.));
    float c = random(i + vec2(0., 1.));
    float d = random(i + vec2(1., 1.));
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

float r31(vec3 u){
	return fract(sin(u.y*125.1 + u.x *125.125 + u.z*525.5215)*115.125235);
} 
float valueNoise(vec3 uv){
    vec3 id = floor(uv);
    vec3 fd = fract(uv);
    fd = smoothstep(0.,1., fd);
    
    float ibl = r31(id + vec3(0,-1,0));
    float ibr = r31(id + vec3(1,-1,0));
    float itl = r31(id + vec3(0));
    float itr = r31(id + vec3(1,0,0));
    
    
    float jbl = r31(id + vec3(0,-1,1));
    float jbr = r31(id + vec3(1,-1,1));
    float jtl = r31(id + vec3(0,0, 1));
    float jtr = r31(id + vec3(1,0, 1));
    
    
    float ibot = mix(ibl, ibr, fd.x); 
    float iup = mix(itl, itr, fd.x);
    float jbot = mix(jbl, jbr, fd.x);
    float jup = mix(jtl, jtr, fd.x);
    
    float i = mix(ibot, iup, fd.y);
    float j = mix(jbot, jup, fd.y);
    
    return mix(i, j, fd.z); 
}

float fbm(vec2 p) { 
	float s = .0;
	float m = .0;
	float a = .5;	
	for(int i = 0; i < 6; i++) {
		s += a * noise(p);
		m += a;
		a *= .5;
		p *= 2.;
	}
	return s / m;
}


			//******** BuffA Code Begins ********


#define mx (20.*_mouse.x/RENDERSIZE.x)
#define TIME (TIME + mx)

#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
#define mx (20.*_mouse.x/RENDERSIZE.x)

#define dmin(a,b) a.x < b.x ? a : b 
#define pmod(p, x) mod(p,x) - 0.5*x

vec3 getWallpaper( vec2 fragCoord )
{
    vec2 p = mod(fragCoord,1.) - 0.5;
	p *= 1.4;
    vec3 col = vec3(0);

	float d = 10e6;
    
    float s = 1.;
    float scale = 1.;
    
    for(int i = 0; i < 5; i++){
        p = -1.0 + 2.0*fract(0.5*p+0.5);
        float dpp = dot(p,p);

        p = abs(p);
        p += vec2(0.06,0.02);
        dpp = clamp(dpp,0.3,0.45);
    	//p /= dpp; 
    	//p = sin(p + t);
    
        float k = s/dpp;
        k = clamp(k,0.,4.);
        
        
		p     *= k;
		scale *= k;
        
    }
    p = abs(p) - 0.1;
    p = abs(p);
   	d = min(p.x, p.y)/scale; 
    
    col += smoothstep(0.01,0.,d);
    
    col = pow(col, vec3(4.));
    
    return col;
}

#define wallW 0.2
#define wallH 0.5


float wp = 0.;
float fw = 0.;
vec2 map(vec3 p){
	vec2 d = vec2(10e6, 1.);
    
    
    //p.x -= mix(0.,exp(abs(p.y)),smoothstep(0.,1.,abs(p.y)*4.))*0.1;
    float lWall = p.x + wallW;
    
    wp = getWallpaper(p.yz*0.25).x;
    
    wp = clamp(wp, 0., 1.);
    
    lWall += wp*0.01;
    
    fw = fbm(p.yz*7.);
    
    lWall += fw *0.01*(1. - pow(wp,40.));
    d.x = min(d.x, lWall);

    lWall = -p.x + wallW;
    
    
    lWall += fw *0.01*(1. - pow(wp,40.));
    lWall += wp*0.01;

    
    d.x = min(d.x,lWall );
    //d.x = min(d.x, -p.y - 0.4);
    
    p = pmod(p,2.);
    
    //d = dmin(d, vec2(length(p.yz) - 0.02,2.));
    
    //d.x = min(d.x, -p.y + wallW);
    d.x *= 0.5;
	return d;
}

vec3 glow = vec3(0);

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e6);
	p = ro;
    hit = false; t = 0.;
    
    for(int i = 0; i < 150 ; i++){
    	d = map(p);
        glow += exp(-d.x*100.);
        if(d.x < 0.001){
        	hit = true;
            break;
        }
        t += d.x;
        p = ro + rd*t;
    }
    
    return d;
}

vec3 getNormal(vec3 p,float sens){
	vec2 t = vec2(sens,0);
	return normalize(map(p).x - vec3(
    	map(p - t.xyy).x,
    	map(p - t.yxy).x,
    	map(p - t.yyx).x
    ));	
}
vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
    uv *= 0.8;
	vec3 dir = normalize(lookAt - ro);
	vec3 right = normalize(cross(vec3(0,1,0), dir));
	vec3 up = normalize(cross(dir, right));
    return normalize(dir + right*uv.x + up*uv.y);
}



vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. + dot(uv,uv)*0.3;
    vec3 col = vec3(0);
	
    vec3 ro = vec3(0);
    ro.y -= 0. + sin(TIME*0.10)*0.9;
    ro.x += 0.14;
    //ro.z += sin(mx)*4.;
    ro.z += TIME*0.2;
    vec3 lookAt = vec3(0,0,ro.z + 9.);
    
    vec3 rd = getRd(ro, lookAt, uv);
    rd.xy*= rot(-0.9 + sin(TIME*0.24)*0.5);
    rd.xz *= rot(0.5 + sin(TIME*0.3)*0.25);
    
    
    
    float t; bool hit; vec3 p;
    
    vec2 d = march(ro, rd, p, t, hit);
    
    if (hit){
        
        float wpsc = smoothstep(0.,1.,wp*1.);
        wpsc = smoothstep(0.,1.,wp*1. + 0.83);
        wpsc = pow(wpsc, 95.);
        
        vec3 albedo = mix(vec3(0.3,0.0 ,0.+ pow(fw,4.)*0.4), vec3(0.99,0.9,0.2), wpsc);
        
        
        vec3 lightCol = vec3(1.2,1.,1.);

        //vec3 L = normalize(l - p);
        vec3 L= normalize(vec3(0.01,0.2,0.7)); vec3 l = L;
        
        vec3 n = getNormal(p, clamp(0.001 + wp*0.008, 0., 1.));
        vec3 H = normalize(L - rd);


        int id = int(d.y);
        float METALNESS = 0.0;
        METALNESS += wp*0.9 + wpsc*0.2;
        float ROUGHNESS = 0.4 + clamp((1. - wpsc)*0.4*fw*2., 0., 0.6);
        
        
        vec3 F0 = vec3(0.03);
        vec3 N = n;
        vec3 V = normalize(ro - p);

        F0 = mix(F0, albedo, METALNESS);

        // calculate per-light radiance
        float distL    = length(l - p)*1.;
        float attenuation = 1. / (distL * distL);
        
        //attenuation = clamp(attenuation, 0., 1.);
        attenuation = 1.;
        vec3 radiance     = lightCol * attenuation;        

        // cook-torrance brdf
        float NDF = DistributionGGX(N, H, ROUGHNESS);   
        float G   = GeometrySmith(N, V, L, ROUGHNESS);      
        vec3 F    = fresnelSchlick(clamp(dot(N, V), 0.0, 1.0), F0);     

        vec3 kS = F;
        vec3 kD = vec3(1.0) - kS;
        kD *= 1.0 - METALNESS;	  

        vec3 numerator    = NDF * G * F;
        float denominator = 4.0 * max(dot(N, V), 0.0) * max(dot(N, L), 0.0);
        vec3 specular     = numerator / max(denominator, 0.001);  

        // add to outgoing radiance Lo
        float NdotL = max(dot(N, L), 0.0); 
        
        col += (kD * albedo / PI + specular) * radiance * NdotL * attenuation; 
        
    }
    
    col = mix(col, vec3(0.3,0.14,0.1)*2., smoothstep(0.,1.,t*0.01));
    
    
    uv.y *= 1.5;
    col *= 1. - (1. - vec3(0.1,0.6,0.1))*smoothstep(0.,1.,dot(uv,uv)*0.9)*0.4;

    col = max(col, 0.);
    col = pow(col, vec3(0.4545));
    
    fragColor = vec4(col,t);
	return fragColor; 
 } 




			//******** BuffB Code Begins ********



vec4 renderPassB() {
	vec4 C = vec4(0.0);
	vec2 fragCoord = _xy;

    C = vec4(0);
    
    C += gaussian(BuffA, fragCoord,R,MIPLVL);
    //C += T(fragCoord);
    
    
	return C; 
 } 


			//******** BuffC Code Begins ********



vec4 renderPassC() {
	vec4 C = vec4(0.0);
	vec2 fragCoord = _xy;

    C = vec4(0);
    
    C += gaussian(BuffB, fragCoord,R,MIPLVL);
    //C += T(fragCoord);
    
    
	return C; 
 } 


			//******** BuffD Code Begins ********



vec4 renderPassD() {
	vec4 C = vec4(0.0);
	vec2 fragCoord = _xy;

    C = vec4(0);
    
    C += gaussian(BuffC, fragCoord,R,MIPLVL);
    //C += T(fragCoord);
    
    C *= smoothstep(0.,1., pow(length(C.xyz)*3.4, 2.));
    
    
    
	return C; 
 } 


// a 2d apollonian-ish fractal used as displacement
// PBR from learnopengl.com

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragColor -= fragColor; // Wondering if this will glitch out for anyone. if it does, please leave a comment.
    // fragColor = vec4(0);
    
    float t = texture(BuffA, fragCoord/R).w;
    float st = 20.;
    
    vec2 uv = (fragCoord - 0.5*R)/R.y;
    float chromAbAmt = smoothstep(0.,1., dot(uv,uv))*0.5;
    // DOF and chromatic abberation
    for(float i = 0.; i < st; i++){
        float sz = 0.00;
        sz += smoothstep(0.,1.,abs(t)*0.02 - 0.01);
        
        //sz += 0.1;
        vec2 c = vec2(
        	sin(tau*i/st)*sz,
        	cos(tau*i/st)*sz
        );
        
        c *= 1. + texture(image30, (fragCoord + i )/R).x*sz*100.;
        
        
    	fragColor += vec4(
            texture(BuffA, (fragCoord + vec2(1)*chromAbAmt)/R + c).x,
            texture(BuffA, (fragCoord - vec2(1,0)*chromAbAmt)/R + c).y,
            texture(BuffA, (fragCoord - vec2(1)*chromAbAmt)/R + c).z,
            0.
        );
    }
    fragColor /= st;
    
    vec4 bloom = texture(BuffD, fragCoord/R, 0.);
    //fragColor = mix(fragColor, bloom, length(bloom.xyz));
    fragColor += bloom*0.5;
    
    fragColor.r *= 1.06;
    fragColor.b *= 1.04;
    fragColor.b *= 0.95;
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderPassB();
	}
	if(PASSINDEX == 2){
		return renderPassC();
	}
	if(PASSINDEX == 3){
		return renderPassD();
	}
	if(PASSINDEX == 4){
		return renderMainImage();
	}
}