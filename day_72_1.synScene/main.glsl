

			//******** Common Code Begins ********


#define TIME (TIME + 2.5)
    #define T TIME
    
    #define SPEED 2.

			//******** BuffA Code Begins ********


#define dmin(a,b) a.x < b.x ? a : b

#define FL_H 0.4
#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))

#define pal(a,b,c,d,e) ((a) + (b)*sin(6.28*((c)*(d) + (e))))

vec3 glow = vec3(0);
vec3 glowB = vec3(0);
vec3 glowC = vec3(0);

vec3 reflAtten = vec3(1);
float randomO(vec2 u){
	return fract(sin(u.y*125.1 + u.x *125.625)*225.5235);
} 

float random(vec2 u){
	return texture(image30, (u/256.)).x;
} 

float noise(vec2 p) {
	vec2 i = ceil(p);
    vec2 f = smoothstep(0.5,1.,fract(p));
    f = smoothstep(0.7,1.,f);
    //f = smoothstep(0.2,1.,f);
    vec2 u = f * f * (3. - 2. * f);
   	float a = random(i);
    float b = random(i + vec2(1., 0.));
    float c = random(i + vec2(0., 1.));
    float d = random(i + vec2(1., 1.));
    return mix(mix(a, b, u.x), mix(c, d, u.x), u.y);
}

vec2 valueNoise(float p){
	vec2 a = texture(image30, vec2(floor(p))/256.).xy;
	vec2 b = texture(image30, vec2(floor(p) + 1.)/256.).xy;
    return mix(a,b,smoothstep(0.,1.,fract(p)));
}

float fbm(vec2 p) { 
    p *= 0.6;
	float s = .0;
	float m = .0;
	float a = .5;	
	for(int i = 0; i < 5; i++) {
		s += a * noise(p);
		m += a;
		a *= .7;
		p *= 2.;
	}
	return s / m;
}

vec3 path (float z){
    z *= 0.5;
	return vec3(
    	sin(z + cos(z*0.7))*0.7,
    	cos(z + cos(z*1.2))*0.6,
        0.
    )*1.;
}

#define pmod(p,x) mod(p,x) - x*0.5
vec2 map(vec3 p){
	vec2 d = vec2(10e6);
	p -= path(p.z);
    
    float n = fbm(p.xz)*1.;
    float m = pow(abs(sin(p.z*0.03)),10.);
    n *= 1. + m*2.;
    n *= 0.5;
    //n = pow(n*1., 4.)*3.;
    
    p.xy *= rot(sin(p.z*0.9 + sin( p.x*2. + p.z*4. + TIME*0.1)*0.9 + p.z*0.1)*0.6);
    
    
    float flTop =(-p.y + FL_H + n*0.7)*0.3;
    d = dmin(d, vec2(flTop,1.));
    
    
    float flBot =(p.y + FL_H + n*0.7)*0.3;
    
    d = dmin(d, vec2(flBot,1.));
    
    //d = dmin(d, vec2(length(p) - 0.1, 2.));
    
    vec3 z = p;
    p.xz = pmod(p.xz, 0.6);
    
    p.y = abs(p.y);
    p.y -= FL_H*0.7 + n*0.7;
    float dBalls = length(p);
    vec3 q = abs(p) - 0.04;
    
    
    float dPipes = max(q.x,q.y);
    float dPipesB = max(q.y,q.z);
    //d = dmin(d, vec2(dPipes, 2.));
    
    float atten = pow(abs(sin(z.z*0.82 + TIME*0.2)), 10.);
    float attenB = pow(abs(sin(z.z*0.92  + sin(z.x + TIME)*0.2 + sin(z.y*3.)*1. + TIME*0.5)), 100.);
    float attenC = pow(abs(sin(z.z*0.71  + sin(z.x + TIME)*0.2 + sin(z.y*3.)*4. + TIME*0.2)), 200.);
    
    vec3 col = pal(0.2,0.3 - attenC,vec3(0.1 + pow(abs(sin(TIME*1.)), 40. )*0.005,2.2,0.3),0.5 + sin(TIME)*0.005,0.5 - attenB*0.6)*pow(0.75+ syn_HighLevel , 2.);
    //vec3 col = pal(0.4,0.6,vec3(0.1 + pow(abs(sin(TIME*1.)), 40.)*0.0,2.2,0.3),0.5 + sin(TIME)*0.01,0.5 );
    
    //vec3 col = pal(0.4,0.6,vec3(0.1 + pow(abs(sin(TIME*1.)), 40.)*0.1,2.2,0.3),0.5 + sin(TIME)*0.01,0.5 - attenB*0.6);

    
    float sc = 60. - atten*55.;
    glowB += exp(-dPipes*sc)*col*reflAtten;
    glowB += exp(-dPipesB*sc)*col*reflAtten;
    //glowC += exp(-dBalls*90.)*colB;
    //glowB -= 0.002/(0.02 + dPipes*dPipes)*0.4;
    d.x *= 0.6;
    return d;
}

vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e6);
	p = ro; t = 0.; hit = false;
//    for (int i = 0; i < 150 ; i++){
    for (int i = 0; i < 125 ; i++){

    	d = map(p);
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

vec3 getRd(vec3 ro, vec3 lookAt, vec2 uv){
	vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0), dir));
    vec3 up = normalize(cross(dir, right));
	return normalize(dir + right*uv.x + up*uv.y);
}


vec3 getNormal(vec3 p){
	vec2 t = vec2(0.01,0);
    return normalize(
    	vec3(
        	map(p - t.xyy).x - map(p + t.xyy).x,
        	map(p - t.yxy).x - map(p + t.yxy).x,
        	map(p - t.yyx).x - map(p + t.yyx).x
        )
    );
}


#define mx (10.*_mouse.x/RENDERSIZE.x)
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    float m = pow(abs(sin(T/2.)), 5.);
    //uv *= 1. - dot(uv,uv)*(1. - pow(m,2.)*1.)*0.4;
    
    
    
    //uv.xy *= rot(0.1)
    vec3 col = vec3(0);
    
    vec3 ro = vec3(0);
    
    ro.z += mx*2.;
    ro.xy += valueNoise(TIME*40.)*(0.01)*m; // camshake
    
    ro.z += TIME*SPEED - sin(T)*SPEED;
    
    ro += path(ro.z);
    
    vec3 lookAt = vec3(0,0,ro.z + 1.);
    
    lookAt += path(lookAt.z);
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    rd.xy *= rot(sin(TIME)*0.1);
    
    
    
    //ro += rd*texture(image30, (uv*200. + TIME*9.)).x*2.;
    
    bool hit; float t; vec3 p;
    
    float bounce;
    
    float firstT = 0.;
    vec2 d;
    for(int i = 0; i < 2     ; i++){
        d = march(ro, rd, p, t, hit);
        vec3 n = getNormal(p);
        
        if(i == 0){
        	firstT = t;
        }
        reflAtten *= 0.53;
           
        rd = reflect(rd, n);
        ro = p + rd*0.1;
    }
    
	
    
    glowB = max(glowB, 0.);
    glowB = pow(glowB, vec3(1./0.45));
    col += glowB*0.0004;
    
    col += glowC*0.004;
    col += glow*0.1;
    col = mix(col, vec3(0.55,0.25,0.2)*0.01, pow(smoothstep(0.,1., firstT*0.08), 2.) );
    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 




vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 uv = fragCoord/RENDERSIZE.xy;
	vec2 uvn = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.xy;
    
    
    //float m = pow(abs(sin(p.z*0.03)),10.);

    // Radial blur
    float steps = 20.;
    float scale = 0.00 + pow(length(uv - 0.5)*1.2,2.7)*0.2;
    //float chromAb = smoothstep(0.,1.,pow(length(uv - 0.5), 0.3))*1.1;
    float chromAb = pow(length(uv - 0.5),1.4)*5.1;
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
    
    
    fragColor = radial*3.; 
    //fragColor = mix(fragColor,smoothstep(0.,1.,fragColor), 0.8);
    //1fragColor *= 18.;
    fragColor = max(fragColor, 0.);
    fragColor = pow(fragColor, vec4(0.4545 + dot(uvn,uvn)*1.7));
    fragColor *= 1. - dot(uvn,uvn)*1.2;
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