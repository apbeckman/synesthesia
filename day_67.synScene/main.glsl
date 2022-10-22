

			//******** Common Code Begins ********

#define MIPLVL 1.

#define PI acos(-1.)
#define tau (2.*PI)
#define R RENDERSIZE.xy
#define t(uv) texture(iChannel0, (uv)/R)

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


			//******** BuffA Code Begins ********


vec3 glow = vec3(0);
vec3 glowB = vec3(0);
#define pi acos(-1.)
    #define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))
float lpscale = 4.;

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

#define dmin(a,b) a.x < b.x ? a : b
vec2 map(vec3 p){
	vec2 d = vec2(10e7);
    
    vec3 q = p;
	float r = length(p);
	p = vec3(log(r), acos(p.z / length(p)), atan(p.y, p.x));

	// Get a scaling factor to compensate for pinching at the poles
	// (there's probably a better way of doing this)
	float xshrink = 1./(abs(p.y-pi)) + 1.0/(abs(p.y)) - 1.0/pi;

	// Scale to fit in the ]-pi,pi] interval
	p *= lpscale;
	// Apply rho-translation, which yields zooming
	p.x -= 0. + TIME;
	
    
	// Turn tiled coordinates into single-tile coordinates
	p = fract(p*0.5) * 2.0 - 1.0;
	p.x *= xshrink;
    
    
    // 
    p = abs(p);	
    p.xz *= rot(0.2*pi);
    float thickness = 0.1 ;
    

    
	// Get cylinder distance
	float ret = max(abs(p.x)- thickness, abs(p.z)- thickness) ;
    
	ret = min(ret, max(abs(p.x)- thickness, abs(p.y)- thickness) );
    
	d.x = ret;
	d.y = 0.;
    
    vec2 dB =vec2(
		ret = min(
            max(abs(p.x)- thickness*1.3, abs(p.z)- thickness*0.3), 
            max(abs(p.x)- thickness*1.3, abs(p.y)- thickness*0.3) )
        , 1.);
	d = dmin(d, dB);
    
    
    glowB += exp(-dB.x*100.);
	// Compensate for all the scaling that's been applied so far
	float mul = r/lpscale/xshrink;
	d.x = d.x * mul ;
    
    
    d.x *= 0.5;
	return d;
}
vec2 march(vec3 ro, vec3 rd, inout vec3 p, inout float t, inout bool hit){
	vec2 d = vec2(10e7);
	hit = false; t = 0.; p = ro;
    
    for(int i = 0; i < 160; i++){
    	d = map(p);
        glow += exp(-d.x*900.);
        if(d.x < 0.00001){
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
	vec3 up = normalize(cross( dir, right));
	return normalize(dir + right*uv.x + up*uv.y);
}

vec3 getNormal(vec3 p){
	vec2 t = vec2(0.00001,0);
    return normalize(
    	vec3(
        	map(p+t.xyy).x - map(p-t.xyy).x,
        	map(p+t.yxy).x - map(p-t.yxy).x,
        	map(p+t.yyx).x - map(p-t.yyx).x
        )
    );
}

vec3 getNormalO(vec3 p){
	vec2 t = vec2(0.00005,0);
	return normalize(map(p).x - vec3(
    	map(p-t.xyy).x,
    	map(p-t.yxy).x,
    	map(p-t.yyx).x
    ));
}



float calcSoftshadow( in vec3 ro, in vec3 rd, in float mint, in float tmax, int technique )
{
	float res = 1.0;
    float t = mint;
    float ph = 1e10; // big, such that y = 0 on the first iteration
    for( int i=0; i<52; i++ )
    {
		float h = map( ro + rd*t ).x;

        // traditional technique
        if( technique==0 )
        {
        	res = min( res, 10.0*h/t );
        }
        // improved technique
        else
        {
            // use this if you are getting artifact on the first iteration, or unroll the
            // first iteration out of the loop
            //float y = (i==0) ? 0.0 : h*h/(2.0*ph); 

            float y = h*h/(2.0*ph);
            float d = sqrt(h*h-y*y);
            res = min( res, 10.0*d/max(0.0,t-y) );
            ph = h;
        }
        
        t += h;
        
        if( res<0.01 || t>tmax ) break;
        
    }
    return clamp( res, 0., 1.0 );
}

#define T (TIME*0.25)
#define mx (10.*_mouse.x/RENDERSIZE.x)
#define my (10.*_mouse.y/RENDERSIZE.x)
vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    uv *= 1. + dot(uv,uv)*(0.56); // squish coords
    
    vec3 col = vec3(0);

    vec3 ro = vec3(0);
    
    ro.z -= 1.5;
    ro.x -= sin(2.*T/pi)*(0.9);
    ro.y -= cos(2.*T/pi)*(0.9);
    
    vec3 lookAt = vec3(0);
    vec3 rd = getRd(ro, lookAt, uv);
    rd.xz *= rot(sin(TIME*0.5)*0.1);
    rd.xy *= rot(sin(TIME*0.5)*0.1);
    //rd.xz *= rot(0.5);
    vec3 p; float t; bool hit;
    vec2 d = march(ro, rd, p, t,hit);
    #define ao(d) clamp(map(p+n*d).x/d,0.,1.)

    vec3 lightColA = vec3(0.1,0.4,0.9);
    vec3 lightColB = vec3(0.25,0.2,0.20)*0.1;
    if(hit){
        // do lighting from two light sources
        vec3 l = normalize(-p);
        
        l.xz *= rot(sin(TIME*0.12 + 0.2)*0.1 + 0.3);
        l.xy *= rot(sin(TIME*0.08 + 0.4)*0.1);
        
        float shad = calcSoftshadow( p,l, 0.1, 100., 1 );
        
        shad = pow(smoothstep(0.,1.,shad*20. ), 2.);

    	vec3 n = getNormal(p);
        vec3 h = normalize(l - rd);
        
        float diff = max(dot(n, l),0.);
        float spec = pow(max(dot(n, h),0.), 20.);
        float fres = pow(max(1. - dot(n, -rd),0.), 5.);
		
        vec3 lightCol = lightColA;
        vec3 albedo = vec3(0.5,0.5,0.7)*1.;
        if(d.y == 1.){
        	albedo = vec3(1.5,1.,1.)*1.;
        }
        
        col += albedo*diff*lightCol + spec*lightCol + fres*lightCol;
        
        col += albedo*0.02*lightCol;
        
        //l = normalize(vec3(-0.4,0.1,-1));
        l = normalize(p);
        h = normalize(l - rd);
        
        diff = max(dot(n, l),0.);
        spec = pow(max(dot(n, h),0.), 20.);
        fres = pow(max(1. - dot(n, -rd),0.), 5.);
		
        lightCol = lightColB;
        
        col += (albedo*diff*lightCol + spec*lightCol + fres*lightCol)*0.9;
        
        col += albedo*0.01*lightCol;
        
        col = col*shad;
    }
    col = mix(col, col + glow*0.03*vec3(0.1,0.4,0.9),
              pow( smoothstep(0.,1.,exp(-length(p)*1.5)), 4.)
             ); 
    col = mix(col,lightColB, smoothstep(0.,1.,t*0.1 ));
    
    col = clamp(col,0., 20.);
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
    
    C *= smoothstep(0.,1., pow(length(C.xyz)*1., 2.));
    
    
    
	return C; 
 } 




// LOG SPHERICAL MAPPING FROM THIS AWESOME ARTICLE
// https://www.osar.fr/notes/logspherical/
// ^^ $$ ^^ $$ ^^

// it's super easy to use to as a tool
// and very interesting to read about

// code and comments in map function also from the author

// shadows from inigo quilez

// color scheme stolen from evvvvil's https://www.shadertoy.com/view/wt33RN

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragColor -= fragColor; 
    // fragColor = vec4(0);
    
    float t = texture(BuffA, fragCoord/R).w;
    float st = 40.;
    
    vec2 uv = (fragCoord - 0.5*R)/R.y;
    float chromAbAmt = smoothstep(0.,1., dot(uv,uv))*1.4;
    
    fragColor += vec4(
        texture(BuffA, (fragCoord + vec2(1)*chromAbAmt)/R).x,
        texture(BuffA, (fragCoord - vec2(1,0)*chromAbAmt)/R ).y,
        texture(BuffA, (fragCoord - vec2(1)*chromAbAmt)/R ).z,
        0.
    );
    
    
    
    vec4 bloom = texture(BuffD, fragCoord/R, 0.);
    //fragColor = mix(fragColor, bloom, length(bloom.xyz));
    fragColor += bloom*0.6;
    
    
    fragColor *= 1. - dot(uv,uv);
    fragColor.r *= 1.06;
    fragColor.b *= 1.04;
    fragColor.b *= 0.95;
    fragColor.xyz *= 1.9;
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