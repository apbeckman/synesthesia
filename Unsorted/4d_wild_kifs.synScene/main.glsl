    

			//******** Common Code Begins ********

// Created by sebastien durand - 2021
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// ----------------------------
// based on [wyatt] kaleidoscope iterative function - https://www.shadertoy.com/view/MdKyRw
// ----------------------------


// Kifs paremeters
#define ITER 17
#define SCALE .766
#define ADD -.75*vec4(3.,-.075,0.,2.3)


// Ray Marching parameters
#define RAY_STEP 100
#define DIST_MAX 5.5

mat4 rot1, mRot;
vec4 closest = vec4(999.,0,0,0);

// 4D adaptation of Kif fractal
vec4 map(vec4 p) {
    vec4 ot = vec4(0);
    float t = .09;
    for (int i = 0; i < ITER; i++) {
        t = t*SCALE;
        p += t*ADD;      
        p.wxyz *= mRot;
        p = abs(p) - t;
        ot += p/t;
    }
    return vec4((length(p)-2.*t),         // Distance
                4.*ot.xwz/(pow(abs(ot.y),1.7)+.01)); // Color (4th Dimension give electric colors tuch !!!)
}

vec4 castRay(vec4 ro, vec4 rd) {
	const float precis = .005;
    float h = precis*2., t = 1.;
	vec4 res;
    closest = vec4(999.,0,0,0);
    for (int i=0; i<RAY_STEP; i++ ) {
        if (abs(h)<precis || t>DIST_MAX) break;
        t += h;
        res = map( ro+rd*t );
        h = res.x;
        if (h < closest.x) // get closest for halo
            closest.x = h;
        closest.yzw += res.yzw; // halo color
    }
    return vec4( t, res.yzw );
}

float softshadow(vec4 ro, vec4 rd, float mint) {
	float res = 1.,
          h,t = mint;
    for( int i=0; i<16; i++ ) {
        h = map( ro + rd*t ).x;
        res = min( res, 7.*h/t );
        t += .028;
    }
    return clamp( res-.6, .0, 1. );
}

vec4 calcNormal(vec4 p) {
    const vec2 e = vec2( 1e-3, 0.);
	return normalize(vec4(
	    map(p+e.xyyy).x - map(p-e.xyyy).x,
	    map(p+e.yxyy).x - map(p-e.yxyy).x,
	    map(p+e.yyxy).x - map(p-e.yyxy).x,
		map(p+e.yyyx).x - map(p-e.yyyx).x
	));
}

float calcAO(vec4 p, vec4 n ){
	float dd, hr, ao = 0., k = 1.;
    vec4 pos; 
    for( int aoi=0; aoi<5; aoi++ ) {
        hr = .01 + .05*float(aoi);
        pos =  n * hr + p;
        ao += -(map(pos).x-hr)*k;
        k *= .75;
    }
    return clamp( 1. - 4.*ao, 0., 1. );
}

vec3 render(vec4 ro, vec4 rd , vec3 backColor, out float d, vec4 lig){ 
    vec3 col;
    vec4 res = castRay(ro,rd);
    float t = res.x;
	vec3 uvw = .85*res.yzw;
    
    if (t<DIST_MAX) {
        vec4 pos = ro + t*rd,
             nor = calcNormal( pos );

		col = vec3(.4) + .6*abs(uvw);
		
        float ao = calcAO( pos, nor ); ao*=ao;
        float dif = clamp( dot( nor, lig ), 0., 1. ),
              bac = clamp( dot( nor, normalize(vec4(-lig.x,0.,-lig.z,0.))), 0., 1. )*clamp(1.-pos.y,0.,1.);

		float sh = 1.;
		if (dif>.02) { 
            sh = softshadow( pos, lig, .025); 
            dif *= sh; 
        }

		vec3 brdf =  .1*vec3(.10,.11,.13)*ao;
             brdf += .2*bac*vec3(.15)*ao;
             brdf += .8*dif*vec3(1,.9,.7);

		float pp = clamp( dot( reflect(rd,nor), lig ), 0., 1. ),
              spe = sh*pow(pp,16.),
              fre = ao*pow( clamp(1.+dot(nor,rd),0.,1.), 2. );

		col = col*brdf + 2.*(.5+.5*col)*spe + .4*fre*(.6+.4*col);
	
    } else {
        col = mix(backColor, clamp(.004*closest.yzw,0.,1.), smoothstep(.42,0.,pow(closest.x,.4)));
    }
    
    d = t;
	return vec3( clamp(col,0.,1.) );
}

// Rotation Matrix to apply to 4D objects
mat4 Rot4(float a, float b, float c) {        
    float c1 = cos(a), s1 = sin(a), 
          c2 = cos(b), s2 = sin(b), 
          c3 = cos(c), s3 = sin(c);	
    return mat4(c2,  s2*s3,   0, -s2*c3,   
                 0,  c1*c3, -s1,  c1*s3,
                 0,  c3*s1,  c1,  s1*s3,
                s2, -c2*s3,   0,  c2*c3);
}


void mainImage2(out vec4 fragColor, vec2 fragCoord, vec2 R, vec2 M, float TIME ) {

    vec2 q = fragCoord.xy/R.xy;
    vec2 p = -1.0+2.0*q;
    p.x *= R.x/R.y;
    
    // Noisy background
    float h = dot(vec3(q,1.),vec3(127.1,311.7,758.5453123));	
	vec3 colorSum = .75*(vec3(.0512) + .05*fract(sin(h)*43758.5453123));
    
    float d = 999.;
    if (length(p)<.92) {
        
        vec2 mo = M.xy/R.xy;	 
        float time = .25*smoothTimeC;
        // Rotations
        mRot = Rot4(.1*time, .351*time+2., .232*time+1.3);
        rot1 = Rot4((smoothTimeC*0.125-3.)/2.031, 1.+(smoothTimeC*0.125-3.)/2.1, .1*smoothTimeC*0.125);
        
        // Camera (real cam4D definition available at: https://www.shadertoy.com/view/4tX3Rn)
        vec4
            ro = vec4(3.2*cos(.124*smoothTimeC + 6.*mo.x+ 1.), .125 + 2.*mo.y, 3.2*sin(.124*smoothTimeC+ 6.*mo.x+1.),0),
            ta = vec4(0),
            cw = normalize( ta-ro ),
            cp = vec4(0,1,0,0),
            cu = normalize(vec4(cross(cw.xyz,cp.xyz),0)),
            cv = normalize(vec4(cross(cu.xyz,cw.xyz),0)),
            rd = normalize( p.x*cu + p.y*cv + 2.5*cw ),
            light = normalize(-cw*.5-cu+cv+.5*cp);
		
        // Rotation of 4D scene
        ro *= rot1;
    	rd *= rot1;
        light *= rot1;
        
        // Render
		colorSum = render( ro, rd, colorSum, d, light);
    }
    
    // Post process
    vec3 col = pow(colorSum.xyz,vec3(.56));
    col *= pow(16.*q.x*q.y*(1.-q.x)*(1.-q.y), .5);    
	fragColor = vec4(col, d);
}



			//******** BuffA Code Begins ********

#define tyme (smoothTime + cos(.125*smoothTime)-17.)

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    mainImage2(fragColor, fragCoord, RENDERSIZE.xy, _mouse.xy, tyme);  
	return fragColor; 
 } 


#define AA
//#define TIME (TIME + cos(.75*TIME)-17.)
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 U = _xy;
 
    O = texture(BuffA, U/RENDERSIZE.xy);

#ifdef AA
    if (O.w < DIST_MAX) {
        vec4 T;         
        for (float x=-.3; x<.4; x+=.6) {              
            for (float y=-.3; y<.4; y+=.6) {  
                mainImage2(T, U+vec2(x,y), RENDERSIZE.xy, _mouse.xy, tyme);  
                O += T;
            }
        }
        O /= 5.;
    }
#endif

	return O; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return renderMainImage();
	}
}