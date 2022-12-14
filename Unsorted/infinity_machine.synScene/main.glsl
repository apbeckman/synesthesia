//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define MAX_STEPS 	100
#define MAX_DIST	150.
#define MIN_DIST	.001
#define EPSILON		.0001

#define PI 3.1415926535897

// Change to 2 to enable antialiasing
#define AA 1

#define ZERO (min(FRAMECOUNT,0))

/**
	Moar transparent marching
	Not sure I've got the motion thing
	right - using some mods but it works!
	100% sure there is some odd code!

	New Camera Rig - way better!
	VR Render updated and it kicks butt!!
*/

vec2 hash2(vec2 p) {
	vec2 o = (p+0.5)/256.;
  	return o;
}

float goldNoise(vec2 coord, float seed) {
	float phi = 1.6180339887498 * 0000.1;
  	float pi2 = PI * 0.1;
  	float sq2 = 1.4142135623730 * 10000.;
  	float temp = fract(
    	sin(dot(coord*(seed+phi), vec2(phi, pi2))) * sq2
  	);
  	return temp;
}

float sdBox(vec3 p,vec3 s) {
    p = abs(p) - s;
    return max(max(p.x,p.y),p.z);
}

float sdSphere(vec3 p, float r) {
    return length(p)-r;    
}

// Single rotation function - return matrix
mat2 r2(float a){ 
  return mat2(cos(a), sin(a), -sin(a), cos(a)); 
}

// iMouse pos function - take in a vec3 like ro
// simple pan and tilt and return that vec3
/*vec3 get_mouse(vec3 ro) {
    float x = iMouse.xy==vec2(0) ? -.2 :
    	(iMouse.y / RENDERSIZE.y * .5 - 0.25) * PI;
    float y = iMouse.xy==vec2(0) ? .0 :
    	-(iMouse.x / RENDERSIZE.x * 1.0 - .5) * PI;
    float z = 0.0;

    ro.zy *= r2(x);
    ro.zx *= r2(y);
    
    return ro;
}
*/
float map(vec3 pos) {
	// set up size for repetition
   	float size = 20.;
	float rep_half = size/2.;
    // get center vec and some movement
    vec3 center = vec3(0.,0., 0.- TIME * 5.);
    vec3 pp = pos-center;
	// make the id's for beams
    vec3 pi = vec3(
        floor((pp + rep_half)/size)
    );
	// make vec and ids for balls
    vec3 tt = pos-center;
    vec3 ti = vec3(
        floor((tt + rep_half)/size)
    );

    // SHAKKADOO!! movement change based on time mod
    // hash/noise vector direction.
    // Using a mod of time < segment flip direction
    if(mod(script_time+script_bass_time*0.25,6.)<2.){
        //up-down
        // fhs gets the hash of the other vec2 
        // its not moving in.
        float fhs = goldNoise(pi.xz, 312.);
        if(fhs>.5){
            tt.y += fract(script_time+script_bass_time*0.25 ) * 20.;
        } else {
            tt.y -= fract(script_time+script_bass_time*0.25) * 20.;
        } 
    } else if (mod(script_time,6.)<4.){
        //left-right
      	float fhs = goldNoise(pi.zy, 312.);
        if(fhs>.5){
            tt.x += fract(script_time+script_bass_time*0.25) * 20.;
        } else {
            tt.x -= fract(script_time+script_bass_time*0.25 ) * 20.;
        }   
    } else {
        //back-forth
      	float fhs = goldNoise(pi.xy, 312.);
        if(fhs>.5){
            tt.z += fract(script_time+script_bass_time*0.25 ) * 20.;
        } else {
            tt.z -= fract(script_time+script_bass_time*0.25) * 20.;
        }   
    }

    pp =  mod(pp+rep_half,size) - rep_half;
	tt =  mod(tt+rep_half,size) - rep_half;
    
    float len = 11.;
    float tx =  1.5;

    //framework
    float d1 = sdBox(abs(pp)-vec3(10.,10.,0.), vec3(tx,tx,len) );
    float d2 = sdBox(abs(pp)-vec3(10.,0.,10.), vec3(tx,len,tx) );
    float d3 = sdBox(abs(pp)-vec3(0.,10.,10.), vec3(len,tx,tx) );
    d1 = min(min(d1,d3),d2);

    //box clips
	d2 = sdBox(abs(pp)-vec3(10.,10.,10.), vec3(3.5) );
    d3 = sdBox(abs(pp)-vec3(10.,10.,10.), vec3(2.5,4.5,2.5) );
    d1 = min(d1,min(d2,d3));
    
    d2 = sdBox(abs(pp)-vec3(10.,10.,10.), vec3(2.5,2.5,4.5) );
    d3 = sdBox(abs(pp)-vec3(10.,10.,10.), vec3(4.5,2.5,2.5) );
    d1 = min(min(d1,d2),d3);
    
    //balls
    float d4 = sdSphere(tt-vec3(0.,0.,0.),2.5);
    if(d4<d1) d1 = d4;
    
    return d1;
}

/**
	transparent marcher
	just go to the maxstep
	like ao plus dist
	fun evil render below
*/
float trans_ray( in vec3 ro, in vec3 rd ) {
    float depth = 0.;
    for (int i = 0; i<MAX_STEPS; i++) {
        vec3 p = ro + depth * rd;
        float dist = map(p);
        dist = max(abs(dist), MIN_DIST);
        // make it look all dark and evil here
        // dist = smoothstep(0.1,.6,max(abs(dist), MIN_DIST));
        depth += abs(dist*.5);
    }
    return depth;
}
    
mat3 get_camera(vec3 ro, vec3 ta, float rotation) {
	vec3 cw = normalize(ta-ro);
	vec3 cp = vec3(sin(rotation), cos(rotation),0.0);
	vec3 cu = normalize( cross(cw,cp) );
	vec3 cv = normalize( cross(cu,cw) );
	return mat3( cu, cv, cw );
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec3 color = vec3(0.);
    vec2 uv;
    // 
    #if AA>1
    for( int m=ZERO; m<AA; m++ )
    for( int n=ZERO; n<AA; n++ )
    {
    	// pixel coordinates
    	vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
		uv = (2. * fragCoord.xy - (RENDERSIZE.xy+o))/RENDERSIZE.y;
	#else    
        uv = (2. * fragCoord.xy - RENDERSIZE.xy )/RENDERSIZE.y;

	#endif
    
	float ftime = TIME * 5.5;
    float ftm = .5 + .5 * sin(TIME*.1);
    
    vec3 ta = vec3(0.0,0.0,0.);
    
    vec3 ro = vec3(0, 0, (script_bass_time*0.125+script_time*0.25)+rate);
    
    mat3 cameraMatrix = get_camera(ro, ta, 0. );
    vec3 rd = cameraMatrix * normalize( vec3(uv.xy, .85) );
    rd.yz = _rotate(rd.yz, lookXY.y*PI);
    rd.xy = _rotate(rd.xy, lookXY.x*PI);
    rd.xz = _rotate(rd.xz, lookZ*PI);
    float dist = trans_ray(ro, rd);
    vec3 col = vec3(0.);
    if(dist>MIN_DIST) {
        vec3 shade = vec3(dist*.005); 
        // fun with .1< and 1.> to filter
        // out the fields and get dark and
        // or washed out effects
        col = 1. * smoothstep(.1,1.,shade); 
        col = pow(col, vec3(0.4545));
        color += col;
    }

    // Color version if thats your thing..
   	 color *= vec3(.4-rd.y,	.8, .8);
    // color *= vec3(.6,	.4, .8-rd.y);
        
    // AA from NuSan
    #if AA>1
    	}
    	color /= float(AA*AA);
	#endif
    
    fragColor = vec4(color,1.0);
    return fragColor;
}

void mainVR( out vec4 fragColor, in vec2 fragCoord, in vec3 fragRayOri, in vec3 fragRayDir ) {
    vec2 uv = (2.*fragCoord-RENDERSIZE.xy)/RENDERSIZE.y;
    vec3 rd = vec3(0.0,0.0,0.);
    vec3 ro = vec3(4.0,4.,-7.);
	vec3 color = vec3(0.);
    float dist = trans_ray(fragRayOri + ro, fragRayDir + rd);
    vec3 col = vec3(0.);
    if(dist>MIN_DIST) {
        vec3 shade = vec3(dist*.005); 
        col = 1. * smoothstep(.1,1.,shade); 
        col = pow(col, vec3(0.4545));
        color += col;
    }

   	color *= vec3(.4-fragRayDir.y,	.8, .8);
  
	fragColor += vec4(color,1.0);
	//return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}