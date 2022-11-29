vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


			//******** Common Code Begins ********

#define ID_GOLD 1.
#define ID_OTHER 2.
#define ID_THIRD 3.

#define stepSz 0.25
#define eps 0.001
#define maxIters 280.
    

#define ambientCol vec3(0.5,0.5,0.5)
#define lightCol vec3(0.9,0.8,0.8)*1.
#define lDir normalize(vec3(-4.5,4,1.))
#define silkCol vec3(0.3,0.22,0.03)*3.
#define sssSilkCol vec3(0.9,0.3,0.3)*0.


#define pi acos(-1.)



#define pmod(p,a) mod(p - 0.5*a,a) - 0.5*a

vec2 gmuv;

vec3 getRd( inout vec3 ro, vec3 lookAt, vec2 uv, float sc){
    vec3 dir = normalize(lookAt - ro);
    vec3 right = normalize(cross(vec3(0,1,0),dir)); 
    vec3 up = normalize(cross(dir, right));    
    //ro += right*uv.x*sc;
    //ro += up*uv.y*sc;
    return normalize(dir + right*uv.x + up*uv.y);
}


float eass(float p, float g) {
    float s = p*0.45;
    for(float i = 0.; i < g; i++){
    	s = smoothstep(0.,1.,s);
    }
    return s;
}

vec3 acesFilm(const vec3 x) {
    const float a = 2.51;
    const float b = 0.03;
    const float c = 2.43;
    const float d = 0.59;
    const float e = 0.14;
    return clamp((x * (a * x + b)) / (x * (c * x + d ) + e), 0.0, 1.0);
}


mat3 getOrthogonalBasis(vec3 direction){
    direction = normalize(direction);
    vec3 right = normalize(cross(vec3(0,1,0),direction));
    vec3 up = normalize(cross(direction, right));
    return mat3(right,up,direction);
}

float cyclicNoise(vec3 p, float time){
    float noise = 0.;
    
    // These are the variables. I renamed them from the original by nimitz
    // So they are more similar to the terms used be other types of noise
    float amp = 1.;
    const float gain = 0.75;
    const float lacunarity = 1.4;
    const int octaves = 4;
    
    
    
    
    const float warp = 0.75;    
    float warpTrk = 1.4 ;
    const float warpTrkGain = 1.25;
    
    // Step 1: Get a simple arbitrary rotation, defined by the direction.
    vec3 seed = vec3(-1,-2.,0.5);
    mat3 rotMatrix = getOrthogonalBasis(seed);
    
    for(int i = 0; i < octaves; i++){
    
        // Step 2: Do some domain warping, Similar to fbm. Optional.
        
        p += sin((p.zxy + vec3(-sin(time)*0.1,time*0.25,time*0.4))*warpTrk - 2.*warpTrk)*warp; 
    
        // Step 3: Calculate a noise value. 
        // This works in a way vaguely similar to Perlin/Simplex noise,
        // but instead of in a square/triangle lattice, it is done in a sine wave.
            
        //noise += sin(dot(cos(p), sin(p.zyx)))*amp;
        
        float f = sin(dot(cos(p), sin(p.zyx)));
        //f = sign(f)*pow(abs(f),1.);
        noise += f*amp;
        
        
        // Step 4: Rotate and scale. 
        
        p *= rotMatrix;
        p *= lacunarity;
        
        warpTrk *= warpTrkGain;
        amp *= gain;
    }
    
    return (noise*0.25 + 0.5);

    //return 1. - abs(noise)*0.5;
}


float cyclicNoiseFog(vec3 p, float time){
    float noise = 0.;
    
    // These are the variables. I renamed them from the original by nimitz
    // So they are more similar to the terms used be other types of noise
    float amp = 1.;
    const float gain = 0.55;
    const float lacunarity = 1.5;
    const int octaves = 4;
    
    
    
    
    const float warp = 0.5;    
    float warpTrk = 1.5 ;
    const float warpTrkGain = .25;
    
    // Step 1: Get a simple arbitrary rotation, defined by the direction.
    vec3 seed = vec3(-1,-2.,0.5);
    mat3 rotMatrix = getOrthogonalBasis(seed);
    
    for(int i = 0; i < octaves; i++){
    
        // Step 2: Do some domain warping, Similar to fbm. Optional.
        
        p += sin((p.zxy + vec3(-sin(time)*0.1,time*0.25,time*0.4))*warpTrk - 2.*warpTrk)*warp; 
        p.xy += MouseXY.xy;
        // Step 3: Calculate a noise value. 
        // This works in a way vaguely similar to Perlin/Simplex noise,
        // but instead of in a square/triangle lattice, it is done in a sine wave.
            
        //noise += sin(dot(cos(p), sin(p.zyx)))*amp;
        
        float f = sin(dot(cos(p), sin(p.zyx)));
        //f = sign(f)*pow(abs(f),1.);
        noise += f*amp;
        
        
        // Step 4: Rotate and scale. 
        
        p *= rotMatrix;
        p *= lacunarity;
        
        warpTrk *= warpTrkGain;
        amp *= gain;
    }
    
    return (noise*0.25 + 0.5);

    //return 1. - abs(noise)*0.5;
}


float opSmoothUnion( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); }

float opSmoothIntersection( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) + k*h*(1.0-h); }





#define pi acos(-1.)

#define tau (2.*pi)
#define pal(a,b,c,d,e) ((a) + (b)*sin((c)*(d) + (e)))

#define rot(x) mat2(cos(x),-sin(x),sin(x),cos(x))

vec3 glow = vec3(0);




float map(vec3 p, float t){
	float d = 10e7;
	
    vec4 q = vec4(p, 1.);
    
    float id = floor(p.z);
    
    
    
    for(float i = 0.; i <0. ; i++){
    	q.xyz = abs(q.xyz);
        q.xy *= rot(0.125*pi + i);
    
        q -= 0.1;
    }
    
    
    float dpp = dot(q.xyz,q.xyz);
    
    
    //    q = q/dpp;
    
    vec4 j;
    
    //q /= dot(q.xyz,q.xyz);
    vec3 b = vec3(2.4, 0.6, 1.);
    
    for(int i = 0; i <4 ; i++){
    	q.xyz = abs(q.xyz);
        q.xy *= rot(0.25*pi);
        q.xz *= rot(1.5);
    	q.x += 0.1;
    }
    
    for(float i = 0.; i < 10.; i++){
    	q.xyz = abs(mod(q.xyz - 0.5*b, b) )- 0.5*b;
        
        float dpp = dot(q.xyz,q.xyz);
        
        
        dpp = clamp(dpp, 0. ,1.54);
        if(i == 2.)
            j = q;
        
        q = q/dpp;
        if(i == 20.)
        	q.xz *= rot(0.7 );
        
    }
    
    q.xyz *= 1. ;
    
    float db = length(j.yx)/q.w - 0.01;
    
        
	float da = length(q.yz)/q.w - 0.02;
    
    float sc = 0.5;
    d = min(db,da);
    //d = da;
    d *= 0.5;
    d += smoothstep(1.,0.,t*.5)*0.7;
    d = abs(d) + 0.003;
    
    //d += exp(-t*4.)*0.7;
    
    
    vec3 c = vec3(1,1.,1.);
    da *= 0.5;
    db *= 0.5;
    da = abs(da) + 0.003;
    db = abs(db) + 0.003;
    glow += 0.9/(0.01 + da*da*1500.)*c;
    glow -= 0.9/(0.01 + db*db*2000.)*c*vec3(0.,0.8,2.);
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
    ro.z += TIME;
    
    float T = TIME*0.2;
    ro.xy += vec2(cos(T), sin(T))*3.;
    
    vec3 lookAt = vec3(0);
    
    lookAt.z = ro.z + 4.;
    
    vec3 rd = getRd(ro, lookAt, uv);
    
    float d;
    vec3 p = ro; float t = 0.; bool hit = false;
    
    for(int i = 0; i < 60; i++){
    	d = map(p, t);
        if(d < 0.001){
        	hit = true;
            //break;
        }
		t += d;
    	p = ro + rd*t;
    }
    
    
    col -= glow*0.001;
    
    col = pow(col, vec3(0.454545));
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}