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
/*
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.75;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.95;
*/

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
        
        p += sin((p.zxy + vec3(-sin(smoothTime)*0.1,smoothTime*0.25,smoothTime*0.4))*warpTrk - 2.*warpTrk)*warp; 
    
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




			//******** BuffA Code Begins ********


#define rot(a) mat2(cos(a),-sin(a),sin(a),cos(a))
vec2 dmin(vec2 d, float db, float dbid){return d.x < db ? d : vec2(db,dbid);}
vec2 dmin(vec2 d, vec2 b){return d.x < b.x ? d : b;}

float wallNoise = 0.;
float wallNoiseB = 0.;

vec2 sdFrac(vec3 po){

    vec4 p = vec4(po,1.);
    
    vec2 d = vec2(10e5);
    
    p.xy -= 1.;
    
    float id = floor(p.z/4. + 0.5);
    p.xyz = pmod(p.xyz,vec3(2.,2.,4.));
    for(int i = 0; i < 4; i++){
        
        p.xyz = abs(p.xyz);
        
        if(p.x > p.y) p.xy = p.yx;
        if(p.x < p.z) p.xz = p.zx;
        if(p.y > p.z) p.yz = p.zy;
        
        p.xyz -= vec3(0.,0.1,0.3 - sin(id*2. )*0.1);
        
        p *= 4.5;
        
        //p.xy *= rot(0.125*pi);
        
    }
    
    
    p.xyz /= p.w;
    
    p.w = 1.;
    
    for(int i = 0; i < 4; i++){
        
        float dpp = dot(p.xyz,p.xyz*1.);
        p.xyz = pmod(p.xyz,vec3(4.3,2.,3.6));
        
        p.xyz = abs(p.xyz);
        
        p /= dpp;
        
        p.x -= 0.06;
        p.y -= 0.03;
        
        p.z -= 0.05;
        
        p.zx *= rot(-0.5*pi);
        p *= 3.;
        
        
    }
    
    p.xyz /= p.w;

    vec3 pp = p.xyz;
    float da = length(p.xy) - 0.03;
    //d = dmin(d, max(abs(p.x),abs(p.y)) - 0.03, ID_GOLD);
    
    p.xz *= rot( 0.*pi);
    p.x += 0.07;
    float db = max(abs(p.x),abs(p.z)) - 0.04;
    db = length(p.xz) - 0.03;
    
    
    pp.xz *= rot(-0.125*pi);
    
    pp -= vec3(-0.,0.01,0.1);
    
    
    pp.yz *= rot(0.25*pi);
    
    float dc = length(pp.yz) - 0.04;
    
    
    d = dmin(d, da, ID_GOLD);
    d = dmin(d,db, ID_OTHER);
    d = dmin(d,dc, ID_THIRD);
    //d = dmin(d, max(abs(p.x),abs(p.z)) - 0.02, ID_OTHER);
    
    //d = min(d,length(p) - 0.2);

    return d;
}


vec2 map(vec3 p){
    vec2 d = vec2(10e5);
    //vec2 dwall = sdWall(p); 
    vec2 dsilk = sdFrac(p);

    d = dmin(d, dsilk);
    
    return d;
}

float softshadow( in vec3 ro, in vec3 rd, float mint, float maxt, float k )
{
    float res = 1.0;
    float ph = 1e20;
    for( float t=mint; t<maxt; )
    {
        float h = map(ro + rd*t).x;
        if( h<0.001 )
            return 0.0;
        float y = h*h/(2.0*ph);
        float d = sqrt(h*h-y*y);
        res = min( res, k*d/max(0.0,t-y) );
        ph = h;
        t += h*0.4;
    }
    return res;
}

vec3 getNormal(vec3 p, float precis){
      vec3 n = vec3(0.0);
    for( int i=0; i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*map(p+e*precis).x;
    }
    return normalize(n);
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
 
    
    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    vec2 muv = (iMouse.xy)/RENDERSIZE.xy;
    
    gmuv = (iMouse.xy - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;
    
    //uv *= 1. + dot(uv,uv)*5.;
    
    muv = muv*2. - 1.;
    if (muv.x < -0.85){
        muv -= muv;
    }
    
    vec3 ro = vec3(0);
    
    ro.z -= 1. - smoothTime*0.2;    
    ro *= 1.9;
    
    
    
    vec3 lookAt = vec3(0);
    lookAt.z = ro.z + 2.;
    
    //lookAt.z += muv.x*0.3;
    lookAt.y += muv.y;
    
    
    lookAt += normalize(vec3(cos(TIME*0.5)*0.3,sin(TIME*0.7)/3,cos(TIME*0.3)/3.5))*0.3;
    
    vec3 rd = getRd(ro, lookAt, uv, 2.);
    
    
    
    //ro -= rd*1.;
    
    vec3 col = vec3(0);
    
    
    vec3 p = ro;
    float t = 0.; bool hit = false;
    vec2 d;
    
    vec3 atten = vec3(1.);
    
    
    vec3 hitCol = vec3(0);
    
    
    
    float i = 0.;
    vec3 reflalbedo = mix(ambientCol, lightCol,0.5);
    for(; i < maxIters; i++){

        p = ro + rd*(t += d.x*stepSz);
        d = map(p);

        if(d.x < eps){
            hit = true;
            break;
        }

    }

    
    vec3 n = getNormal(p,eps);
    n = normalize(n - lDir*cyclicNoise(p*220.,0.)*0.3);
        
                    
    {
        #define AO(n,a) clamp(map(p + normalize(mix(n,lDir,0.23))*a).x/a, 0., 1.)
        
        float ao = AO(n,0.3)*AO(n,0.5)*AO(n,0.04)*AO(n,.14)*1.;
         
        //float SSS = SSS(vec3(-1.,.5,0.9),1.9)*1.*SSS(vec3(-4.,2.,0.2),5.)*1.;//*AO(lDir,0.5)*AO(-n,0.9);

        vec3 halfV = normalize(lDir - rd);
        float fres = pow(1.-max(dot(-rd, n),0.0001),3.);
        float spec = pow(max(dot(n,halfV),0.),8.);
        float diff = dot(n,lDir);
                     
        vec3 albedo = vec3(0.);
        vec3 sssalbedo = vec3(0.);
        
        float shad = softshadow( p, lDir, 0.01, 10., 2.);
        
        //float shad = 1.;
        shad = min(shad,diff);
        
        vec3 colSilk = vec3(0);
        {
            map(p);
            albedo = silkCol;
            
            albedo -= albedo;
            
            vec3 r = reflect(rd,n);
            
            float fact = pow(length(sin(n*2.)*0.5 + 0.5),4.)/sqrt(3.); 
            float factb = pow(length(sin(r*3. + n*5. + 4. + smoothTime*0.)*0.5 + 0.5)/sqrt(3.),0.4); 
            float factc = pow(length(sin(r*5. + n*4. + 4. + smoothTime*2.2)*0.5 + 0.5),1.)/sqrt(3.); 
            
            factb = clamp(factb,0.,1.);
            factc = clamp(factc,0.,1.);
            albedo += mix( silkCol, vec3(0) + ambientCol*0.04 + silkCol*silkCol*silkCol*0.45,factb);;
            albedo = mix( albedo, vec3(0) + ambientCol*0.02 + silkCol*silkCol*0.2,pow(factc,1.)*1.);;
            
            
            
            albedo += lightCol*45.*(1.-pow(factb,.005))*(1. + silkCol*.3);
            colSilk = albedo*lightCol;

            //colSilk = mix(colSilk, reflalbedo, clamp(fres + spec,0.,1.)*ao*0.3);
            ao = max(ao,0.);
            shad = max(shad,0.);
            
            
            colSilk = mix((colSilk + colSilk*ambientCol)*0., colSilk, ao);
            //colSilk = mix((colSilk + colSilk*ambientCol)*0.5,colSilk, shad);
            
            
        }
        
        vec3 colOther = vec3(0);
        
        
        {
            map(p);
            
            albedo -= albedo;
            vec3 r = reflect(rd,n);
            
            float fact = pow(length(sin(r*4. - n*15.)*0.5 + 0.5)/sqrt(2.),0.2); 
            float factb = pow(length(sin(r*4. + n*6. + 0. + smoothTime*0.)*0.5 + 0.5),1.5)/sqrt(3.); 
            float factc = pow(length(sin(r*1. + n*5. + 4. + smoothTime*0.)*0.5 + 0.5),1.)/sqrt(3.); 
            
            factb = clamp(factb,0.,1.);
            factc = clamp(factc,0.,1.);
            
            
            vec3 otherCol = vec3(0.7,0.5,0.15)*0.1;
            albedo += mix( otherCol, vec3(0.3,0.6,0.4)*0.00,fact);;
            //albedo = mix( albedo, vec3(0) + ambientCol*0. + otherCol*otherCol*0.2,pow(factc,1.)*1.);;
            
            
            
            albedo += (1.-pow(factb,0.4))*vec3(1.8,1.2,0.55)*0.05;
            
            colOther = albedo;

            ao = max(ao,0.);
            shad = max(shad,0.);
            
            
            colOther = mix((colOther + colOther*ambientCol)*0.5, colOther, ao);
            colOther = mix((colOther + colOther*ambientCol)*0.5,colOther, shad);
            
            
        }
        
        
        vec3 colThird = vec3(0);
        
        {
            map(p);
            
            albedo -= albedo;
            vec3 r = reflect(rd,n);
            
            float fact = pow(length(sin(r*11. - n*2.)*0.5+ 0.5)/sqrt(3.),0.07); 
            float factb = pow(length(sin(r*4. + n*6. + 0. + smoothTime*0.)*0.5 + 0.5),1.5)/sqrt(3.); 
            float factc = pow(length(sin(r*1. + n*5. + 4. + smoothTime*0.)*0.5 + 0.5),1.)/sqrt(3.); 
            
            factb = clamp(factb,0.,1.);
            factc = clamp(factc,0.,1.);
            
            
            vec3 otherCol = vec3(0.5,0.5,0.4)*0.76;
            albedo += mix( otherCol, vec3(0.5,0.3,0.)*0.003,fact);;
            //albedo = mix( albedo, vec3(0) + ambientCol*0. + otherCol*otherCol*0.2,pow(factc,1.)*1.);;
            
            
            
            albedo += (1.-pow(factb,0.05))*vec3(1.8,1.2,0.55)*0.03;
            
            colThird = albedo;

            ao = max(ao,0.);
            shad = max(shad,0.);
            
            
            colThird = mix((colThird + colThird*ambientCol)*0.3, colThird, ao);
            colThird = mix((colThird + colThird*ambientCol)*0.4,colThird, shad);
            
            
        }
        
        hitCol += colSilk*float(floor(d.y) == ID_GOLD);

        hitCol += colOther*float(floor(d.y) == ID_OTHER);

        hitCol += colThird*float(floor(d.y) == ID_THIRD);

        
        
        //hitCol = mix(hitCol, hitCol + hitCol*ambientCol, shad);
    
    }
    
    if(hit)
        col += hitCol*atten;
    
    //col = mix(col, vec3(1),smoothstep(0.,1.,t*0.05 - 0.6));
    p = ro + rd*min(t,10.);
    
    float fogFact = smoothstep(0.,1.,i/maxIters*3.75 - 0.4);
    float fogFactB = smoothstep(0.,1.,length(t)*0.13 - .2);
    
    
    
    //fogFact *= 1. - cyclicNoise(p*2.,1.);
    
    
    vec3 fogCol = vec3(0.1,0.5,0.5)*0.02;
    
    //fogCol = mix(fogCol, fogCol*fogCol*fogCol*0.3, cyclicNoiseFog(normalize(rd)*14.,TIME*0.5)*pow(fogFact,4.2));
    fogCol = mix(fogCol, fogCol*fogCol*fogCol*0.3, cyclicNoiseFog(normalize(rd)*14. + vec3(0,0,smoothTimeB),smoothTimeB*0.5)*pow(fogFactB,4.2));
    
    
    
    col = mix(col, fogCol,fogFactB);
    
    
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


			//******** BuffB Code Begins ********


float FXAAamt = 0.5;




vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 p = fragCoord.xy/RENDERSIZE.xy;
    vec2 pp = 1. / RENDERSIZE.xy;
    vec4 color = texture(BuffA, vec2(fragCoord.xy * pp));
    vec3 luma = vec3(.299, 0.587, 0.114);
    
    float lumaNW = dot(texture(BuffA, (fragCoord.xy + vec2(-1.0, -1.0)*FXAAamt) * pp).xyz, luma);
    float lumaNE = dot(texture(BuffA, (fragCoord.xy + vec2(1.0, -1.0)*FXAAamt) * pp).xyz, luma);
    float lumaSW = dot(texture(BuffA, (fragCoord.xy + vec2(-1.0, 1.0)*FXAAamt) * pp).xyz, luma);
    float lumaSE = dot(texture(BuffA, (fragCoord.xy + vec2(1.0, 1.0)*FXAAamt) * pp).xyz, luma);
    float lumaM  = dot(color.xyz,  luma);
    float lumaMin = min(lumaM, min(min(lumaNW, lumaNE), min(lumaSW, lumaSE)));
    float lumaMax = max(lumaM, max(max(lumaNW, lumaNE), max(lumaSW, lumaSE)));

    vec2 dir = vec2(-((lumaNW + lumaNE) - (lumaSW + lumaSE)), ((lumaNW + lumaSW) - (lumaNE + lumaSE)));

    float dirReduce = max((lumaNW + lumaNE + lumaSW + lumaSE) *
                          (0.25 * (1.0/8.0)), (1.0/128.0));

    float rcpDirMin = 2.5 / (min(abs(dir.x), abs(dir.y)) + dirReduce);
    dir = min(vec2(8.0, 8.0),
              max(vec2(-8.0, -8.0),
              dir * rcpDirMin)) * pp;

    vec3 rgbA = 0.5 * (
        texture(BuffA, fragCoord.xy * pp + dir * (1.0 / 3.0 - 0.5)).xyz +
        texture(BuffA, fragCoord.xy * pp + dir * (2.0 / 3.0 - 0.5)).xyz);
    vec3 rgbB = rgbA * 0.5 + 0.25 * (
        texture(BuffA, fragCoord.xy * pp + dir * -0.5).xyz +
        texture(BuffA, fragCoord.xy * pp + dir * 0.5).xyz);

    float lumaB = dot(rgbB, luma);
    if ((lumaB < lumaMin) || (lumaB > lumaMax)){
        fragColor = vec4(rgbA, color.w);
    } else {
        fragColor = vec4(rgbB, color.w);
    }

	return fragColor; 
 } 


// wip


// shadows and smoothops from iq

// FXAA from mudlord?

// cyclic noise from nimitz

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv = (fragCoord - 0.5*RENDERSIZE.xy)/RENDERSIZE.y;

    vec3 col = vec3(0);

    col = texture(BuffB,fragCoord/RENDERSIZE.xy).xyz;
    
    
    col = pow(col,vec3(1.03,1.03,0.95 + dot(uv,uv)*0.02));
    
    
    col *= 3.3;
    col = mix(col,smoothstep(0.,1.,col*vec3(0.95,1.,1.)*1.),0.5);
    
    col = mix(acesFilm(col), col, 0.);
    
    
    col = mix(col,col*col*col*0.5,dot(uv,uv));
    //col *= 1. - dot(uv,uv*0.4)*2.;
    
    col = pow(col,vec3(0.454545));
    
    fragColor = vec4(col,1.0);
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
		return renderMainImage();
	}
}