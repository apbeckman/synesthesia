//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 





			//******** BuffA Code Begins ********



// created by florian berger (flockaroo) - 2016

// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.



// single pass CFD

// ---------------

// this is some "computational flockarooid dynamics" ;)

// the self-advection is done purely rotational on all scales. 

// therefore i dont need any divergence-free velocity field. 

// with stochastic sampling i get the proper "mean values" of rotations 

// over time for higher order scales.

//

// try changing "RotNum" for different accuracies of rotation calculation

// for even RotNum uncomment the line #define SUPPORT_EVEN_ROTNUM



#define RotNum 7 

//#define SUPPORT_EVEN_ROTNUM



#define Res  RENDERSIZE

#define Res1 RENDERSIZE



#define keyTex iChannel3

#define KEY_I texture(keyTex,vec2((105.5-32.0)/256.0,(0.5+0.0)/3.0)).x

// implementation of MurmurHash (https://sites.google.com/site/murmurhash/) for a 

// single unsigned integer.



uint hash(uint x, uint seed) {

    const uint m = 0x5bd1e995U;

    uint hash = seed;

    // process input

    uint k = x;

    k *= m;

    k ^= k >> 24;

    k *= m;

    hash *= m;

    hash ^= k;

    // some final mixing

    hash ^= hash >> 13;

    hash *= m;

    hash ^= hash >> 15;

    return hash;

}



// implementation of MurmurHash (https://sites.google.com/site/murmurhash/) for a  

// 3-dimensional unsigned integer input vector.



uint hash(uvec3 x, uint seed){

    const uint m = 0x5bd1e995U;

    uint hash = seed;

    // process first vector element

    uint k = x.x; 

    k *= m;

    k ^= k >> 24;

    k *= m;

    hash *= m;

    hash ^= k;

    // process second vector element

    k = x.y; 

    k *= m;

    k ^= k >> 24;

    k *= m;

    hash *= m;

    hash ^= k;

    // process third vector element

    k = x.z; 

    k *= m;

    k ^= k >> 24;

    k *= m;

    hash *= m;

    hash ^= k;

	// some final mixing

    hash ^= hash >> 13;

    hash *= m;

    hash ^= hash >> 15;

    return hash;

}





vec3 gradientDirection(uint hash) {

    switch (int(hash) & 15) { // look at the last four bits to pick a gradient direction

    case 0:

        return vec3(1, 1, 0);

    case 1:

        return vec3(-1, 1, 0);

    case 2:

        return vec3(1, -1, 0);

    case 3:

        return vec3(-1, -1, 0);

    case 4:

        return vec3(1, 0, 1);

    case 5:

        return vec3(-1, 0, 1);

    case 6:

        return vec3(1, 0, -1);

    case 7:

        return vec3(-1, 0, -1);

    case 8:

        return vec3(0, 1, 1);

    case 9:

        return vec3(0, -1, 1);

    case 10:

        return vec3(0, 1, -1);

    case 11:

        return vec3(0, -1, -1);

    case 12:

        return vec3(1, 1, 0);

    case 13:

        return vec3(-1, 1, 0);

    case 14:

        return vec3(0, -1, 1);

    case 15:

        return vec3(0, -1, -1);

    }

}



float interpolate(float value1, float value2, float value3, float value4, float value5, float value6, float value7, float value8, vec3 t) {

    return mix(

        mix(mix(value1, value2, t.x), mix(value3, value4, t.x), t.y),

        mix(mix(value5, value6, t.x), mix(value7, value8, t.x), t.y),

        t.z

    );

}



vec3 fade(vec3 t) {

    // 6t^5 - 15t^4 + 10t^3

	return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);

}



float _perlinNoise(vec3 position, uint seed) {

    vec3 floorPosition = floor(position);

    vec3 fractPosition = position - floorPosition;

    uvec3 cellCoordinates = uvec3(floorPosition);

    float value1 = dot(gradientDirection(hash(cellCoordinates, seed)), fractPosition);

    float value2 = dot(gradientDirection(hash((cellCoordinates + uvec3(1, 0, 0)), seed)), fractPosition - vec3(1, 0, 0));

    float value3 = dot(gradientDirection(hash((cellCoordinates + uvec3(0, 1, 0)), seed)), fractPosition - vec3(0, 1, 0));

    float value4 = dot(gradientDirection(hash((cellCoordinates + uvec3(1, 1, 0)), seed)), fractPosition - vec3(1, 1, 0));

    float value5 = dot(gradientDirection(hash((cellCoordinates + uvec3(0, 0, 1)), seed)), fractPosition - vec3(0, 0, 1));

    float value6 = dot(gradientDirection(hash((cellCoordinates + uvec3(1, 0, 1)), seed)), fractPosition - vec3(1, 0, 1));

    float value7 = dot(gradientDirection(hash((cellCoordinates + uvec3(0, 1, 1)), seed)), fractPosition - vec3(0, 1, 1));

    float value8 = dot(gradientDirection(hash((cellCoordinates + uvec3(1, 1, 1)), seed)), fractPosition - vec3(1, 1, 1));

    return interpolate(value1, value2, value3, value4, value5, value6, value7, value8, fade(fractPosition));

}



float perlinNoise(vec3 position, int frequency, int octaveCount, float persistence, float lacunarity, uint seed) {

    float value = 0.0;

    float amplitude = 1.0;

    float currentFrequency = float(frequency);

    uint currentSeed = seed;

    for (int i = 0; i < octaveCount; i++) {

        currentSeed = hash(currentSeed, 0x0U); // create a new seed for each octave

        value += _perlinNoise(position * currentFrequency, currentSeed) * amplitude;

        amplitude *= persistence;

        currentFrequency *= lacunarity;

    }

    return value;

}

float Noise = perlinNoise(vec3(_xy/RENDERSIZE, TIME * 0.25), 1, 6, 0.5, 2.0, uint(6));

float ang = 2.0*3.141592653589793238462643383279502884197/float(RotNum);

mat2 m = mat2(cos(ang),sin(ang),-sin(ang),cos(ang));

mat2 mh = mat2(cos(ang*0.5),sin(ang*0.5),-sin(ang*0.5),cos(ang*0.5));



vec4 randS(vec2 uv)

{

    return texture(image30,uv*Res.xy/Res1.xy)-vec4(0.25);

}



float getRot(vec2 pos, vec2 b)

{

    vec2 p = b;

    float rot=0.0;

    for(int i=0;i<RotNum;i++)

    {

        rot+=dot(texture(BuffA,fract((pos+p)/Res.xy)).xy-vec2(0.5),p.yx*vec2(1+_noise(sin(TIME)),-1+_noise(cos(TIME))));
        //rot+=dot(texture(BuffA,fract((pos+p)/Res.xy)).xy-vec2(0.5),p.yx*vec2(-2,2));

        p = m*p;

    }

    return rot/float(RotNum)/dot(b,b);

}



vec4 renderPassA() {

	vec4 fragColor = vec4(0.0);

	vec2 fragCoord = _xy;

    fragCoord -= _uvc;

    vec2 pos = fragCoord.xy;
    
    float rnd = randS(vec2(float(FRAMECOUNT)/Res.x,0.5/Res1.y)).x;

    

    vec2 b = vec2(cos(ang*rnd),sin(ang*rnd));

    vec2 v=vec2(0);

    float bbMax=0.7*Res.y; bbMax*=bbMax;

    for(int l=0;l<20;l++)

    {

        if ( dot(b,b) > bbMax ) break;

        vec2 p = b;

        for(int i=0;i<RotNum;i++)

        {

#ifdef SUPPORT_EVEN_ROTNUM

            v+=p.yx*getRot(pos+p,-mh*b);

#else

            // this is faster but works only for odd RotNum

            v+=p.yx*getRot(pos+p,b);

#endif

            p = m*p;

        }

        b*=2.0;

    }

    

    fragColor=texture(BuffA,fract((pos+v*vec2(-1,1)*(.1250+Impact*(syn_MidLevel*0.125+0.25* syn_BassLevel)))/(Res.xy+_uvc*Noise)));

    // /fragColor.rgb += 0.01*(_loadUserImageAsMask().rgb);



//    /fragColor *=1.0+0.00125* texture(syn_UserImage,fragCoord.xy/(Res.xy))*syn_Intensity;

    // add a little "motor" in the center

    vec2 scr=(fragCoord.xy/Res.xy)*2.0-vec2(1.0);

    fragColor.xy += ((0.001)*scr.xy / (dot(scr,scr)/0.1+0.3));

    if(FRAMECOUNT<=4 || Reset>0.5) 

          if (_exists(syn_UserImage)){

            

        fragColor=texture(syn_UserImage,fragCoord.xy/(Res.xy));

          }

    else{

    fragColor=texture(image5,fragCoord.xy/(Res.xy));

    }
    fragColor = mix(fragColor, _edgeDetectSobel(syn_UserImage)*0.2, 0+syn_Hits*0.0025);

    return fragColor; 

 } 







// created by florian berger (flockaroo) - 2016

// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.



// single pass CFD

// ---------------

// this is some "computational flockarooid dynamics" ;)

// the self-advection is done purely rotational on all scales. 

// therefore i dont need any divergence-free velocity field. 

// with stochastic sampling i get the proper "mean values" of rotations 

// over time for higher order scales.

//

// try changing "RotNum" for different accuracies of rotation calculation

// for even RotNum uncomment the line #define SUPPORT_EVEN_ROTNUM



float getVal(vec2 uv)

{

    return length(texture(BuffA,uv).xyz);

}

    

vec2 getGrad(vec2 uv,float delta)

{

    vec2 d=vec2(delta,0);

    return vec2(

        getVal(uv+d.xy)-getVal(uv-d.xy),

        getVal(uv+d.yx)-getVal(uv-d.yx)

    )/delta;

}



vec4 renderMainImage() {

	vec4 fragColor = vec4(0.0);

	vec2 fragCoord = _xy;



	vec2 uv = fragCoord.xy / RENDERSIZE.xy;

    vec3 n = vec3(getGrad(uv,1.0/RENDERSIZE.y),250.0);

    //n *= n;

    n=normalize(n);

    fragColor=vec4(n,1);

    vec3 light = normalize(vec3(1,1,2));

    float diff=clamp(dot(n,light),0.5,1.0)*(1.0+0.125*syn_HighHits);

    float spec=clamp(dot(reflect(light,n),vec3(0,0,-1)),0.0,1.0);

    spec=pow(spec,128.0)*2.5;

    //spec=0.0;

	fragColor = texture(BuffA,uv)*vec4(diff)+vec4(spec);
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