

			//******** Common Code Begins ********


// #define USE_CINEMATIC_MODE      // Uncomment this line for a more cinematic view (camera will side-scroll)
// #define USE_BUMPY_STREAMS_MODE  // Uncomment this line to make streams bumpy (sausages-like)
// #define USE_GENERATION_SEED 123 // Uncomment this line to use a fixed generation seed (then reset the simulation to apply the changes)
    
const int nParticles = 20;
const float particlesSize = 8.0;
const float collisionDamping = 0.5;
const float streamsFadingExp = 0.001;
 float gravityStrength = (1.6+syn_BassLevel) / particlesSize;

const vec3 ambientLightDir = normalize(vec3(1.0, 2.0, 0.0));
const vec3 ambientLightCol = vec3(1.1, 1.0, 0.9);
const vec3 backgroundColor = vec3(0.01);
const float streamsGlossExp = 20.0;
const float spotlightsGlare = 0.0;

#ifdef USE_BUMPY_STREAMS_MODE
#define particlesSize mix(particlesSize, particlesSize * 0.5, (1.0 + sin(1.85 + TIME * 11.93805208)) * 0.5)
#endif

#ifdef USE_GENERATION_SEED
#define generationSeed float(USE_GENERATION_SEED) // a fixed seed will generate the same output (in respect of the viewport size)
#else
#define generationSeed vec4(2019.0, 1.0, 1.0, TIME).w // if no custom seed is provided, POSIX time is used instead (producing different results every time)
#endif

const ivec2 cameraVelocity =
#ifdef USE_CINEMATIC_MODE
ivec2(1, 0);
#else
ivec2(0);
#endif

// Buf A: particles positions and inertia
// Buf B: scene albedo  (accumulated)
// Buf C: scene normals (accumulated)
// Image: final compositing


			//******** BuffA Code Begins ********

// Compute Physics (Verlet Integration)

float rand(in vec2 co) {
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

vec2 randVec2(in vec2 co) {
	return vec2(rand(co.xy + generationSeed * 0.0001), rand(-co.yx + generationSeed * 0.0001));
}

vec2 randNrm2(in vec2 fragCoord)
{
	vec2 n = vec2(-1.0) + randVec2(fragCoord) * 2.0;
    
    float l = length(n);   
    if(l <= 0.0000001) n = vec2(0.0, (l = 1.0));
    
    return (n / l);
}

void initParticle(in vec2 fragCoord, inout vec2 particlePrevPosition, inout vec2 particleCurrPosition)
{
	particleCurrPosition = randVec2(fragCoord) * RENDERSIZE.xy;
    particlePrevPosition = particleCurrPosition - randNrm2(fragCoord) * particlesSize * 0.0625;
}

vec2 getParticlePosition(in int particleID)
{
    int BuffA_width = int(RENDERSIZE.x);
	ivec2 particleCoord = ivec2(particleID % BuffA_width, particleID / BuffA_width);
    
    return texelFetch(BuffA, particleCoord, 0).xy;
}

vec2 computeGravitation(in int particleID, in vec2 particlePosition)
{
    vec2 acceleration = vec2(0.0);
        
	for(int i = 0; i < nParticles; ++i) if(i != particleID)
    {
        vec2 v = (getParticlePosition(i) - particlePosition);
        float d = length(v);
        
        if(d > 0.0000001) acceleration += (v / d) / pow(max(d, particlesSize * 2.0) * gravityStrength, 2.0);
    }
    
    return acceleration;
}

void solveCollisions(inout vec2 particlePrevPosition, inout vec2 particleCurrPosition)
{
    vec2 particleInertia = (particleCurrPosition - particlePrevPosition);
    
	if(particleCurrPosition.x < particlesSize || particleCurrPosition.x > RENDERSIZE.x - particlesSize)
    {
    	particleCurrPosition.x = clamp(particleCurrPosition.x, particlesSize, RENDERSIZE.x - particlesSize);
        particlePrevPosition.x = particleCurrPosition.x + particleInertia.x * collisionDamping;
    }
    
    if(particleCurrPosition.y < particlesSize || particleCurrPosition.y > RENDERSIZE.y - particlesSize)
    {
    	particleCurrPosition.y = clamp(particleCurrPosition.y, particlesSize, RENDERSIZE.y - particlesSize);
        particlePrevPosition.y = particleCurrPosition.y + particleInertia.y * collisionDamping;
    }
}

vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    int particleID = int(floor(fragCoord.x) + RENDERSIZE.x * floor(fragCoord.y));
    if(particleID >= nParticles) return fragColor;
    
    vec4 particleData = texelFetch(BuffA, ivec2(fragCoord), 0);
    vec2 particlePrevPosition = particleData.zw;
    vec2 particleCurrPosition = particleData.xy;
     
    if(FRAMECOUNT <= 1) initParticle(fragCoord, particlePrevPosition, particleCurrPosition);
   
    vec2 particleAcceleration = computeGravitation(particleID, particleCurrPosition);
    vec2 particleInertia = particleCurrPosition - particlePrevPosition;
    vec2 particleVelocity = particleInertia + particleAcceleration;
    
    particlePrevPosition = particleCurrPosition;
    particleCurrPosition += particleVelocity;
    
    solveCollisions(particlePrevPosition, particleCurrPosition);
    
    fragColor = vec4(particleCurrPosition, particlePrevPosition);
	return fragColor; 
 } 



			//******** BuffB Code Begins ********

// Compute Scene Albedo
/*
vec2 getParticlePosition(in int particleID)
{
    int BuffA_width = int(RENDERSIZE.x);
	ivec2 particleCoord = ivec2(particleID % BuffA_width, particleID / BuffA_width);
    
    return texelFetch(BuffA, particleCoord, 0).xy;
}
*/
vec3 getParticleColor(in vec2 p) {
    return normalize(vec3(0.1) + texture(image16, p * 0.1 + TIME * 0.005).rgb);
}

vec4 renderPassB() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
 
    fragColor = texelFetch(BuffB, ivec2(fragCoord) + cameraVelocity, 0);
    fragColor.a *= (1.0 - streamsFadingExp);
        
	for(int i = 0; i < nParticles; ++i)
    {
        vec2 particlePos = getParticlePosition(i);
        vec3 particleCol = getParticleColor(particlePos);
        
        float alpha = smoothstep(particlesSize, particlesSize * 0.5, distance(fragCoord, particlePos));
        fragColor = mix(fragColor, vec4(particleCol , 1.0), alpha);
    }
	return fragColor; 
 } 



			//******** BuffC Code Begins ********

// Compute Scene Normals
/*
vec2 getParticlePosition(in int particleID)
{
    int BuffA_width = int(RENDERSIZE.x);
	ivec2 particleCoord = ivec2(particleID % BuffA_width, particleID / BuffA_width);
    
    return texelFetch(BuffA, particleCoord, 0).xy;
}
*/
vec4 renderPassC() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;
 
    fragColor = texelFetch(BuffC, ivec2(fragCoord) + cameraVelocity, 0);
    
	for(int i = 0; i < nParticles; ++i)
    {
        vec2 v = fragCoord - getParticlePosition(i);
        
        float l = length(v);
        float alpha = smoothstep(particlesSize, particlesSize * 0.5, l);
        
        float z = sqrt(abs(particlesSize * particlesSize - l * l));
        vec3 n = normalize(vec3(v, z));

        fragColor = mix(fragColor, vec4(n, 1.0), alpha);
    }
	return fragColor; 
 } 



// Final Compositing (Deferred Lighting + Bloom)
/*
vec2 getParticlePosition(in int particleID)
{
    int BuffA_width = int(RENDERSIZE.x);
	ivec2 particleCoord = ivec2(particleID % BuffA_width, particleID / BuffA_width);
    
    return texelFetch(BuffA, particleCoord, 0).xy;
}
*/
vec3 computeLighting( in vec3 surfaceAlbedo,
                      in vec3 surfaceNormal,
                      in float surfaceGloss,
                      in vec3 lightCol,
                      in vec3 lightDir,
                      in float lightSpec,
                      in float lightAmb )
{
    float dot_n  = clamp(dot(surfaceNormal, lightDir), 0.0, 1.0);
    
    vec3 diffuse  = lightCol * surfaceAlbedo * clamp(dot_n, lightAmb, 1.0);
    vec3 specular = lightCol * float(dot_n > 0.0) * pow(clamp(dot(reflect(-lightDir, surfaceNormal), vec3(0.0, 0.0, 1.0)), 0.0, 1.0), surfaceGloss);
    
    return diffuse + specular * lightSpec;
}

vec3 computeSpotLight( in vec3 surfaceAlbedo,
                       in vec3 surfaceNormal,
                       in float surfaceGloss,
                       in vec3 surfacePos,  
                       in vec3 lightCol,
                       in vec3 lightPos,
                       in float lightRadius )
{
    vec3 lightVec = lightPos - surfacePos;
    float contribution = 1.0 / max(dot(lightVec, lightVec) * 0.08 / (lightRadius * lightRadius), 1.0);
    
    return computeLighting(surfaceAlbedo, surfaceNormal, surfaceGloss, lightCol, normalize(lightVec), 0.066667 * surfaceGloss, 0.0) * contribution;
}

vec3 computeLightGlow(in vec3 position, in vec3 lightCol, in vec3 lightPos, in float lightRadius)
{
    vec3 glare = spotlightsGlare * lightCol * smoothstep(lightRadius * 10.0, 0.0, length((lightPos.xy - position.xy) * vec2(1.0, 16.0)));
    vec3 innerGlow = vec3(0.8) * smoothstep(lightRadius, lightRadius * 0.5, distance(lightPos.xy, position.xy));
    vec3 outerGlow = 0.25 * lightCol * smoothstep(lightRadius * 2.5, 0.0, distance(lightPos.xy, position.xy));
  
    return innerGlow + outerGlow + glare;
}

vec3 computeVignetting(in vec2 fragCoord, in vec3 src) // https://www.shadertoy.com/view/4lSXDm
{
	vec2 coord = ((fragCoord.xy / RENDERSIZE.xy) - 0.5) * (RENDERSIZE.x / RENDERSIZE.y) * 2.0;
    float rf = sqrt(dot(coord, coord)) * 0.25;
    float rf2_1 = rf * rf + 1.0;
    
	return src * pow((1.0 / (rf2_1 * rf2_1)), 2.24);
}    

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec4 albedo = texelFetch(BuffB, ivec2(fragCoord), 0);
	vec3 normal = normalize(texelFetch(BuffC, ivec2(fragCoord), 0).xyz);
    vec3 position = vec3(fragCoord, -(1.0 - albedo.a) * 384.0 / particlesSize); // fake Z-depth from fade level
        
    fragColor = vec4(vec3(0.0), albedo.a); 
    fragColor.rgb += computeLighting(albedo.rgb, normal, streamsGlossExp, ambientLightCol, ambientLightDir, 0.5, 0.175);
    
    for(int i = 0; i < nParticles; ++i)
    {
        vec3 particlePos = vec3(getParticlePosition(i), 0.0);
        vec3 particleCol = texelFetch(BuffB, ivec2(particlePos.xy), 0).rgb;
            
        fragColor.rgb += computeSpotLight(albedo.rgb, normal, streamsGlossExp, position, particleCol, particlePos, particlesSize);
    }
    
    fragColor.rgb = 1.25 * fragColor.rgb - vec3(0.075);
    fragColor.rgb = mix(backgroundColor, fragColor.rgb, min(fragColor.a * .25, 1.0));
    fragColor.rgb = computeVignetting(fragCoord, fragColor.rgb);
    
    for(int i = 0; i < nParticles; ++i)
    {
        vec3 particlePos = vec3(getParticlePosition(i), 0.0);
        vec3 particleCol = texelFetch(BuffB, ivec2(particlePos.xy), 0).rgb;
        
        fragColor.rgb += computeLightGlow(position, particleCol, particlePos, particlesSize);
    }
    
    fragColor = vec4(pow(fragColor.rgb, vec3(1.0 / 2.24)), 1.0); // gamma correction
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
		return renderMainImage();
	}
}