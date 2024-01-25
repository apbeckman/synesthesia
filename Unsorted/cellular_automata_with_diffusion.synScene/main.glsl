

			//******** BuffA Code Begins ********

bool fc = FRAMECOUNT <= 1 || Reset != 0.0;
float growthFactor = pow((syn_BassLevel*0.35)+(syn_MidLevel*0.35)+syn_Intensity*0.3, 2.0)*(0.5+0.35*syn_Intensity);
// vec4 mixedEdgeCol = mix( _loadMedia(), mediaEdges, .1+0.5*edge_mix)*media_impact*0.25;
float mediahits_mix =  mix(media_impact, media_impact*(.125+0.875*syn_BassLevel), media_hits);

vec4 mediaEdges = texture(media_pass_fx, _uv);
bool mediaOn = _exists(syn_UserImage);
float media_lum = sin(PI*0.75*length(mediaEdges))*mediahits_mix;
// float paintsize = _smin(70.*paint_size, 60.*paint_size, 1.0);
float paintsize = 60.*paint_size;



// All components are in the range [0…1], including hue.
vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

vec3 rgb2hsv(vec3 c)
{
    vec4 K = vec4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
    vec4 p = mix(vec4(c.bg, K.wz), vec4(c.gb, K.xy), step(c.b, c.g));
    vec4 q = mix(vec4(p.xyw, c.r), vec4(c.r, p.yzx), step(p.x, c.r));

    float d = q.x - min(q.w, q.y);
    float e = 1.0e-10;
    return vec3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

int xorshift(in int value) // For random number generation
{
    // Xorshift*32
    // Based on George Marsaglia's work: http://www.jstatsoft.org/v08/i14/paper
    value ^= value << 13;
    value ^= value >> 17;
    value ^= value << 5;
    return value;
}

int nextInt(inout int seed)  // RNG
{
    seed = xorshift(seed);
    return seed;
}

float nextFloat(inout int seed) // RNG
{
    seed = xorshift(seed);
    return abs(fract(float(seed) / 31416.592898653));
}

float nextFloat(inout int seed, in float max) // RNG
{
    return nextFloat(seed) * max;
}

vec3 nearby(in vec2 fragCoord, in vec2 offset) // For getting pixels from prev. frame
{
    vec2 uv = (fragCoord + offset)/RENDERSIZE.xy;
    return vec3(texture(BuffA, uv));
}


float hueDiff(in float a, in float b) // Finds the shortest difference between two hues
{	
    float diff=fract(a)-fract(b);
    if(abs(diff)<0.5)
        return diff;
    else
        return diff - 1.*sign(diff);
}

float checkfunction(in float near_hue, in float prev_hue) // used to determine the likelyhood of a pixel moving
{
    return (sign(hueDiff(near_hue,prev_hue)) + 4.*hueDiff(near_hue,prev_hue))/abs(hueDiff(near_hue,prev_hue));
}



const float threshold = 0.1; // Amount of randomness



vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    int rngSeed = int(fragCoord.x) + int(fragCoord.y) * int(RENDERSIZE.x) + int(TIME * 1000.0);
	
    float mouseDown = _mouse.z;
    vec2 uvMouse = vec2(_mouse)/RENDERSIZE.xy;

    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/RENDERSIZE.xy;
    vec3 previous = rgb2hsv(vec3(texture(BuffA, uv)));
    vec3 next = previous;
    
    
    if(next[2]>0.05)    // Slowly fade un-updated pixels to black 
    	next[2]*=0.99;
    else
    {
    	if(nextFloat(rngSeed, 1.0)>0.9) // Seed some randomness to stop it from dying
        {
         	next[2]=1.;
            next[0]=nextFloat(rngSeed, 1.0);
            next[1]=1.;
                
        }
    }
    
    
    if(nextFloat(rngSeed, 1.0) > threshold) // randomly update pixels
    {
        
        // Adjacent pixles
        vec3 up = rgb2hsv(nearby(fragCoord, vec2(0.0, -1.0)));
        vec3 down = rgb2hsv(nearby(fragCoord, vec2(0.0, 1.0)));
        vec3 left = rgb2hsv(nearby(fragCoord, vec2(-1.0, 0.0)));
        vec3 right = rgb2hsv(nearby(fragCoord, vec2(1.0, 0.0)));
        
        // NOTE: everything is now [hue, sat, val] rather than [red, green, blue]
        
        // Weights of this pixel becoming the same color as one of the adjacents
        float upweight = up[1]*(abs(hueDiff(up[0],previous[0]))-0.3)*(0.7-abs(hueDiff(up[0],previous[0])))*(0.1+0.6*up[2])*checkfunction(up[0],previous[0]);
        float downweight = down[1]*(abs(hueDiff(down[0],previous[0]))-0.3)*(0.7-abs(hueDiff(down[0],previous[0])))*(0.1+0.6*down[2])*checkfunction(down[0],previous[0]);
        float leftweight = left[1]*(abs(hueDiff(left[0],previous[0]))-0.3)*(0.7-abs(hueDiff(left[0],previous[0])))*(0.1+0.6*left[2])*checkfunction(left[0],previous[0]);
        float rightweight = right[1]*(abs(hueDiff(right[0],previous[0]))-0.3)*(0.7-abs(hueDiff(right[0],previous[0])))*(0.1+0.6*right[2])*checkfunction(right[0],previous[0]);
        float _max = max(max(upweight,downweight),max(leftweight,rightweight));
		
        // Only update if the best option is strong enough
        if(_max>((previous[2]))*1.7)
        {
            if( _max == upweight)
                next = up;
            
            if (_max == downweight)
                next = down;

            if (_max == leftweight)
                next = left;

            if (_max == rightweight)
                next = right;
            // Inject some saturation and lightness to new pixles
            if(next[2]<1.)
                next[2]+=0.5+0.5*(abs(hueDiff(next[0],previous[0]))-1.);
            if(next[1]<1.)
                next[1]+=0.0003; 
        }
        else // make pixels diffuse if they can't update. Stops single-pixel islands
        {
            vec3 up_rgb = (nearby(fragCoord, vec2(0.0, -1.0)));
            vec3 down_rgb = (nearby(fragCoord, vec2(0.0, 1.0)));
            vec3 left_rgb = (nearby(fragCoord, vec2(-1.0, 0.0)));
            vec3 right_rgb = (nearby(fragCoord, vec2(1.0, 0.0)));
            vec3 prev_rgb = vec3(texture(BuffA, uv));
            next = rgb2hsv(prev_rgb*0.9 + 0.1*(up_rgb+down_rgb+left_rgb+right_rgb)/4.);
        }
    }

	// use the webcam for the first frame
    if((FRAMECOUNT < 5) || mouseDown>0.5) 
    {
        next = rgb2hsv(vec3(texture(media_pass_fx,uv)));
        previous = next;
        
    }
    

    // Output to screen
    fragColor = vec4(hsv2rgb(next),1.0);

	return fragColor; 
 } 





vec4 mediaPass() {
    vec4 media = vec4(0.);
    vec2 U = _xy;
    media = _loadMedia();
    // media = _loadMedia()*media_impact*2.;
    return media;
 } 
vec4 mediaPassFX() {
    vec4 media = vec4(0.);
    vec4 media_edge_mix = mix(_edgeDetectSobel(media_pass), texture(media_pass, _uv), 1.0-edge_mix*0.975-0.0125)*media_impact;
    return media_edge_mix;
 } 
vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // Output to screen
    vec2 uv = fragCoord/RENDERSIZE.xy;
    fragColor = texture(BuffA, uv);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderPassA();
	}
	if(PASSINDEX == 1){
		return mediaPass();
	}
	if(PASSINDEX == 2){
		return mediaPassFX();
	}
	if(PASSINDEX == 3){
		return renderMainImage();
	}
}