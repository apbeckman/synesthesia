

// Infinite Cube Spiral by Anthony Hall
// Inspired by "Infinite Cube Zoom" by bitless:
// https://www.shadertoy.com/view/3lc3zH


// Transforms the distinct cubes into one continuous spiral.
// When disabled, it's pretty much a less fancy version of
// bitless' shader. I can't take any credit for the concept,
// I just really like how it looks without the spiral.
#define SPIRAL

// Supersamples pixels near the center. At high resolutions,
// this doesn't impact performance very much
#define SS

const float pi = radians(180.0);
const float twoPi = radians(360.0);

// Time in seconds to zoom through a full cube
const float zoomPeriod = 6.0;

// Globals set in mainImage
float time;
float base;

const vec3 toLight = normalize(vec3(3, 4, 2));

// Procedural color palette
vec3 palette(float t)
{
    t *= twoPi;
    
    vec2 colors = vec2(
        0.5 + 0.5 * cos(t - 1.0),
        0.6 - 0.4 * cos(t + 1.9));

    return colors.xxy;
}

// Returns the scene color
vec3 image(vec2 point)
{
    float r = length(point);
    float theta = atan(point.y, point.x);
    
    // Mapping which of the 6 slices we are in
    float region = floor((6.0 * theta + pi) / twoPi);
    
    // Set each slice to be repeating from -30 to 30 deg
    theta = (mod(6.0 * theta + pi, twoPi) - pi) / 6.0;
    
#ifdef SPIRAL

    // The spiral is achieved by offsetting each slice's radius
    // by 1/3 of a layer, then rotating each slice so the lines
    // connect to both neighbors. We can exactly solve for this
    // rotation using the law of cosines:
    // a^2 = b^2 + c^2 - 2bc cos A
    // Here, a is the ccw side and b is the clockwise side.
    
    float a = pow(base, -1.0 / 6.0);
    float b = 1.0 / a;
    float c = sqrt(a*a + b*b - a*b); // cos C = 0.5
    
    float A = acos((a*a - b*b - c*c) / (-2.0*b*c));
    float correction = pi/3.0 - A;
    theta -= correction;
    
    float level = log(r * cos(theta)) / log(base) - 2.0 * fract(smoothTime / zoomPeriod) + region / 3.0;

#else

    float level = log(r * cos(theta)) / log(base) - 2.0 * fract(smoothTime / zoomPeriod);

#endif

    // Calculate which direction this surface is facing for lighting.
    // The face index corresponds to the axis the normal is on
    float invert = floor(mod(level, 2.0));
    region = mod(region + 3.0 * invert - 2.0, 6.0);
    float face = floor(region / 2.0);

    float diffuse = toLight[int(face)];
    diffuse = 0.4 + 0.6 * diffuse;

    // Unfortunately, this palette only yields good colors for
    // 0 < t < 0.5
    vec3 color = palette((fract(level) + face) / 6.0);
    color *= diffuse;
    
    return color;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    time = TIME;

// I find that different base values are ideal for spiral vs non-spiral
#ifdef SPIRAL

    if (_mouse.z > 0.0)
        base = exp2(0.35 + 0.8 * _mouse.y / RENDERSIZE.y);
    else
        base = 1.52 + 0.08 * cos(0.18 * smoothTime);
    
#else

    if (_mouse.z > 0.0)
        base = exp2(0.4 + 1.6 * _mouse.y / RENDERSIZE.y);
    else
        base = 1.85 + 0.35 * cos(0.18 * smoothTime);

#endif

    vec3 color = vec3(0.0);
    
#ifdef SS
    
    // Since the level of detail is logarithmic, it is dependent
    // on the number of pixels away from the center, regardless
    // of resolution. Thus, supersampling is done based on absolute
    // pixel distance rather than percentage of the viewport.
    
    vec2 absPixelCoord = abs(fragCoord - 0.5 * RENDERSIZE.xy);
    float bound = max(absPixelCoord.x, absPixelCoord.y);
    
    // This samples 2x2 within 120 pixels of the center, and
    // 3x3 within 12 pixels. 120 is definitely overkill for
    // strictly avoiding aliasing, I just find the jaggies
    // visually unappealing when the threshold is smaller.
    
    int samples = 1 + int(bound < 120.0) + int(bound < 12.0);
    float increment = 1.0 / float(samples);
    float offset = 0.5 * increment - 0.5;
    
    for (int x = 0; x < samples; x++)
    {
        for (int y = 0; y < samples; y++)
        {
            vec2 ssFragCoord = fragCoord + vec2(x, y) * increment + offset;
            vec2 ssPoint = (2.0 * ssFragCoord - RENDERSIZE.xy) / RENDERSIZE.y;
            color += image(ssPoint);
        }
    }
    color /= float(samples * samples);

#else
    
    vec2 point = (2.0 * fragCoord - RENDERSIZE.xy) / RENDERSIZE.y; 
    color = image(point);
    
#endif

    fragColor = vec4(color, 1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}