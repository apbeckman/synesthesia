

			//******** Common Code Begins ********


// Dave Hoskins
// https://www.shadertoy.com/view/4djSRW
vec3 hash33(vec3 p3)
{
	p3 = fract(p3 * vec3(.1031, .1030, .0973));
    p3 += dot(p3, p3.yxz+33.33);
    return fract((p3.xxy + p3.yxx)*p3.zyx);
}

// Inigo Quilez
// https://iquilezles.org/articles/distfunctions
float smin( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); }

// rotation matrix
mat2 rot(float a) { return mat2(cos(a),-sin(a),sin(a),cos(a)); }



			//******** BuffA Code Begins ********

// taste of noise 11 by leon denise 2021/10/26
// thanks to Inigo Quilez, David Hoskins, NuSan, Fabrice Neyret and many others
// licensed under hippie love conspiracy

// global variable
vec3 rng;

// geometry
float map (vec3 p)
{
    float t = smoothTime * 0.0075;
    float t2 = smoothTimeC * 0.0075;
    
    // parameters
    vec3 angle = vec3(1,2,3)+rng.x*.1+sin(length(p))+p*.1;
    float size = (0.01*(1.0+basshits))*sin(rng.x*3.14)*rng.y;
    float range = .6;
    
    // geometric iteration
    const float count = 12.0;
    float a = 1.0;
    float scene = 1000.;
    for (float index = 0.0; index < count; ++index)
    {        
        // rotate
        p.yx *= rot((angle.z*(1.0-t))/a);
        p.xz *= rot((angle.y*(1.0-t2))/a);
        p.yz *= rot((angle.x*(1.0-t))/a);
        
        // fold
        vec3 ppp = p;
        p.x = abs(p.x)-range*a;
        
        // add sdf object
        scene = smin(scene, length(p)-size, 0.5*a);
        
        // falloff
        a /= 1.4;
    }
        
    return scene;
}

// Inigo Quilez (https://www.shadertoy.com/view/Xds3zN)
float getAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float h = 0.01 + 0.1*float(i)/4.0;
        float d = map( pos + h*nor );
        occ += (h-d)*sca;
        sca *= 0.95;
        if( occ>0.35 ) break;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 );// * (0.5+0.5*nor.y);
}

vec3 color (vec3 pos, vec3 ray, vec3 normal)
{    
    // lighting
    vec3 rf = reflect(ray, normal);
    float ld = dot(rf, vec3(0,0,1))*0.5+0.5;
    vec3 ld2 = vec3(0.875,0.722,1.000) * sqrt(ld);
    ld = dot(rf, normalize(vec3(0,1,0)))*0.5+0.5;
    vec3 light = vec3(0.580,0.918,1.000) * pow(ld,10.);

    // color palette by Inigo Quilez (https://iquilezles.org/articles/palettes)
    vec3 tint = .5+.5*cos(vec3(0, .3, .6)*4.+length(pos)*3.-1.);
    
    // ambient occlusion by Inigo Quilez (https://www.shadertoy.com/view/Xds3zN)
    float ao = mix(1., getAO(pos, normal), .9);

    // compositing
    return (tint + ld2 + light) * ao * .5;
}


vec4 renderPassA() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // reset color
    fragColor = vec4(0,0,0,1);
    
    // white noise
    vec3 seed = vec3(gl_FragCoord.xy, TIME);
    rng = hash33(seed);
    
    // pixel coordinates
    vec2 uv = (fragCoord.xy - RENDERSIZE.xy * 0.5) / RENDERSIZE.y;
    
    // blur edges
    vec3 rng3 = hash33(seed+78.);
    vec2 dof = vec2(cos(rng3.x*1.28),sin(rng3.x*1.28))*rng3.y;
    uv += dof*pow(length(uv), 8.0)*.5;
    
    // camera coordinates
    vec3 eye = vec3(0,0,-3);
    vec3 ray = normalize(vec3(uv, 2));
    vec3 pos = eye + ray * 1.;
    
    // raymarching
    const float count = 30.;
    bool hit = false;
    float index;
    for (index = 0.; index < count; ++index)
    {
        float dist = map(pos);
        if (dist < 0.001)
        {
            hit = true;
            break;
        }
        dist *= 0.9 + 0.1 * rng.z;
        pos += ray * dist;
    }
    
    if (hit)
    {
        // compute normal by NuSan (https://www.shadertoy.com/view/3sBGzV)
        vec2 off=vec2(0.001,0);
        vec3 normal = normalize(map(pos)-vec3(map(pos-off.xyy), map(pos-off.yxy), map(pos-off.yyx)));

        // coloring
        float shade = 1.-index/count;
        fragColor.rgb = color(pos, ray, normal) * shade;
    }
    
    // feedback with fade out
    vec4 frame = texture(BuffA, gl_FragCoord.xy/RENDERSIZE.xy);
    fragColor.rgb = max(fragColor.rgb, frame.rgb - 0.002);
	return fragColor; 
 } 








// taste of noise 11 by leon denise 2021/10/26
// thanks to Inigo Quilez, David Hoskins, NuSan, Fabrice Neyret and many others
// licensed under hippie love conspiracy

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    fragColor = texture(BuffA, gl_FragCoord.xy/RENDERSIZE.xy);
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